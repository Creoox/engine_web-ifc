/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

#pragma once

#if defined(DEBUG_DUMP_SVG) || defined(DUMP_CSG_MESHES)
#include "../../test/io_helpers.h"
#include "../../test/dumpToThree.h"
#endif

#include "../nurbs.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <CDT.h>
#include <tinynurbs/tinynurbs.h>
#include <spdlog/spdlog.h>
#include "geometryutils.h"
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace webifc::geometry
{
    constexpr double DEG_FULL = 360.0;
    constexpr double DEG_HALF = 180.0;
    constexpr double EPS_ANGLE = 1e-9;

    inline double VectorToAngle3D(double x, double y)
    {
        if (std::abs(x) < EPS_MINISCULE && std::abs(y) < EPS_MINISCULE)
            return 0.0;

        double angleRad = std::atan2(x, y);
        return glm::degrees(angleRad);
    }

    inline void NormalizeAngle(double& angle)
    {
        angle = std::fmod(angle, DEG_FULL);
        if (angle < 0.0)
            angle += DEG_FULL;
    }

    inline void UnwrapAngles(std::vector<double>& angles)
    {
        if (angles.size() < 2)
            return;

        for (size_t i = 1; i < angles.size(); ++i)
        {
            double delta = angles[i] - angles[i - 1];

            if (delta > DEG_HALF)
                angles[i] -= DEG_FULL;
            else if (delta < -DEG_HALF)
                angles[i] += DEG_FULL;
        }
    }

    inline void AppendGeometry(IfcGeometry& geometry, const bimGeometry::Geometry& geom)
    {
        for (int i = 0; i < geom.numFaces; ++i)
        {
            auto f = geom.GetFace(i);
            geometry.AddFace( geom.GetPoint(f.i0), geom.GetPoint(f.i1), geom.GetPoint(f.i2));
        }
    }

    inline std::unordered_map<int, std::vector<glm::dvec3>> GroupBoundsByIndex(const std::vector<IfcBound3D>& bounds)
    {
        std::unordered_map<int, std::vector<glm::dvec3>> groups;
        for (const auto& bound : bounds)
        {
            for (size_t i = 0; i < bound.curve.points.size(); ++i)
            {
                groups[bound.curve.indices[i]].push_back( bound.curve.points[i]);
            }
        }

        return groups;
    }

    inline std::pair<double, double> ComputeAngleInterval(const std::vector<glm::dvec3>& points, const glm::dvec3& origin,
            const glm::dvec3& axisX, const glm::dvec3& axisY)
    {
        std::vector<double> angles;
        angles.reserve(points.size());

        for (const auto& p : points)
        {
            glm::dvec3 v = p - origin;
            double dx = glm::dot(axisX, v);
            double dy = glm::dot(axisY, v);

            double a = VectorToAngle3D(dx, dy);
            NormalizeAngle(a);
            angles.push_back(a);
        }

        if (angles.empty())
            return { 0.0, 0.0 };

        UnwrapAngles(angles);

        auto [minIt, maxIt] = std::minmax_element(
            angles.begin(), angles.end());

        return { *minIt, *maxIt };
    }

    inline void TriangulateRevolution( IfcGeometry& geometry, const std::vector<IfcBound3D>& bounds, const IfcSurface& surface, double circleSegments)
    {
        spdlog::debug("[TriangulateRevolution()]");

#if defined(DEBUG_DUMP_SVG) || defined(DUMP_CSG_MESHES)
        for (size_t ii = 0; ii < bounds.size(); ++ii)
        {
            const auto& bound = bounds[ii];
            if (bound.curve.points.empty()) continue;
            std::string fileName = "faceBounds" + std::to_string(ii) + ".html";
            webifc::dump::DumpCurveToHtml(bound.curve.points, fileName);
        }
#endif

        // Extract revolution axis coordinate system
        glm::dvec3 center = glm::dvec3(surface.RevolutionSurface.Direction[3]);
        glm::dvec3 vecX = glm::normalize(glm::dvec3(surface.RevolutionSurface.Direction[0]));
        glm::dvec3 vecY = glm::normalize(glm::dvec3(surface.RevolutionSurface.Direction[1]));
        glm::dvec3 vecZ = glm::normalize(glm::dvec3(surface.RevolutionSurface.Direction[2]));

        // ── Build profile in meridional (radius, height) space ───────────
        const auto& profilePts = surface.RevolutionSurface.Profile.curve.points;
        if (profilePts.size() < 2) return;

        struct MeridPt { double r; double h; };
        std::vector<MeridPt> merid;
        std::vector<double> arcLen;  // cumulative arc length along profile
        merid.reserve(profilePts.size());
        arcLen.reserve(profilePts.size());

        for (const auto& p : profilePts)
        {
            glm::dvec3 rel = p - center;
            double dx = glm::dot(vecX, rel);
            double dy = glm::dot(vecY, rel);
            double dz = glm::dot(vecZ, rel);
            merid.push_back({std::sqrt(dx * dx + dy * dy), dz});
        }
        arcLen.push_back(0.0);
        for (size_t i = 1; i < merid.size(); i++)
        {
            double dr = merid[i].r - merid[i-1].r;
            double dh = merid[i].h - merid[i-1].h;
            arcLen.push_back(arcLen.back() + std::sqrt(dr*dr + dh*dh));
        }
        double totalArc = arcLen.back();
        if (totalArc < 1e-12) return;

        // Interpolate profile at arc-length s → (radius, height)
        auto profileAt = [&](double s) -> MeridPt
        {
            if (s <= arcLen.front()) return merid.front();
            if (s >= arcLen.back())  return merid.back();
            auto it = std::lower_bound(arcLen.begin(), arcLen.end(), s);
            size_t idx = static_cast<size_t>(std::distance(arcLen.begin(), it));
            if (idx == 0) idx = 1;
            double t = (s - arcLen[idx-1]) / (arcLen[idx] - arcLen[idx-1]);
            return {
                merid[idx-1].r + t * (merid[idx].r - merid[idx-1].r),
                merid[idx-1].h + t * (merid[idx].h - merid[idx-1].h)
            };
        };

        // UV → 3D: u = angle (degrees), v = arc-length along profile
        auto uvTo3D = [&](double u, double v) -> glm::dvec3
        {
            auto mp = profileAt(v);
            double uRad = glm::radians(u);
            return center + mp.h * vecZ
                 + std::sin(uRad) * mp.r * vecX
                 + std::cos(uRad) * mp.r * vecY;
        };

        // ── Step 1: Project bounds into UV space as separate closed loops ────
        //   u = revolution angle (degrees)
        //   v = arc-length along profile in meridional plane
        struct UVLoop { std::vector<glm::dvec2> points; };
        std::vector<UVLoop> uvLoops;

        for (const auto& bound : bounds)
        {
            if (bound.curve.points.size() < 3) continue;
            std::vector<double> angles, vParams;
            for (const auto& p : bound.curve.points)
            {
                glm::dvec3 rel = p - center;
                double dx = glm::dot(vecX, rel);
                double dy = glm::dot(vecY, rel);
                double dz = glm::dot(vecZ, rel);
                double dd = std::sqrt(dx * dx + dy * dy);

                double u = VectorToAngle3D(dx, dy);
                NormalizeAngle(u);
                angles.push_back(u);

                // Find closest point on profile in meridional plane
                double bestS = 0.0;
                double bestDistSq = 1e30;
                for (size_t i = 0; i + 1 < merid.size(); i++)
                {
                    double r0 = merid[i].r, h0 = merid[i].h;
                    double sr = merid[i+1].r - r0, sh = merid[i+1].h - h0;
                    double segLenSq = sr*sr + sh*sh;
                    double t = 0.0;
                    if (segLenSq > 1e-24)
                    {
                        t = ((dd - r0)*sr + (dz - h0)*sh) / segLenSq;
                        t = std::max(0.0, std::min(1.0, t));
                    }
                    double pr = r0 + t * sr, ph = h0 + t * sh;
                    double distSq = (dd - pr)*(dd - pr) + (dz - ph)*(dz - ph);
                    if (distSq < bestDistSq)
                    {
                        bestDistSq = distSq;
                        bestS = arcLen[i] + t * (arcLen[i+1] - arcLen[i]);
                    }
                }
                vParams.push_back(bestS);
            }
            UnwrapAngles(angles);

            UVLoop loop;
            for (size_t i = 0; i < angles.size(); i++)
                loop.points.push_back({angles[i], vParams[i]});

            // Remove duplicate closing vertex(es)
            while (loop.points.size() >= 2)
            {
                const auto& f = loop.points.front();
                const auto& l = loop.points.back();
                if (std::abs(f.x - l.x) < 1e-8 && std::abs(f.y - l.y) < 1e-8)
                    loop.points.pop_back();
                else
                    break;
            }

            if (loop.points.size() >= 3)
                uvLoops.push_back(std::move(loop));
        }
        if (uvLoops.empty()) return;

        // ── Step 2: Compute U subdivision threshold ──────────────────────────
        // Only the circumferential (U) direction is curved on the surface.
        double minU = 1e30, maxU = -1e30;
        for (const auto& loop : uvLoops)
            for (const auto& uv : loop.points)
            {
                minU = std::min(minU, uv.x);
                maxU = std::max(maxU, uv.x);
            }
        double arcDeg = maxU - minU;
        int uSteps = std::max(static_cast<int>(circleSegments * arcDeg / 360.0), 4);
        double maxEdgeAngle = arcDeg / uSteps;

        // Scale V so that CDT sees comparable coordinate ranges in U and V,
        // avoiding degenerate triangulations from extreme aspect ratios.
        double vScale = (totalArc > 1e-12 && arcDeg > 1e-9)
            ? (arcDeg / totalArc) : 1.0;

        // ── Step 3: Build CDT vertices and constrained edges ─────────────────
        std::vector<CDT::V2d<double>> cdtVerts;
        std::vector<CDT::Edge> cdtEdges;

        for (const auto& loop : uvLoops)
        {
            size_t n = loop.points.size();
            std::vector<uint32_t> loopIndices;

            for (size_t i = 0; i < n; i++)
            {
                const auto& pi = loop.points[i];
                const auto& pj = loop.points[(i + 1) % n];

                // Add this boundary vertex (v scaled for CDT)
                loopIndices.push_back(static_cast<uint32_t>(cdtVerts.size()));
                cdtVerts.push_back(CDT::V2d<double>::make(pi.x, pi.y * vScale));

                // Subdivide edge if circumferential span is large
                double dU = std::abs(pj.x - pi.x);
                int nSub = (maxEdgeAngle > 1e-9)
                    ? static_cast<int>(std::ceil(dU / maxEdgeAngle))
                    : 1;
                for (int k = 1; k < nSub; k++)
                {
                    double t = static_cast<double>(k) / nSub;
                    double mu = pi.x + t * (pj.x - pi.x);
                    double mv = pi.y + t * (pj.y - pi.y);
                    loopIndices.push_back(static_cast<uint32_t>(cdtVerts.size()));
                    cdtVerts.push_back(CDT::V2d<double>::make(mu, mv * vScale));
                }
            }

            // Close the loop with constrained edges
            for (size_t i = 0; i < loopIndices.size(); i++)
                cdtEdges.push_back(CDT::Edge(
                    loopIndices[i],
                    loopIndices[(i + 1) % loopIndices.size()]));
        }

        if (cdtVerts.size() < 3) return;

        // ── Step 4: Constrained Delaunay Triangulation ───────────────────────
        CDT::RemoveDuplicatesAndRemapEdges(cdtVerts, cdtEdges);

        // Remove degenerate edges (both endpoints identical after dedup)
        cdtEdges.erase(
            std::remove_if(cdtEdges.begin(), cdtEdges.end(),
                [](const CDT::Edge& e){ return e.v1() == e.v2(); }),
            cdtEdges.end());

        CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::AsProvided);
        cdt.insertVertices(cdtVerts);
        cdt.insertEdges(cdtEdges);
        cdt.eraseSuperTriangle();

        // ── Step 5: Map UV triangles to 3D (undo vScale for arc-length) ──────
        // Use point-in-polygon on triangle centroids to keep only interior
        // triangles (more robust than eraseOuterTrianglesAndHoles).
        auto insideUVPoly = [&](double pu, double pv) -> bool
        {
            int crossings = 0;
            for (const auto& loop : uvLoops)
            {
                size_t n = loop.points.size();
                for (size_t i = 0, j = n - 1; i < n; j = i++)
                {
                    double yi = loop.points[i].y, yj = loop.points[j].y;
                    if ((yi > pv) != (yj > pv))
                    {
                        double xCross = loop.points[j].x
                            + (pv - yj) / (yi - yj)
                            * (loop.points[i].x - loop.points[j].x);
                        if (pu < xCross) crossings++;
                    }
                }
            }
            return (crossings & 1) != 0;
        };

        double vScaleInv = (std::abs(vScale) > 1e-30) ? (1.0 / vScale) : 1.0;
        for (const auto& tri : cdt.triangles)
        {
            const auto& v0 = cdt.vertices[tri.vertices[0]];
            const auto& v1 = cdt.vertices[tri.vertices[1]];
            const auto& v2 = cdt.vertices[tri.vertices[2]];

            // Centroid in original (unscaled) UV space
            double cu = (v0.x + v1.x + v2.x) / 3.0;
            double cv = (v0.y + v1.y + v2.y) / 3.0 * vScaleInv;
            if (!insideUVPoly(cu, cv)) continue;

            geometry.AddFace(
                uvTo3D(v0.x, v0.y * vScaleInv),
                uvTo3D(v1.x, v1.y * vScaleInv),
                uvTo3D(v2.x, v2.y * vScaleInv));
        }

#ifdef _DEBUG
        if (geometry.numFaces > 300) {
            std::cout << "geometry.numFaces: " << geometry.numFaces;
        }
#endif
    }

    inline void TriangulateCylindricalSurface( IfcGeometry& geometry, const std::vector<IfcBound3D>& bounds, const IfcSurface& surface, double numCircleSegments)
    {
        spdlog::debug("[TriangulateCylindricalSurface()]");

#if defined(DEBUG_DUMP_SVG) || defined(DUMP_CSG_MESHES)
        for (size_t ii = 0; ii < bounds.size(); ++ii )
        {
            const auto& bound = bounds[ii];
            if (bound.curve.points.empty()) continue;
            std::vector<double> angles, heights;
            std::string fileName = "faceBounds" + std::to_string(ii) + ".html";
            webifc::dump::DumpCurveToHtml(bound.curve.points, fileName);
        }
#endif

        glm::dvec3 origin = surface.transformation[3];
        glm::dvec3 axisX = glm::normalize(surface.transformation[0]);
        glm::dvec3 axisY = glm::normalize(surface.transformation[1]);
        glm::dvec3 axisZ = glm::normalize(surface.transformation[2]);
        double radius = surface.CylinderSurface.Radius;

        // UV → 3D (matches RevolveCylinder: sin(u)*axisX + cos(u)*axisY)
        auto uvTo3D = [&](double u, double v) -> glm::dvec3
        {
            double uRad = glm::radians(u);
            return origin + v * axisZ
                 + std::sin(uRad) * radius * axisX
                 + std::cos(uRad) * radius * axisY;
        };

        // ── Step 1: Project bounds into UV space as separate closed loops ────
        //   u = circumferential angle (degrees)
        //   v = axial distance along axisZ
        struct UVLoop { std::vector<glm::dvec2> points; };
        std::vector<UVLoop> uvLoops;

        for (const auto& bound : bounds)
        {
            if (bound.curve.points.size() < 3) continue;
            std::vector<double> angles, heights;
            for (const auto& p : bound.curve.points)
            {
                glm::dvec3 rel = p - origin;
                double u = VectorToAngle3D(glm::dot(axisX, rel), glm::dot(axisY, rel));
                NormalizeAngle(u);
                angles.push_back(u);
                heights.push_back(glm::dot(axisZ, rel));
            }
            UnwrapAngles(angles);

            UVLoop loop;
            for (size_t i = 0; i < angles.size(); i++)
                loop.points.push_back({angles[i], heights[i]});

            // Remove duplicate closing vertex(es)
            while (loop.points.size() >= 2)
            {
                const auto& f = loop.points.front();
                const auto& l = loop.points.back();
                if (std::abs(f.x - l.x) < 1e-8 && std::abs(f.y - l.y) < 1e-8)
                    loop.points.pop_back();
                else
                    break;
            }

            if (loop.points.size() >= 3)
                uvLoops.push_back(std::move(loop));
        }
        if (uvLoops.empty()) return;

        // ── Step 2: Compute U subdivision threshold ──────────────────────────
        // Only the circumferential (U) direction is curved; axial (V) is straight.
        // Edges with large U span need subdivision to approximate the curvature.
        double minU = 1e30, maxU = -1e30;
        for (const auto& loop : uvLoops)
            for (const auto& uv : loop.points)
            {
                minU = std::min(minU, uv.x);
                maxU = std::max(maxU, uv.x);
            }
        double arcDeg = maxU - minU;
        int uSteps = std::max(static_cast<int>(numCircleSegments * arcDeg / 360.0), 4);
        double maxEdgeAngle = arcDeg / uSteps;

        // ── Step 3: Build CDT vertices and constrained edges ─────────────────
        // Each bound becomes a closed loop of constrained edges.  Edges with
        // large U span are subdivided to keep the chord error small.
        std::vector<CDT::V2d<double>> cdtVerts;
        std::vector<CDT::Edge> cdtEdges;

        for (const auto& loop : uvLoops)
        {
            size_t n = loop.points.size();
            std::vector<uint32_t> loopIndices;

            for (size_t i = 0; i < n; i++)
            {
                const auto& pi = loop.points[i];
                const auto& pj = loop.points[(i + 1) % n];

                // Add this boundary vertex
                loopIndices.push_back(static_cast<uint32_t>(cdtVerts.size()));
                cdtVerts.push_back(CDT::V2d<double>::make(pi.x, pi.y));

                // Subdivide edge if circumferential span is large
                double dU = std::abs(pj.x - pi.x);
                int nSub = (maxEdgeAngle > 1e-9)
                    ? static_cast<int>(std::ceil(dU / maxEdgeAngle))
                    : 1;
                for (int k = 1; k < nSub; k++)
                {
                    double t = static_cast<double>(k) / nSub;
                    loopIndices.push_back(static_cast<uint32_t>(cdtVerts.size()));
                    cdtVerts.push_back(CDT::V2d<double>::make(
                        pi.x + t * (pj.x - pi.x),
                        pi.y + t * (pj.y - pi.y)));
                }
            }

            // Close the loop with constrained edges
            for (size_t i = 0; i < loopIndices.size(); i++)
                cdtEdges.push_back(CDT::Edge(
                    loopIndices[i],
                    loopIndices[(i + 1) % loopIndices.size()]));
        }

        if (cdtVerts.size() < 3) return;

        // ── Step 4: Constrained Delaunay Triangulation ───────────────────────
        CDT::RemoveDuplicatesAndRemapEdges(cdtVerts, cdtEdges);

        // Remove degenerate edges (both endpoints identical after dedup)
        cdtEdges.erase(
            std::remove_if(cdtEdges.begin(), cdtEdges.end(),
                [](const CDT::Edge& e){ return e.v1() == e.v2(); }),
            cdtEdges.end());

        CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::AsProvided);
        cdt.insertVertices(cdtVerts);
        cdt.insertEdges(cdtEdges);
        cdt.eraseSuperTriangle();

        // ── Step 5: Map UV triangles to 3D ───────────────────────────────────
        // Use point-in-polygon on triangle centroids to keep only interior
        // triangles (more robust than eraseOuterTrianglesAndHoles).
        auto insideUVPoly = [&](double pu, double pv) -> bool
        {
            int crossings = 0;
            for (const auto& loop : uvLoops)
            {
                size_t n = loop.points.size();
                for (size_t i = 0, j = n - 1; i < n; j = i++)
                {
                    double yi = loop.points[i].y, yj = loop.points[j].y;
                    if ((yi > pv) != (yj > pv))
                    {
                        double xCross = loop.points[j].x
                            + (pv - yj) / (yi - yj)
                            * (loop.points[i].x - loop.points[j].x);
                        if (pu < xCross) crossings++;
                    }
                }
            }
            return (crossings & 1) != 0;
        };

        for (const auto& tri : cdt.triangles)
        {
            const auto& v0 = cdt.vertices[tri.vertices[0]];
            const auto& v1 = cdt.vertices[tri.vertices[1]];
            const auto& v2 = cdt.vertices[tri.vertices[2]];

            // Centroid in original UV space
            double cu = (v0.x + v1.x + v2.x) / 3.0;
            double cv = (v0.y + v1.y + v2.y) / 3.0;
            if (!insideUVPoly(cu, cv)) continue;

            geometry.AddFace(
                uvTo3D(v0.x, v0.y),
                uvTo3D(v1.x, v1.y),
                uvTo3D(v2.x, v2.y));
        }
    }

    inline void TriangulateExtrusion( IfcGeometry& geometry, const std::vector<IfcBound3D>&, const IfcSurface& surface)
    {
        spdlog::debug("[TriangulateExtrusion()]");

        double len = surface.ExtrusionSurface.Length;
        glm::dvec3 dir = surface.ExtrusionSurface.Direction;

        if (!surface.ExtrusionSurface.Profile.isComposite)
        {
            const std::vector<glm::dvec3>& pts = surface.ExtrusionSurface.Profile.curve.points;
            auto geom = bimGeometry::Extrude(pts, dir, len);
            AppendGeometry(geometry, geom);
        }
        else
        {
            for (const auto& profile : surface.ExtrusionSurface.Profile.profiles)
            {
                const std::vector<glm::dvec3>& pts = profile.curve.points;
                auto geom = bimGeometry::Extrude(pts, dir, len);
                AppendGeometry(geometry, geom);
            }
        }
    }

    inline void TriangulateBspline( IfcGeometry& geometry, const std::vector<IfcBound3D>& bounds, const IfcSurface& surface, double scaling)
    {
        spdlog::debug("[TriangulateBspline()]");
        Nurbs nurbs{ geometry, bounds, surface, scaling };
        nurbs.fill_geometry();
    }
}