/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

#pragma once

#ifdef _DEBUG
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

    inline void TriangulateRevolution( IfcGeometry& geometry, const std::vector<IfcBound3D>& bounds, const IfcSurface& surface, double numRots)
    {
        spdlog::debug("[TriangulateRevolution()]");
        auto groups = GroupBoundsByIndex(bounds);

        std::vector<glm::dvec3> centroids;
        for (const auto& [_, pts] : groups)
        {
            if (pts.empty())
                continue;

            glm::dvec3 sum(0.0);
            for (const auto& p : pts)
                sum += p;

            centroids.push_back(sum / double(pts.size()));
        }

        if (centroids.empty())
            return;

        glm::dmat4 transform;
        transform[3] = surface.RevolutionSurface.Direction[3];
        transform[0] = glm::normalize(surface.RevolutionSurface.Direction[0]);
        transform[1] = glm::normalize(surface.RevolutionSurface.Direction[1]);
        transform[2] = glm::normalize(surface.RevolutionSurface.Direction[2]);

        auto [startDeg, endDeg] = ComputeAngleInterval( centroids, glm::dvec3(transform[3]), glm::dvec3(transform[0]), glm::dvec3(transform[1]));
        auto geom = bimGeometry::Revolution( transform, startDeg, endDeg, surface.RevolutionSurface.Profile.curve.points, numRots);

        AppendGeometry(geometry, geom);
    }

    inline void TriangulateCylindricalSurface( IfcGeometry& geometry, const std::vector<IfcBound3D>& bounds, const IfcSurface& surface, double numCircleSegments)
    {
        spdlog::debug("[TriangulateCylindricalSurface()]");

#ifdef _DEBUG
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

        CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::AsProvided);
        cdt.insertVertices(cdtVerts);
        cdt.insertEdges(cdtEdges);
        cdt.eraseOuterTrianglesAndHoles();

        // ── Step 5: Map UV triangles to 3D ───────────────────────────────────
        for (const auto& tri : cdt.triangles)
        {
            const auto& v0 = cdt.vertices[tri.vertices[0]];
            const auto& v1 = cdt.vertices[tri.vertices[1]];
            const auto& v2 = cdt.vertices[tri.vertices[2]];
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