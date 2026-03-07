/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

#pragma once

#include "../nurbs.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <cstdint>
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

        glm::dvec3 origin = surface.transformation[3];
        glm::dvec3 axisX = glm::normalize(surface.transformation[0]);
        glm::dvec3 axisY = glm::normalize(surface.transformation[1]);
        glm::dvec3 axisZ = glm::normalize(surface.transformation[2]);
        double radius = surface.CylinderSurface.Radius;

        // ── Step 1: project boundary into cylinder UV space ──────────────────────
        //   u = circumferential angle (degrees, same convention as RevolveCylinder)
        //   v = axial distance along axisZ
        std::vector<glm::dvec2> uvBoundary;
        for (const auto& bound : bounds)
        {
            if (bound.curve.points.empty()) continue;
            std::vector<double> angles, heights;
            for (const auto& p : bound.curve.points)
            {
                glm::dvec3 rel = p - origin;
                double dx = glm::dot(axisX, rel);
                double dy = glm::dot(axisY, rel);
                double dz = glm::dot(axisZ, rel);
                double u = VectorToAngle3D(dx, dy);
                NormalizeAngle(u);
                angles.push_back(u);
                heights.push_back(dz);
            }
            UnwrapAngles(angles);
            for (size_t i = 0; i < angles.size(); i++)
                uvBoundary.push_back({angles[i], heights[i]});
        }
        if (uvBoundary.size() < 3) return;

        // ── Step 2: UV bounding box ───────────────────────────────────────────────
        double minU = uvBoundary[0].x, maxU = uvBoundary[0].x;
        double minV = uvBoundary[0].y, maxV = uvBoundary[0].y;
        for (const auto& uv : uvBoundary)
        {
            minU = std::min(minU, uv.x);  maxU = std::max(maxU, uv.x);
            minV = std::min(minV, uv.y);  maxV = std::max(maxV, uv.y);
        }

        // ── Step 3: UV → 3D (matches RevolveCylinder: sin(u)*axisX + cos(u)*axisY) ──
        auto uvTo3D = [&](double u, double v) -> glm::dvec3
        {
            double uRad = glm::radians(u);
            return origin + v * axisZ
                 + std::sin(uRad) * radius * axisX
                 + std::cos(uRad) * radius * axisY;
        };

        // ── Step 4: scanline tessellation ────────────────────────────────────────
        // For each U column we intersect the UV boundary polygon with the vertical
        // line u = const to find the exact V extent.
        //
        // Cylinder is straight in the V (axial) direction, so one row of quads
        // (vSteps=1) is geometrically exact — V subdivision cannot improve accuracy.
        // Only U needs subdivision because the surface is curved circumferentially.
        double arcDeg = maxU - minU;
        int uSteps = std::max(static_cast<int>(numCircleSegments * arcDeg / 360.0), 4);
        double uStep = arcDeg / uSteps;

        // Collect distinct V values from boundary vertices that lie on
        // near-vertical (axial) edges.  These are "side" points that adjacent
        // faces may share, so the strip must pass through them to avoid T-junctions.
        std::vector<double> sideVals;
        for (size_t i = 0, n = uvBoundary.size(); i < n; i++)
        {
            size_t j = (i + 1) % n;
            double du = std::abs(uvBoundary[j].x - uvBoundary[i].x);
            double dv = std::abs(uvBoundary[j].y - uvBoundary[i].y);
            // Edge is "near-vertical" if its U span is tiny relative to V span
            if (du < 1e-6 * (dv + 1e-12))
            {
                sideVals.push_back(uvBoundary[i].y);
                sideVals.push_back(uvBoundary[j].y);
            }
        }
        std::sort(sideVals.begin(), sideVals.end());
        sideVals.erase(std::unique(sideVals.begin(), sideVals.end(),
                        [](double a, double b){ return std::abs(a - b) < 1e-9; }),
                       sideVals.end());

        auto getVRangeAtU = [&](double u) -> std::pair<double,double>
        {
            std::vector<double> crossings;
            crossings.reserve(4);
            const size_t n = uvBoundary.size();
            for (size_t i = 0, j = n - 1; i < n; j = i++)
            {
                double ui = uvBoundary[i].x, uj = uvBoundary[j].x;
                double vi = uvBoundary[i].y, vj = uvBoundary[j].y;
                if (std::abs(ui - uj) < 1e-9) continue;       // skip vertical edges
                if ((ui > u) == (uj > u))      continue;       // edge does not straddle u
                double t = (u - uj) / (ui - uj);
                if (t < -1e-9 || t > 1.0 + 1e-9) continue;
                crossings.push_back(vj + t * (vi - vj));
            }
            if (crossings.size() < 2) return {0.0, -1.0};     // outside polygon
            auto mm = std::minmax_element(crossings.begin(), crossings.end());
            return {*mm.first, *mm.second};
        };

        for (int iu = 0; iu < uSteps; iu++)
        {
            double u0 = minU + iu       * uStep;
            double u1 = minU + (iu + 1) * uStep;

            auto [vA0, vB0] = getVRangeAtU(u0);
            auto [vA1, vB1] = getVRangeAtU(u1);
            auto [vAm, vBm] = getVRangeAtU((u0 + u1) * 0.5);

            if (vAm > vBm) continue;                           // column outside polygon
            // Fallback for u0/u1 that sit exactly on a vertical boundary edge
            if (vA0 > vB0) { vA0 = vAm; vB0 = vBm; }
            if (vA1 > vB1) { vA1 = vAm; vB1 = vBm; }

            // Build V subdivision for this column: [vA..vB] plus any side
            // vertices that fall inside the range.
            double vLo0 = std::min(vA0, vA1), vHi0 = std::max(vB0, vB1);
            std::vector<double> vLevels;
            vLevels.push_back(0.0);  // t=0 (bottom)
            for (double sv : sideVals)
            {
                if (sv > vLo0 + 1e-9 && sv < vHi0 - 1e-9)
                {
                    // Convert to parametric t in [0,1] relative to each side
                    // Left side:  v = vA0 + t*(vB0-vA0)  →  t = (sv-vA0)/(vB0-vA0)
                    // Right side: v = vA1 + t*(vB1-vA1)  →  t = (sv-vA1)/(vB1-vA1)
                    // Use average t so the row aligns across the column
                    double rangeL = vB0 - vA0, rangeR = vB1 - vA1;
                    double tL = (rangeL > 1e-12) ? (sv - vA0) / rangeL : 0.5;
                    double tR = (rangeR > 1e-12) ? (sv - vA1) / rangeR : 0.5;
                    double tAvg = (tL + tR) * 0.5;
                    if (tAvg > 1e-6 && tAvg < 1.0 - 1e-6)
                        vLevels.push_back(tAvg);
                }
            }
            vLevels.push_back(1.0);  // t=1 (top)
            std::sort(vLevels.begin(), vLevels.end());
            vLevels.erase(std::unique(vLevels.begin(), vLevels.end(),
                            [](double a, double b){ return std::abs(a - b) < 1e-9; }),
                          vLevels.end());

            for (size_t iv = 0; iv + 1 < vLevels.size(); iv++)
            {
                double t0 = vLevels[iv];
                double t1 = vLevels[iv + 1];

                double v00 = vA0 + t0 * (vB0 - vA0);
                double v01 = vA0 + t1 * (vB0 - vA0);
                double v10 = vA1 + t0 * (vB1 - vA1);
                double v11 = vA1 + t1 * (vB1 - vA1);

                geometry.AddFace(uvTo3D(u0, v00), uvTo3D(u0, v01), uvTo3D(u1, v10));
                geometry.AddFace(uvTo3D(u1, v10), uvTo3D(u0, v01), uvTo3D(u1, v11));
            }
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