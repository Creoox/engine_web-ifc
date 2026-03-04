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

    inline void TriangulateCylindricalSurface( IfcGeometry& geometry, const std::vector<IfcBound3D>& bounds, const IfcSurface& surface, double numRots)
    {
        spdlog::debug("[TriangulateCylindricalSurface()]");

        auto groups = GroupBoundsByIndex(bounds);

        glm::dvec3 origin = surface.transformation[3];
        glm::dvec3 axisX = glm::normalize(surface.transformation[0]);
        glm::dvec3 axisY = glm::normalize(surface.transformation[1]);
        glm::dvec3 axisZ = glm::normalize(surface.transformation[2]);

        double minZ = std::numeric_limits<double>::max();
        double maxZ = -std::numeric_limits<double>::max();

        std::vector<glm::dvec3> allPoints;
        for (const auto& [_, pts] : groups)
        {
            for (const auto& p : pts)
            {
                glm::dvec3 v = p - origin;
                double dz = glm::dot(axisZ, v);

                minZ = std::min(minZ, dz);
                maxZ = std::max(maxZ, dz);

                allPoints.push_back(p);
            }
        }

        if (allPoints.empty())
            return;

        auto [startDeg, endDeg] = ComputeAngleInterval( allPoints, origin, axisX, axisY);

        auto geom = bimGeometry::RevolveCylinder( surface.transformation, startDeg, endDeg,
            minZ, maxZ, numRots, surface.CylinderSurface.Radius);

        AppendGeometry(geometry, geom);
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