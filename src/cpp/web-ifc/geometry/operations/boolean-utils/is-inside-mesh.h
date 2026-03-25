/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

/*
    This function was revised on 2024-03-17.  Further review is needed once we have an overview
    of the entire package.

    Parameter normal is assumed to be normalised.
    Function computeNormal returns a normalised vector.

    The constants TOLERANCE_PLANE_DEVIATION and toleranceParallel are defined in eps.h.
*/
#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <stack>
#include "intersect-ray-tri.h"
#include "geometry.h"
#include "bvh.h"

using Vec = glm::dvec3;

namespace fuzzybools
{
    enum class MeshLocation
    {
        INSIDE,
        OUTSIDE,
        BOUNDARY
    };

    struct InsideResult
    {
        MeshLocation loc;
        Vec normal;
    };

    inline InsideResult isInsideMesh(
        const Vec &pt,
        Vec normal,
        const Geometry &g,
        BVH &bvh,
        Vec dir = Vec(1.0, 1.1, 1.4),
        bool UNION = false)
    {
        int winding = 0;
        dir = dir + Vec(0.02, 0.01, 0.04); // Randomly changing the normal to create a truly random direction for raytrace
        InsideResult result;
        result.loc = MeshLocation::BOUNDARY;
        result.normal = glm::dvec3(0);

        bool hasResult = bvh.IntersectRay(pt, dir, [&](uint32_t i) -> bool
                                          {
                Face f = g.GetFace(i);
                const Vec a = g.GetPoint(f.i0);
                const Vec b = g.GetPoint(f.i1);
                const Vec c = g.GetPoint(f.i2);

                Vec intersection;
                double distance;
                double d_plane;
                bool hasIntersection = intersect_ray_triangle(pt, pt + dir, a, b, c, intersection, distance, d_plane, true);
                if (hasIntersection)
                {
                    Vec otherNormal = computeNormal(a, b, c);  // normalised
                    double d = glm::dot(otherNormal, dir);
                    double dn = glm::dot(otherNormal, normal);
                    if (std::fabs(d_plane) < _TOLERANCE_PLANE_DEVIATION)
                    {
                        if (dn > 1.0 - toleranceParallel)
                        {
                            // The normals point in the same direction, which means that the boundary is an inside boundary.
                            result.loc = MeshLocation::BOUNDARY;
                            result.normal = normal;
                            return true;
                        }
                        else if (dn < -1.0 + toleranceParallel)
                        {
                            // The normals point in opposite directions, which means that the boundary is an outside boundary.
                            if(!UNION)
                            {
                                result.loc = MeshLocation::OUTSIDE;
                                result.normal = normal;
                                return true;
                            }
                            else
                            {
                                result.loc = MeshLocation::BOUNDARY;
                                result.normal = normal;
                                return true;
                            }
                        }
                        else
                        {
                            result.loc = MeshLocation::BOUNDARY;
                            result.normal = otherNormal;
                            return true;
                        }
                    }

                    winding++;
                }
                else
                {
                    // Edge-hit post-check.
                    //
                    // intersect_ray_triangle returns false for three reasons:
                    //   (A) ray nearly parallel to face  (|NdotDir| < toleranceParallelTight)
                    //   (B) face too far behind ray      (t < -_TOLERANCE_BACK_DEVIATION_DISTANCE)
                    //   (C) inside-outside test failed   (pt projects outside triangle)
                    //
                    // Case C with t ≈ 0 is the "edge-hit" case: pt lies exactly on the plane
                    // of this face (and of an adjacent face), but the floating-point cross-
                    // product in the inside-outside test gives valdot = −ε for BOTH faces,
                    // so hasIntersection stays false.  Without this fix the winding counter
                    // never sees the near face; instead it counts the far face of the other
                    // solid, giving winding=1 -> INSIDE -> face wrongly dropped -> open mesh.
                    //
                    // Distinguish case C from A/B by re-computing NdotDir and t, then confirm
                    // with a barycentric test that pt actually sits on the triangle boundary.
                    {
                        const Vec n_unnorm = glm::cross(b - a, c - a);
                        const double NdotDir = glm::dot(n_unnorm, dir);
                        if (std::fabs(NdotDir) >= toleranceParallelTight)          // not case A
                        {
                            const double t_eh = -glm::dot(n_unnorm, pt - a) / NdotDir;
                            if (t_eh >= -_TOLERANCE_BACK_DEVIATION_DISTANCE)        // not case B
                            {
                                // Case C: the face was in range but the inside-outside test
                                // rejected pt.  Only treat this as an edge-hit when pt is
                                // truly near the face plane (planeDist < tolerance).
                                // When t_eh ≠ 0 the orthogonal projection proj differs from
                                // the ray-intersection point p by −t_eh·dir_tangent; for
                                // oblique rays against a tilted face this divergence can make
                                // the barycentric test pass even when p is outside the triangle.
                                const Vec faceNormal = computeNormal(a, b, c);
                                const double planeDist = std::abs(glm::dot(faceNormal, pt - a));
                                if (planeDist < _TOLERANCE_PLANE_DEVIATION)
                                {
                                    const Vec proj = pt - glm::dot(pt - a, faceNormal) * faceNormal;
                                    const Vec v0 = b - a, v1 = c - a, v2 = proj - a;
                                    const double d00   = glm::dot(v0, v0);
                                    const double d01   = glm::dot(v0, v1);
                                    const double d11   = glm::dot(v1, v1);
                                    const double d20   = glm::dot(v2, v0);
                                    const double d21   = glm::dot(v2, v1);
                                    const double denom = d00 * d11 - d01 * d01;
                                    if (std::abs(denom) > 1e-20)
                                    {
                                        const double bu = (d11 * d20 - d01 * d21) / denom;
                                        const double bv = (d00 * d21 - d01 * d20) / denom;
                                        const double bw = 1.0 - bu - bv;
                                        // Allow a small negative tolerance so that a point
                                        // sitting exactly on an edge (theoretically bary = 0)
                                        // passes despite floating-point rounding.
                                        constexpr double BARY_EPS = 1e-8;
                                        if (bu >= -BARY_EPS && bv >= -BARY_EPS && bw >= -BARY_EPS)
                                        {
                                            // Confirmed edge-hit: apply coplanar classification.
                                            const double dn = glm::dot(faceNormal, normal);
                                            if (dn > 1.0 - toleranceParallel)
                                            {
                                                result.loc = MeshLocation::BOUNDARY;
                                                result.normal = normal;
                                                return true;
                                            }
                                            else if (dn < -1.0 + toleranceParallel)
                                            {
                                                result.loc = UNION ? MeshLocation::BOUNDARY
                                                                   : MeshLocation::OUTSIDE;
                                                result.normal = normal;
                                                return true;
                                            }
                                            else
                                            {
                                                result.loc = MeshLocation::BOUNDARY;
                                                result.normal = faceNormal;
                                                return true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Continue to search.
                return false; });

        if (hasResult)
        {
            return result;
        }

        result.loc = winding % 2 == 1 ? MeshLocation::INSIDE : MeshLocation::OUTSIDE;

        return result;
    }
}
