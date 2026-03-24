/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include <glm/glm.hpp>
#include <cstdint>

#include "util.h"
#include "is-inside-mesh.h"
#include "geometry.h"
#include "bvh.h"

namespace fuzzybools
{
    static void doubleClipSingleMesh(Geometry& mesh, BVH& bvh1, BVH& bvh2, Geometry& result)
    {  
        #ifdef CSG_DEBUG_OUTPUT
            std::vector<std::vector<glm::dvec2>> edgesPrinted;
        #endif

        // Hash-map based O(1) duplicate-triangle detection, replaces the O(N) linear scan.
        //
        // hashDouble: MurmurHash finalizer applied to the bit pattern of a double so that bitwise-identical 
        // values always hash the same.
        // canonicalTriHash: cyclic-rotation-invariant hash of (a,b,c) -
        //   invariant under (a,b,c)->(b,c,a)->(c,a,b) but NOT reflection, matching the three cyclic equalities
        // tested by the dedup check.
        auto hashDouble = [](double x) -> std::size_t {
            std::uint64_t bits = 0;
            std::memcpy(&bits, &x, sizeof(bits));
            bits ^= bits >> 33;
            bits *= 0xff51afd7ed558ccdULL;
            bits ^= bits >> 33;
            return static_cast<std::size_t>(bits);
        };
        auto hashPt = [&](const glm::dvec3& v) -> std::size_t {
            std::size_t h = hashDouble(v.x);
            h ^= hashDouble(v.y) + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= hashDouble(v.z) + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        };
        auto canonicalTriHash = [&](const glm::dvec3& pa,
                                    const glm::dvec3& pb,
                                    const glm::dvec3& pc) -> std::size_t {
            const std::size_t ha = hashPt(pa), hb = hashPt(pb), hc = hashPt(pc);
            std::size_t h0, h1, h2;
            if      (ha <= hb && ha <= hc) { h0 = ha; h1 = hb; h2 = hc; }
            else if (hb <= ha && hb <= hc) { h0 = hb; h1 = hc; h2 = ha; }
            else                           { h0 = hc; h1 = ha; h2 = hb; }
            std::size_t h = h0;
            h ^= h1 + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= h2 + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        };
        // canonical-triangle-hash -> face indices already seen
        std::unordered_map<std::size_t, std::vector<int>> triHashMap;
        triHashMap.reserve(mesh.data / 4 + 16);
        // Track which faces were accepted (added to result)
        std::unordered_set<int> acceptedFaces;

        for(auto &plane: mesh.planes)
        {
            result.hasPlanes = true;
            result.planes.push_back(plane);
        }

        for (uint32_t i = 0; i < mesh.data; i++)
        {
            bool doit = false;
            Face tri = mesh.GetFace(i);
            glm::dvec3 a = mesh.GetPoint(tri.i0);
            glm::dvec3 b = mesh.GetPoint(tri.i1);
            glm::dvec3 c = mesh.GetPoint(tri.i2);

            auto aabb = mesh.GetFaceBox(i);

            if (!aabb.intersects(bvh2.box))
            {
                // Why is this commented?

                // when subtracting, if box is outside the second operand, its guaranteed to remain
                // result.AddFace(a, b, c);
                //continue;
            }
            else if (!aabb.intersects(bvh1.box))
            {
                // when subtracting, if box is outside the first operand, it won't remain ever
                //continue;
            }

            bool doNext = true;

            // O(1) expected: hash-map lookup
            // Skip duplicates ONLY if the earlier occurrence was accepted.
            // If it was rejected (both-boundary-same-dir), allow the
            // duplicate to be classified -- it may come from the other
            // operand and get a different result.
            {
                const std::size_t triKey = canonicalTriHash(a, b, c);
                auto it = triHashMap.find(triKey);
                if (it != triHashMap.end())
                {
                    for (int prevIdx : it->second)
                    {
                        Face tri_temp = mesh.GetFace(prevIdx);
                        glm::dvec3 at = mesh.GetPoint(tri_temp.i0);
                        glm::dvec3 bt = mesh.GetPoint(tri_temp.i1);
                        glm::dvec3 ct = mesh.GetPoint(tri_temp.i2);

                        if((equals(at,a, EPS_MINISCULE) &&  equals(bt,b, EPS_MINISCULE) && equals(ct,c, EPS_MINISCULE))
                        || (equals(at,b, EPS_MINISCULE) &&  equals(bt,c, EPS_MINISCULE) && equals(ct,a, EPS_MINISCULE))
                        || (equals(at,c, EPS_MINISCULE) &&  equals(bt,a, EPS_MINISCULE) && equals(ct,b, EPS_MINISCULE)))
                        {
                            // Only skip if the earlier copy was ACCEPTED
                            if (acceptedFaces.count(prevIdx))
                            {
                                doNext = false;
                            }
                            break;
                        }
                    }
                }
                if (doNext)
                    triHashMap[triKey].push_back(static_cast<int>(i));
            }

            if(!doNext)
            {
                continue;
            }

            glm::dvec3 n = computeNormal(a, b, c);

            glm::dvec3 triCenter = (a + b * 1.02 + c * 1.03) * 1.0 / 3.05; // Using true centroid could cause issues (#540)

            auto isInsideTarget = MeshLocation::INSIDE;

            Vec raydir = computeNormal(a, b, c);

            // This is an example about how to debug specific triangles in specific boolean operations
            // if ((i == 49 || i == 53 || i == 84) && _BOOLSTATUS == 66)
            // {
            //     doit = false; // This assignation is useless just to add some content
            // }

            auto isInside1Loc = isInsideMesh(triCenter, n, *bvh1.ptr, bvh1, raydir);
            auto isInside2Loc = isInsideMesh(triCenter, n, *bvh2.ptr, bvh2, raydir);

            Vec extraDir1 = glm::normalize(raydir + Vec(0.02,0.01,0.04));
            Vec extraDir2 = glm::normalize(raydir + Vec(0.20,-0.1,0.40));

            auto isInside1Loc_B = isInsideMesh(triCenter, n, *bvh1.ptr, bvh1, extraDir1);
            auto isInside2Loc_B = isInsideMesh(triCenter, n, *bvh2.ptr, bvh2, extraDir1);

            if(isInside1Loc.loc != isInside1Loc_B.loc)
            {
                auto isInside1Loc_C = isInsideMesh(triCenter, n, *bvh1.ptr, bvh1, extraDir2);
                if(isInside1Loc_C.loc == isInside1Loc_B.loc){isInside1Loc = isInside1Loc_B;}
                else if(isInside1Loc_B.loc != isInside1Loc_C.loc && isInside1Loc.loc != isInside1Loc_C.loc)
                {
                    isInside1Loc = isInside1Loc_B;
                }
            }

            if(isInside2Loc.loc != isInside2Loc_B.loc)
            {
                auto isInside2Loc_C = isInsideMesh(triCenter, n, *bvh2.ptr, bvh2, extraDir2);
                if(isInside2Loc_C.loc == isInside2Loc_B.loc){isInside2Loc = isInside2Loc_B;}
                else if(isInside2Loc_B.loc != isInside2Loc_C.loc && isInside2Loc.loc != isInside2Loc_C.loc)
                {
                    isInside2Loc = isInside2Loc_B;
                }
            }

            auto isInside1 = isInside1Loc.loc;
            auto isInside2 = isInside2Loc.loc;

            if (isInside1 == MeshLocation::OUTSIDE && isInside2 == MeshLocation::OUTSIDE)
            {
            }
            if (isInside1 != MeshLocation::BOUNDARY && isInside2 != MeshLocation::BOUNDARY)
            {
            }
            else if (isInside1 == MeshLocation::BOUNDARY && isInside2 == MeshLocation::BOUNDARY)
            {
                auto dot = glm::dot(isInside1Loc.normal, isInside2Loc.normal);
                if (dot < 0)
                {
                    // Opposing normals: A-B boundary face, keep
                    result.AddFace(a, b, c, tri.pId);
                    doit = true;
                }
                else
                {
                    // Same-direction normals: coplanar faces of A and B.
                    // This face should be KEPT if there is solid (A minus B)
                    // material behind it (in the -normal direction).
                    // Probe a point slightly behind the face and test:
                    //   in A AND NOT in B  -->  solid behind  -->  keep
                    //   in A AND in B      -->  void behind   -->  discard
                    //   not in A           -->  exterior      -->  discard
                    // Use multiple probe distances to handle thin flanges.
                    bool keepFace = false;
                    Vec probeN = glm::length(n) > 0.5 ? n : isInside1Loc.normal;
                    double probeDists[] = {2e-3, 1e-2, 5e-2};
                    for (double pd : probeDists)
                    {
                        Vec probeP = triCenter - probeN * pd;
                        auto pA = isInsideMesh(probeP, probeN, *bvh1.ptr, bvh1, raydir);
                        auto pB = isInsideMesh(probeP, probeN, *bvh2.ptr, bvh2, raydir);
                        if (pA.loc == MeshLocation::INSIDE && pB.loc != MeshLocation::INSIDE)
                        {
                            keepFace = true;
                            break;
                        }
                        if (pA.loc == MeshLocation::INSIDE && pB.loc == MeshLocation::INSIDE)
                            break; // void behind, discard
                    }
                    if (keepFace)
                    {
                        result.AddFace(a, b, c, tri.pId);
                        doit = true;
                    }
                }
            }
            else
            {
                if (isInside2 == MeshLocation::INSIDE || isInside1 == MeshLocation::OUTSIDE)
                {
                }
                else
                {
                    if (isInside2 == MeshLocation::BOUNDARY && isInside1 == MeshLocation::INSIDE)
                    {
                        if (glm::dot(n, isInside2Loc.normal) < 0)
                        {
                            result.AddFace(a, b, c, tri.pId);
                            doit = true;
                        }
                        else
                        {
                            result.AddFace(b, a, c, tri.pId);
                            doit = true;
                        }

                    }
                    else if (isInside1 == MeshLocation::BOUNDARY)
                    {
                        if (glm::dot(n, isInside1Loc.normal) < 0)
                        {
                            result.AddFace(b, a, c, tri.pId);
                            doit = true;
                        }
                        else
                        {
                            result.AddFace(a, b, c, tri.pId);
                            doit = true;
                        }

                    }
                    else
                    {
                        result.AddFace(a, b, c, tri.pId);
                        doit = true;

                    }
                }
            }

            if (doit)
                acceptedFaces.insert(static_cast<int>(i));

            #ifdef CSG_DEBUG_OUTPUT
            #endif
        }

        #ifdef CSG_DEBUG_OUTPUT
            // DumpSVGLines(edgesPrinted, L"final_tri.html");
        #endif

        for (uint32_t i = mesh.data; i < mesh.numFaces; i++)
        {
            Face tri = mesh.GetFace(i);
            glm::dvec3 a = mesh.GetPoint(tri.i0);
            glm::dvec3 b = mesh.GetPoint(tri.i1);
            glm::dvec3 c = mesh.GetPoint(tri.i2);
            result.AddFace(a, b, c, tri.pId);
        }
    }

    static void doubleClipSingleMesh2(Geometry& mesh, BVH& bvh1, BVH& bvh2, Geometry& result)
    {
        for(auto &plane: mesh.planes)
        {
            result.hasPlanes = true;
            result.planes.push_back(plane);
        }

        for (uint32_t i = 0; i < mesh.data; i++)
        {
            Face tri = mesh.GetFace(i);
            glm::dvec3 a = mesh.GetPoint(tri.i0);
            glm::dvec3 b = mesh.GetPoint(tri.i1);
            glm::dvec3 c = mesh.GetPoint(tri.i2);

            glm::dvec3 n = computeNormal(a, b, c);

            auto area = areaOfTriangle(a, b, c);

            glm::dvec3 triCenter = (a + b * 2.0 + c * 3.0) * 1.0 / 6.0; // Using true centroid could cause issues (#540)

            auto isInsideTarget = MeshLocation::INSIDE;

            Vec raydir = computeNormal(a, b, c);

            auto isInside1Loc = isInsideMesh(triCenter, n, *bvh1.ptr, bvh1, raydir, true);
            auto isInside2Loc = isInsideMesh(triCenter, n, *bvh2.ptr, bvh2, raydir, true);

            Vec extraDir1 = glm::normalize(Vec(1.1, 1.4, 1.2));
            Vec extraDir2 = glm::normalize(Vec(-2.1, 1.4, -3.2));

            auto isInside1Loc_B = isInsideMesh(triCenter, n, *bvh1.ptr, bvh1, extraDir1, true);
            auto isInside2Loc_B = isInsideMesh(triCenter, n, *bvh2.ptr, bvh2, extraDir1, true);

            if(isInside1Loc.loc != isInside1Loc_B.loc)
            {
                auto isInside1Loc_C = isInsideMesh(triCenter, n, *bvh1.ptr, bvh1, extraDir2, true);
                if(isInside1Loc_C.loc == isInside1Loc_B.loc){isInside1Loc = isInside1Loc_B;}
                else if(isInside1Loc_B.loc != isInside1Loc_C.loc && isInside1Loc.loc != isInside1Loc_C.loc)
                {
                    isInside1Loc = isInside1Loc_B;
                }
            }

            if(isInside2Loc.loc != isInside2Loc_B.loc)
            {
                auto isInside2Loc_C = isInsideMesh(triCenter, n, *bvh2.ptr, bvh2, extraDir2, true);
                if(isInside2Loc_C.loc == isInside2Loc_B.loc){isInside2Loc = isInside2Loc_B;}
                else if(isInside2Loc_B.loc != isInside2Loc_C.loc && isInside2Loc.loc != isInside2Loc_C.loc)
                {
                    isInside2Loc = isInside2Loc_B;
                }
            }

            auto isInside1 = isInside1Loc.loc;
            auto isInside2 = isInside2Loc.loc;

            if (isInside1 == MeshLocation::OUTSIDE && isInside2 == MeshLocation::OUTSIDE)
            {
                // both outside, no dice, should be impossible though
            }
            else if (isInside1 == MeshLocation::INSIDE || isInside2 == MeshLocation::INSIDE)
            {
                // we only keep boundaries, no dice
            }
            else if (isInside1 == MeshLocation::BOUNDARY && isInside2 == MeshLocation::BOUNDARY)
            {
                // both boundary, no dice if normals are opposite direction
                auto dot = glm::dot(isInside1Loc.normal, isInside2Loc.normal);

                // since both are on the boundary, and we're sampling the center of the tri, these two faces are coplanar
                // hence we can test dot < 0 to see if normals point the opposite way
                if (dot > 0)
                {
                    // normals face away from eachother, we can keep this face
                    // furthermore, since the first operand is the first added, we don't flip
                    result.AddFace(a, b, c, tri.pId);
                }
            }
            else if (isInside1 == MeshLocation::BOUNDARY && isInside2 == MeshLocation::OUTSIDE)
            {
                // either is a boundary, keep
                if (glm::dot(n, isInside1Loc.normal) < 0)
                {
                    result.AddFace(b, a, c, tri.pId);
                }
                else
                {
                    result.AddFace(a, b, c, tri.pId);
                }
            }
            else if (isInside2 == MeshLocation::BOUNDARY && isInside1 == MeshLocation::OUTSIDE)
            {
                // either is a boundary, keep
                if (glm::dot(n, isInside2Loc.normal) < 0)
                {
                    result.AddFace(b, a, c, tri.pId);
                }
                else
                {
                    result.AddFace(a, b, c, tri.pId);
                }
            }
            else
            {
                // neither a boundary, neither inside, neither outside, nothing left
            }
        }

        for (uint32_t i = mesh.data; i < mesh.numFaces; i++)
        {
            Face tri = mesh.GetFace(i);
            glm::dvec3 a = mesh.GetPoint(tri.i0);
            glm::dvec3 b = mesh.GetPoint(tri.i1);
            glm::dvec3 c = mesh.GetPoint(tri.i2);
            result.AddFace(a, b, c, tri.pId);
        }
    }

    static Geometry clipJoin(Geometry& mesh, BVH bvh1, BVH bvh2)
    {
        Geometry resultingMesh;

        doubleClipSingleMesh2(mesh, bvh1, bvh2, resultingMesh);

        return resultingMesh;
    }

    static Geometry clipSubtract(Geometry& mesh, BVH bvh1, BVH bvh2)
    {

        Geometry resultingMesh;

        doubleClipSingleMesh(mesh, bvh1, bvh2, resultingMesh);

        return resultingMesh;
    }
}