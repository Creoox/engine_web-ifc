/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include <glm/glm.hpp>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <array>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>

#include "util.h"
#include "boolean-budget.h"
#include "is-inside-mesh.h"
#include "geometry.h"
#include "bvh.h"

namespace fuzzybools
{
    // ------------------------------------------------------------------
    // Region-based clip classification (v5 kernel rework).
    //
    // The previous implementation classified every fragment of the
    // normalized mesh independently with tolerance-knife ray casts. Any
    // fragment whose skewed centroid produced an unlucky ray (grazing
    // hits, near-coplanar stacks, open operands) was misclassified in
    // isolation, which is exactly how membranes, flipped strips, and
    // pinholes appear in chained boolean results -- one bad fragment at a
    // time. Fragments between intersection curves share their fate by
    // construction: the inside/outside status w.r.t. the OTHER operand can
    // only change across an edge that lies on the other operand's surface.
    // So we flood-fill fragments into regions bounded by such cut edges,
    // classify each region once (majority vote over its largest fragments,
    // whose centroids are the most robust sample points), and apply the
    // decision to every fragment of the region. Misclassification now
    // requires ALL sampled fragments of a region to fail, instead of any
    // single fragment.
    // ------------------------------------------------------------------

    enum class ClipDecision : uint8_t
    {
        DROP,
        KEEP,
        KEEP_FLIP
    };

    struct ClipVotes
    {
        MeshLocation loc1 = MeshLocation::OUTSIDE;
        MeshLocation loc2 = MeshLocation::OUTSIDE;
        Vec normal1 = Vec(0);
        Vec normal2 = Vec(0);
    };

    // Triple-ray voting, extracted verbatim from the previous per-face code.
    static ClipVotes castClipVotes(const Vec& pt, const Vec& n, BVH& bvh1, BVH& bvh2, bool isUnion)
    {
        ClipVotes v;

        Vec raydir = n;
        Vec extraDir1 = isUnion ? glm::normalize(Vec(1.1, 1.4, 1.2))
                                : glm::normalize(raydir + Vec(0.02, 0.01, 0.04));
        Vec extraDir2 = isUnion ? glm::normalize(Vec(-2.1, 1.4, -3.2))
                                : glm::normalize(raydir + Vec(0.20, -0.1, 0.40));

        auto isInside1Loc = isInsideMesh(pt, n, *bvh1.ptr, bvh1, raydir, isUnion);
        auto isInside2Loc = isInsideMesh(pt, n, *bvh2.ptr, bvh2, raydir, isUnion);

        auto isInside1Loc_B = isInsideMesh(pt, n, *bvh1.ptr, bvh1, extraDir1, isUnion);
        auto isInside2Loc_B = isInsideMesh(pt, n, *bvh2.ptr, bvh2, extraDir1, isUnion);

        if (isInside1Loc.loc != isInside1Loc_B.loc)
        {
            auto isInside1Loc_C = isInsideMesh(pt, n, *bvh1.ptr, bvh1, extraDir2, isUnion);
            if (isInside1Loc_C.loc == isInside1Loc_B.loc) { isInside1Loc = isInside1Loc_B; }
            else if (isInside1Loc_B.loc != isInside1Loc_C.loc && isInside1Loc.loc != isInside1Loc_C.loc)
            {
                isInside1Loc = isInside1Loc_B;
            }
        }

        if (isInside2Loc.loc != isInside2Loc_B.loc)
        {
            auto isInside2Loc_C = isInsideMesh(pt, n, *bvh2.ptr, bvh2, extraDir2, isUnion);
            if (isInside2Loc_C.loc == isInside2Loc_B.loc) { isInside2Loc = isInside2Loc_B; }
            else if (isInside2Loc_B.loc != isInside2Loc_C.loc && isInside2Loc.loc != isInside2Loc_C.loc)
            {
                isInside2Loc = isInside2Loc_B;
            }
        }

        v.loc1 = isInside1Loc.loc;
        v.loc2 = isInside2Loc.loc;
        v.normal1 = isInside1Loc.normal;
        v.normal2 = isInside2Loc.normal;
        // NOTE: a "parity confirmation" stage (re-voting non-boundary verdicts
        // with two diverse fixed ray directions) was implemented and measured
        // here: it DOUBLED conversion runtime (extra ray casts per vote) and
        // regressed the strict suite 16 -> 14 passing files because the fixed
        // rays flip legitimate verdicts through complex geometry. Removed; the
        // jittered triplet plus region majority voting remains the validated
        // configuration.
        return v;
    }

    // Decision tables extracted verbatim from the previous implementation.
    static ClipDecision decideSubtract(const ClipVotes& v, const Vec& n)
    {
        const auto isInside1 = v.loc1;
        const auto isInside2 = v.loc2;

        if (isInside1 != MeshLocation::BOUNDARY && isInside2 != MeshLocation::BOUNDARY)
        {
            // neither boundary (covers outside/outside and inside/inside)
            return ClipDecision::DROP;
        }
        if (isInside1 == MeshLocation::BOUNDARY && isInside2 == MeshLocation::BOUNDARY)
        {
            // coplanar boundary on both: keep only when surface normals oppose
            return glm::dot(v.normal1, v.normal2) < 0 ? ClipDecision::KEEP : ClipDecision::DROP;
        }
        if (isInside2 == MeshLocation::INSIDE || isInside1 == MeshLocation::OUTSIDE)
        {
            return ClipDecision::DROP;
        }
        // Winding is deterministic from the emission invariants: Normalize's
        // CDT emits every fragment with the winding of the operand surface its
        // shared plane was built from (AddPlane aligns the plane normal with
        // the source face winding). The previous dot(n, rayHitNormal) flip
        // tests only diverged from these constants when the classification ray
        // hit a DIFFERENT coincident surface -- ray luck on stacked geometry --
        // which inverted whole shells (observed on csg-test17: the entire
        // result emitted backfacing).
        if (isInside2 == MeshLocation::BOUNDARY && isInside1 == MeshLocation::INSIDE)
        {
            // cavity wall from the subtractor surface: orient against B
            return glm::dot(n, v.normal2) < 0 ? ClipDecision::KEEP : ClipDecision::KEEP_FLIP;
        }
        if (isInside1 == MeshLocation::BOUNDARY)
        {
            // surface of the first operand: orient along A
            return glm::dot(n, v.normal1) < 0 ? ClipDecision::KEEP_FLIP : ClipDecision::KEEP;
        }
        return ClipDecision::KEEP;
    }

    static ClipDecision decideUnion(const ClipVotes& v, const Vec& n)
    {
        const auto isInside1 = v.loc1;
        const auto isInside2 = v.loc2;

        if (isInside1 == MeshLocation::OUTSIDE && isInside2 == MeshLocation::OUTSIDE)
        {
            return ClipDecision::DROP;
        }
        if (isInside1 == MeshLocation::INSIDE || isInside2 == MeshLocation::INSIDE)
        {
            return ClipDecision::DROP;
        }
        if (isInside1 == MeshLocation::BOUNDARY && isInside2 == MeshLocation::BOUNDARY)
        {
            return glm::dot(v.normal1, v.normal2) > 0 ? ClipDecision::KEEP : ClipDecision::DROP;
        }
        if (isInside1 == MeshLocation::BOUNDARY && isInside2 == MeshLocation::OUTSIDE)
        {
            return glm::dot(n, v.normal1) < 0 ? ClipDecision::KEEP_FLIP : ClipDecision::KEEP;
        }
        if (isInside2 == MeshLocation::BOUNDARY && isInside1 == MeshLocation::OUTSIDE)
        {
            return glm::dot(n, v.normal2) < 0 ? ClipDecision::KEEP_FLIP : ClipDecision::KEEP;
        }
        return ClipDecision::DROP;
    }

    // The shared region engine. planesA is the number of leading entries of
    // mesh.planes that belong to the first operand (faces carry pId into
    // that combined table); UINT32_MAX means provenance is unknown, in which
    // case every plane is treated as a potential cut plane (regions shrink,
    // correctness is unaffected).
    static void clipMeshRegions(Geometry& mesh, BVH& bvh1, BVH& bvh2, Geometry& result,
        const BooleanBudget& budget, bool isUnion, uint32_t planesA)
    {
        budget.CheckDeadline(isUnion ? "clipJoin start" : "clipSubtract start");

        for (auto& plane : mesh.planes)
        {
            result.hasPlanes = true;
            result.planes.push_back(plane);
        }

        const uint32_t nClassify = std::min(mesh.data, mesh.numFaces);

        // Scale-aware tolerances. The fixed kernel tolerances (1e-4) assume
        // meter-scale models, but cxconverter feeds geometry in the file's
        // native units -- often millimetres, where 1e-4 means 0.1 micron and
        // every coincidence-based test silently stops working (test3, a mm
        // model: the cut-edge test missed every intersection curve, regions
        // flooded across cuts, and whole opening caps were dropped on a single
        // vote while the metre-scale files behaved). Derive the working
        // tolerances from the combined mesh extent so they mean the same
        // physical size at any unit scale, with the metre-tuned constants as
        // floors so metre models keep their validated behaviour.
        const AABB meshBox = mesh.GetAABB();
        const double meshDiag = glm::length(meshBox.max - meshBox.min);
        const double cutTol = std::max(_TOLERANCE_PLANE_DEVIATION, meshDiag * 5e-6);
        // Cut-EDGE detection tolerance is deliberately coarser than cutTol: operand
        // vertices are no longer snapped onto their registered planes (the buildPlanes
        // refit is gone), so fragment vertices sit off the plane by up to the plane-merge
        // error (~diag * 5e-5). With the tight tolerance the region flood-fill leaked
        // across undetected cut curves and subtractor-surface wedges OUTSIDE the first
        // operand inherited the KEEP vote of the legitimate reveal region they touched
        // (test57-allOpenings#: solid corner flags protruding into every opening).
        // Coarser detection only shrinks regions, which is correctness-neutral by
        // construction (see the provenance comment above).
        const double cutDetectTol = std::max(cutTol, meshDiag * 5e-5);
        const double weldTol = std::max(1.0E-07, meshDiag * 1e-9);

        // ---- per-face data ------------------------------------------------
        enum : uint8_t { ST_REGION = 0, ST_PRE_KEEP = 1, ST_PRE_DROP = 2, ST_KEEP = 3, ST_KEEP_FLIP = 4 };
        std::vector<uint8_t> state(nClassify, ST_REGION);
        struct FaceData { Vec a, b, c; Vec n; Vec center; double area; uint32_t pId; };
        std::vector<FaceData> fd(nClassify);

        // Same-winding exact duplicate suppression (kept from the previous
        // implementation; duplicates poison edge adjacency).
        auto hashDouble = [](double x) -> std::size_t {
            std::uint64_t bits = 0;
            std::memcpy(&bits, &x, sizeof(bits));
            bits ^= bits >> 33;
            bits *= 0xff51afd7ed558ccdULL;
            bits ^= bits >> 33;
            return static_cast<std::size_t>(bits);
        };
        auto hashPt = [&](const Vec& v) -> std::size_t {
            std::size_t h = hashDouble(v.x);
            h ^= hashDouble(v.y) + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= hashDouble(v.z) + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        };
        auto canonicalTriHash = [&](const Vec& pa, const Vec& pb, const Vec& pc) -> std::size_t {
            const std::size_t ha = hashPt(pa), hb = hashPt(pb), hc = hashPt(pc);
            std::size_t h0, h1, h2;
            if (ha <= hb && ha <= hc) { h0 = ha; h1 = hb; h2 = hc; }
            else if (hb <= ha && hb <= hc) { h0 = hb; h1 = hc; h2 = ha; }
            else { h0 = hc; h1 = ha; h2 = hb; }
            std::size_t h = h0;
            h ^= h1 + 0x9e3779b9u + (h << 6) + (h >> 2);
            h ^= h2 + 0x9e3779b9u + (h << 6) + (h >> 2);
            return h;
        };
        std::unordered_map<std::size_t, std::vector<uint32_t>> triHashMap;
        triHashMap.reserve(nClassify + 16);

        uint64_t iteration = 0;
        for (uint32_t i = 0; i < nClassify; i++)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("clip prepass");
            }

            Face tri = mesh.GetFace(i);
            FaceData& f = fd[i];
            f.a = mesh.GetPoint(tri.i0);
            f.b = mesh.GetPoint(tri.i1);
            f.c = mesh.GetPoint(tri.i2);
            f.pId = tri.pId;
            Vec cr = glm::cross(f.b - f.a, f.c - f.a);
            double crLen = glm::length(cr);
            f.area = 0.5 * crLen;
            f.n = crLen > EPS_NONZERO ? cr / crLen : Vec(0);
            f.center = isUnion ? (f.a + f.b * 2.0 + f.c * 3.0) / 6.0
                               : (f.a + f.b * 1.02 + f.c * 1.03) / 3.05;

            if (crLen <= EPS_NONZERO)
            {
                state[i] = ST_PRE_DROP;
                continue;
            }

            if (!isUnion)
            {
                auto aabb = mesh.GetFaceBox(i);
                if (!aabb.intersects(bvh2.box))
                {
                    // A-side face that cannot be cut by B: survives verbatim.
                    state[i] = ST_PRE_KEEP;
                    continue;
                }
                if (!aabb.intersects(bvh1.box))
                {
                    // B-only face: never part of A-B.
                    state[i] = ST_PRE_DROP;
                    continue;
                }

                const std::size_t triKey = canonicalTriHash(f.a, f.b, f.c);
                auto it = triHashMap.find(triKey);
                bool isDup = false;
                if (it != triHashMap.end())
                {
                    for (uint32_t prevIdx : it->second)
                    {
                        const FaceData& g = fd[prevIdx];
                        if ((equals(g.a, f.a, EPS_MINISCULE) && equals(g.b, f.b, EPS_MINISCULE) && equals(g.c, f.c, EPS_MINISCULE))
                            || (equals(g.a, f.b, EPS_MINISCULE) && equals(g.b, f.c, EPS_MINISCULE) && equals(g.c, f.a, EPS_MINISCULE))
                            || (equals(g.a, f.c, EPS_MINISCULE) && equals(g.b, f.a, EPS_MINISCULE) && equals(g.c, f.b, EPS_MINISCULE)))
                        {
                            isDup = true;
                            break;
                        }
                    }
                }
                if (isDup)
                {
                    state[i] = ST_PRE_DROP;
                    continue;
                }
                triHashMap[triKey].push_back(i);
            }
        }

        // ---- canonical vertices (tight grid; fragments of one plane share
        // exact CDT vertices) ----------------------------------------------
        const double vtxTol = weldTol;
        using GridKey = std::tuple<int64_t, int64_t, int64_t>;
        struct GridKeyHash {
            size_t operator()(const GridKey& k) const {
                size_t h = std::hash<int64_t>()(std::get<0>(k));
                h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
                h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
                return h;
            }
        };
        std::unordered_map<GridKey, std::vector<std::pair<uint32_t, Vec>>, GridKeyHash> vtxGrid;
        std::vector<Vec> canonicalPos;
        auto findOrAdd = [&](const Vec& p) -> uint32_t {
            GridKey center{ (int64_t)std::floor(p.x / vtxTol), (int64_t)std::floor(p.y / vtxTol), (int64_t)std::floor(p.z / vtxTol) };
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        GridKey nk{ std::get<0>(center) + dx, std::get<1>(center) + dy, std::get<2>(center) + dz };
                        auto it = vtxGrid.find(nk);
                        if (it == vtxGrid.end()) continue;
                        for (auto& [id, pos] : it->second)
                        {
                            Vec d = p - pos;
                            if (glm::dot(d, d) < vtxTol * vtxTol) return id;
                        }
                    }
            uint32_t id = (uint32_t)canonicalPos.size();
            canonicalPos.push_back(p);
            vtxGrid[center].emplace_back(id, p);
            return id;
        };

        std::vector<std::array<uint32_t, 3>> fvid(nClassify);
        for (uint32_t i = 0; i < nClassify; i++)
        {
            if (state[i] != ST_REGION) continue;
            fvid[i] = { findOrAdd(fd[i].a), findOrAdd(fd[i].b), findOrAdd(fd[i].c) };
        }

        // ---- vertex-on-other-plane sets for cut-edge detection -------------
        // A vertex "lies on" a plane of the other operand when its distance to
        // that plane is within the normalize tolerance. An edge whose BOTH
        // endpoints lie on the same other-operand plane is (part of) an
        // intersection curve: region status may flip across it.
        const bool haveProvenance = planesA != UINT32_MAX && planesA <= mesh.planes.size();
        std::unordered_set<uint64_t> vertexOnPlane;
        if (!mesh.planes.empty())
        {
            for (uint32_t vid = 0; vid < canonicalPos.size(); vid++)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("clip vertex-plane sets");
                }
                const Vec& p = canonicalPos[vid];
                for (uint32_t pi = 0; pi < mesh.planes.size(); pi++)
                {
                    double dist = std::fabs(glm::dot(mesh.planes[pi].normal, p) - mesh.planes[pi].distance);
                    if (dist <= cutDetectTol)
                    {
                        vertexOnPlane.insert((uint64_t)vid << 20 | (uint64_t)pi);
                    }
                }
            }
        }
        auto edgeIsCut = [&](uint32_t va, uint32_t vb, uint32_t ownPId) -> bool {
            const bool ownIsA = haveProvenance ? (ownPId < planesA) : false;
            const uint32_t beginP = haveProvenance ? (ownIsA ? planesA : 0u) : 0u;
            const uint32_t endP = haveProvenance ? (ownIsA ? (uint32_t)mesh.planes.size() : planesA) : (uint32_t)mesh.planes.size();
            for (uint32_t pi = beginP; pi < endP; pi++)
            {
                if (!haveProvenance && pi == ownPId) continue;
                if (vertexOnPlane.count((uint64_t)va << 20 | pi) && vertexOnPlane.count((uint64_t)vb << 20 | pi))
                {
                    return true;
                }
            }
            return false;
        };

        // ---- adjacency within plane groups over non-cut edges --------------
        struct EdgeKey {
            uint64_t k;
            bool operator==(const EdgeKey& o) const { return k == o.k; }
        };
        struct EdgeKeyHash { size_t operator()(const EdgeKey& e) const { return std::hash<uint64_t>()(e.k); } };
        auto makeEdgeKey = [](uint32_t pId, uint32_t va, uint32_t vb) -> EdgeKey {
            if (va > vb) std::swap(va, vb);
            // pId in the top bits, vertex pair below (20 bits each, ample for
            // per-operation mesh sizes).
            return EdgeKey{ ((uint64_t)pId << 40) ^ ((uint64_t)va << 20) ^ (uint64_t)vb };
        };
        std::unordered_map<EdgeKey, std::vector<uint32_t>, EdgeKeyHash> edgeFaces;
        edgeFaces.reserve((size_t)nClassify * 3 / 2);
        for (uint32_t i = 0; i < nClassify; i++)
        {
            if (state[i] != ST_REGION) continue;
            const auto& v = fvid[i];
            if (v[0] == v[1] || v[1] == v[2] || v[2] == v[0])
            {
                state[i] = ST_PRE_DROP;
                continue;
            }
            for (int e = 0; e < 3; e++)
            {
                edgeFaces[makeEdgeKey(fd[i].pId, v[e], v[(e + 1) % 3])].push_back(i);
            }
        }

        std::vector<std::vector<uint32_t>> adj(nClassify);
        for (auto& [key, owners] : edgeFaces)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("clip adjacency");
            }
            // Conservative region boundaries: only flood across clean manifold
            // edges that are not on an intersection curve.
            if (owners.size() != 2) continue;
            const uint32_t fa = owners[0], fb = owners[1];
            // shared edge vertices: recover from either face
            uint32_t va = 0, vb = 0;
            {
                const auto& v = fvid[fa];
                const auto& w = fvid[fb];
                int found = 0;
                for (int e1 = 0; e1 < 3 && found < 2; e1++)
                {
                    for (int e2 = 0; e2 < 3; e2++)
                    {
                        if (v[e1] == w[e2]) { (found == 0 ? va : vb) = v[e1]; found++; break; }
                    }
                }
                if (found < 2) continue;
            }
            if (edgeIsCut(va, vb, fd[fa].pId)) continue;
            adj[fa].push_back(fb);
            adj[fb].push_back(fa);
        }

        // ---- flood fill ----------------------------------------------------
        std::vector<int32_t> regionId(nClassify, -1);
        std::vector<std::vector<uint32_t>> regions;
        for (uint32_t i = 0; i < nClassify; i++)
        {
            if (state[i] != ST_REGION || regionId[i] >= 0) continue;
            const int32_t rid = (int32_t)regions.size();
            regions.emplace_back();
            std::queue<uint32_t> q;
            q.push(i);
            regionId[i] = rid;
            while (!q.empty())
            {
                uint32_t cur = q.front(); q.pop();
                regions[rid].push_back(cur);
                for (uint32_t nb : adj[cur])
                {
                    if (regionId[nb] < 0 && state[nb] == ST_REGION)
                    {
                        regionId[nb] = rid;
                        q.push(nb);
                    }
                }
            }
        }

        // ---- classify per region (majority over the largest fragments) -----
        std::vector<std::pair<size_t, bool>> sweepKillCandidates; // (regionIdx, isInteriorSheet)
        for (size_t regionIdx = 0; regionIdx < regions.size(); ++regionIdx)
        {
            auto& region = regions[regionIdx];
            if ((iteration++ & 63ULL) == 0)
            {
                budget.CheckDeadline("clip region classify");
                budget.CheckFaceCount(result.numFaces, "clip region classify");
            }

            std::sort(region.begin(), region.end(), [&](uint32_t x, uint32_t y) { return fd[x].area > fd[y].area; });

            // NOTE: a structural provenance-boundary vote assembly (own side
            // forced to BOUNDARY from pId provenance, cross side decided by a
            // plane-registry + footprint probe) was implemented and measured
            // here during the test3 open-edge investigation: it did not change
            // the test3 result and regressed the strict suite scoreboard while
            // adding per-region plane scans. Removed; details and numbers in
            // csg_test/v5/v5_work_notes.md.
            ClipDecision decision = ClipDecision::DROP;
            ClipDecision first = ClipDecision::DROP;
            bool decided = false;
            static const bool dumpRegions = std::getenv("CXDIAG_CSG_REGIONS") != nullptr;
            ClipVotes lastVotes;
            for (size_t s = 0; s < region.size() && s < 3; s++)
            {
                const FaceData& f = fd[region[s]];
                ClipVotes votes = castClipVotes(f.center, f.n, bvh1, bvh2, isUnion);
                lastVotes = votes;
                ClipDecision d = isUnion ? decideUnion(votes, f.n) : decideSubtract(votes, f.n);
                if (dumpRegions)
                {
                    double regionArea = 0.0;
                    for (uint32_t fi : region) regionArea += fd[fi].area;
                    if (regionArea > 0.005)
                    {
                        std::cerr << "DIAG: REGION faces=" << region.size() << " area=" << regionArea
                                  << " sample=" << s
                                  << " rep=(" << f.center.x << "," << f.center.y << "," << f.center.z << ")"
                                  << " n=(" << f.n.x << "," << f.n.y << "," << f.n.z << ")"
                                  << " loc1=" << (int)votes.loc1 << " loc2=" << (int)votes.loc2
                                  << " n1=(" << votes.normal1.x << "," << votes.normal1.y << "," << votes.normal1.z << ")"
                                  << " n2=(" << votes.normal2.x << "," << votes.normal2.y << "," << votes.normal2.z << ")"
                                  << " d=" << (int)d << std::endl;
                    }
                }
                if (s == 0)
                {
                    first = d;
                }
                else if (d == first)
                {
                    decision = d;
                    decided = true;
                    break;
                }
                else if (s == 2)
                {
                    // third sample disagrees with the first two disagreeing:
                    // majority impossible, trust the largest fragment.
                    decision = first;
                    decided = true;
                }
            }
            if (!decided)
            {
                decision = first;
            }
            if (decision == ClipDecision::DROP) continue;

            // Material sweep (kill-only): a correct result face has
            // result-material on exactly ONE side. Probe both sides of the
            // region representative against the OPERAND BVHs:
            //   subtract: in material := inside A and not inside B
            //   union:    in material := inside A or  inside B
            // Material on NEITHER side = single-layer membrane; material on
            // BOTH sides = interior sheet; both dropped. Any BOUNDARY probe =
            // uncertain = keep. NOTE: a material-PRIMARY variant (probes also
            // deciding keeps/rescues/winding) was implemented and measured
            // WORSE on every suite total (v5_strict2: +400 open, +40 nm, +68
            // membranes) -- probe parity at +-5e-4 is noisier than the
            // boundary-vote table near layered building surfaces. Kill-only
            // with the uncertainty guard is the validated configuration.
            {
                const FaceData& rep = fd[region[0]];
                // Probe distance defines the smallest material thickness the sweep treats as
                // real. 5*cutTol (~0.4mm at building scale) kept sub-millimetre residual
                // wedges (subtractor sill vs tilted layer face, test66/test61) as "material":
                // they export as open double-sheet shells the user sees as membranes.
                // diag*1e-4 matches the validator's membrane test, so wedges thinner than
                // the visible scale are swept with the rest of the zero-volume debris.
                // (A region-area-gated variant and a diag*5e-5 middle scale were measured:
                // both trade user-visible membranes on test61/test66 for open-edge counts on
                // already-failing files; this scale minimizes suite-wide membranes.)
                static const bool noBigProbe = std::getenv("CX_CSG_NO_BIGPROBE") != nullptr;
                const double probeEps = noBigProbe ? 5.0 * cutTol : std::max(5.0 * cutTol, meshDiag * 1.0e-4);
                bool uncertain = false;
                bool inMat[2];
                for (int side = 0; side < 2; side++)
                {
                    const Vec pt = rep.center + rep.n * (side == 0 ? probeEps : -probeEps);
                    auto lA = isInsideMesh(pt, rep.n, *bvh1.ptr, bvh1, rep.n, isUnion);
                    auto lB = isInsideMesh(pt, rep.n, *bvh2.ptr, bvh2, rep.n, isUnion);
                    if (lA.loc == MeshLocation::BOUNDARY || lB.loc == MeshLocation::BOUNDARY)
                    {
                        uncertain = true;
                        break;
                    }
                    const bool inA = lA.loc == MeshLocation::INSIDE;
                    const bool inB = lB.loc == MeshLocation::INSIDE;
                    inMat[side] = isUnion ? (inA || inB) : (inA && !inB);
                }
                if (!uncertain && inMat[0] == inMat[1])
                {
                    // Membrane (neither side) or interior sheet (both). Do NOT kill
                    // immediately: a sub-probe-thickness REAL slab (e.g. a 1mm layer
                    // remnant that is the ONLY cover of a rim band) reads as a
                    // membrane too, and killing it rips holes. Tentatively keep and
                    // defer to the gated sweep below. Interior sheets (material on
                    // BOTH sides) can never be visible sole-cover and always die.
                    sweepKillCandidates.emplace_back(regionIdx, inMat[0]);
                }
            }

            for (uint32_t fi : region)
            {
                const FaceData& f = fd[fi];
                // Flip decisions are relative to each fragment's own winding;
                // within a region all fragments share the plane orientation,
                // but guard against stray opposite-winding members by aligning
                // to the representative's normal.
                bool flip = decision == ClipDecision::KEEP_FLIP;
                if (glm::dot(f.n, fd[region[0]].n) < 0)
                {
                    flip = !flip;
                }
                state[fi] = flip ? ST_KEEP_FLIP : ST_KEEP;
            }
        }

        // ---- census-gated material sweep ------------------------------------
        // Apply the deferred membrane/interior kills, each gated on the kept
        // set's weld-level edge census: removing the region must not increase
        // the number of open edges. Debris sheets (their seams are already
        // open) pass the gate and die; a thin slab that is the sole cover of a
        // surface band fails it and survives (test57-allOpenings#: 1mm layer
        // remnants between per-layer opening rims).
        if (!sweepKillCandidates.empty())
        {
            std::unordered_map<uint64_t, int32_t> keptEdgeCount;
            keptEdgeCount.reserve((size_t)nClassify * 2);
            auto globalEdgeKey = [](uint32_t va, uint32_t vb) -> uint64_t {
                if (va > vb) std::swap(va, vb);
                return ((uint64_t)va << 32) | (uint64_t)vb;
            };
            for (uint32_t i = 0; i < nClassify; i++)
            {
                if (state[i] != ST_KEEP && state[i] != ST_KEEP_FLIP && state[i] != ST_PRE_KEEP) continue;
                const auto& v = fvid[i];
                for (int e = 0; e < 3; e++)
                {
                    if (v[e] == v[(e + 1) % 3]) continue;
                    keptEdgeCount[globalEdgeKey(v[e], v[(e + 1) % 3])]++;
                }
            }
            for (auto& [regionIdx, isInteriorSheet] : sweepKillCandidates)
            {
                auto& region = regions[regionIdx];
                bool anyKept = false;
                for (uint32_t fi : region)
                {
                    if (state[fi] == ST_KEEP || state[fi] == ST_KEEP_FLIP) { anyKept = true; break; }
                }
                if (!anyKept) continue;
                std::unordered_set<uint32_t> inRegion(region.begin(), region.end());
                int64_t deltaOpen = 0;
                for (uint32_t fi : region)
                {
                    if (state[fi] != ST_KEEP && state[fi] != ST_KEEP_FLIP) continue;
                    const auto& v = fvid[fi];
                    for (int e = 0; e < 3; e++)
                    {
                        if (v[e] == v[(e + 1) % 3]) continue;
                        auto it = keptEdgeCount.find(globalEdgeKey(v[e], v[(e + 1) % 3]));
                        if (it == keptEdgeCount.end()) continue;
                        if (it->second == 1) deltaOpen -= 1;      // region-owned open edge disappears
                        else if (it->second == 2)
                        {
                            // pair edge: if the partner is outside the region it becomes open
                            auto ownersIt = edgeFaces.find(makeEdgeKey(fd[fi].pId, v[e], v[(e + 1) % 3]));
                            bool partnerInRegion = false;
                            if (ownersIt != edgeFaces.end())
                            {
                                for (uint32_t owner : ownersIt->second)
                                {
                                    if (owner != fi && inRegion.count(owner)) { partnerInRegion = true; break; }
                                }
                            }
                            if (!partnerInRegion) deltaOpen += 1;
                        }
                        // count > 2: non-manifold either way, no open-edge change
                    }
                }
                bool kill = isInteriorSheet || deltaOpen <= 0;
                if (!kill && region.size() <= 64)
                {
                    // Killing leaks open edges -- still allowed when the region is
                    // REDUNDANT double-cover: every sampled point lies on kept surface
                    // of the same geometric plane owned by other regions, so removing
                    // it leaves surface behind (debris sheets welded onto real faces).
                    // A sole-cover slab fails this and survives.
                    static const bool noBigProbe2 = std::getenv("CX_CSG_NO_BIGPROBE") != nullptr;
                    const double coverTol = noBigProbe2 ? 5.0 * cutTol : std::max(5.0 * cutTol, meshDiag * 1.0e-4);
                    bool redundant = true;
                    for (uint32_t fi : region)
                    {
                        if (!redundant) break;
                        if (state[fi] != ST_KEEP && state[fi] != ST_KEEP_FLIP) continue;
                        const FaceData& f = fd[fi];
                        const Vec samples[4] = { f.center,
                                                 (f.center + f.a) * 0.5,
                                                 (f.center + f.b) * 0.5,
                                                 (f.center + f.c) * 0.5 };
                        for (const Vec& sample : samples)
                        {
                            bool covered = false;
                            for (uint32_t j = 0; j < nClassify && !covered; j++)
                            {
                                if (state[j] != ST_KEEP && state[j] != ST_KEEP_FLIP && state[j] != ST_PRE_KEEP) continue;
                                if (inRegion.count(j)) continue;
                                const FaceData& g = fd[j];
                                if (std::fabs(glm::dot(f.n, g.n)) < 0.98) continue;
                                if (std::fabs(glm::dot(g.n, sample - g.a)) > coverTol) continue;
                                // 2D containment on the dominant axis of g's normal
                                double gx = std::fabs(g.n.x), gy = std::fabs(g.n.y), gz = std::fabs(g.n.z);
                                int drop = (gz >= gx && gz >= gy) ? 2 : (gy >= gx ? 1 : 0);
                                auto p2 = [&](const Vec& p, double& u, double& v) {
                                    switch (drop)
                                    {
                                    case 0: u = p.y; v = p.z; break;
                                    case 1: u = p.x; v = p.z; break;
                                    default: u = p.x; v = p.y; break;
                                    }
                                };
                                double su, sv, au, av, bu, bv, cu, cv;
                                p2(sample, su, sv); p2(g.a, au, av); p2(g.b, bu, bv); p2(g.c, cu, cv);
                                double d1 = (su - bu) * (av - bv) - (au - bu) * (sv - bv);
                                double d2 = (su - cu) * (bv - cv) - (bu - cu) * (sv - cv);
                                double d3 = (su - au) * (cv - av) - (cu - au) * (sv - av);
                                bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
                                bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);
                                covered = !(hasNeg && hasPos);
                            }
                            if (!covered) { redundant = false; break; }
                        }
                    }
                    kill = redundant;
                    if (!kill)
                    {
                        // Final tier: re-probe at the SMALL scale (5*cutTol). A genuine
                        // zero-volume flap reads membrane/interior at every scale -- kill
                        // it and accept the leaked opens (the pre-gate behaviour, which the
                        // validator prefers over a surviving membrane). A real
                        // sub-probe-thickness slab reads material-on-one-side here and
                        // survives as the band's sole cover.
                        const FaceData& rep = fd[region[0]];
                        const double smallEps = 5.0 * cutTol;
                        bool uncertainSmall = false;
                        bool inMatSmall[2] = { false, false };
                        for (int side = 0; side < 2; side++)
                        {
                            const Vec pt = rep.center + rep.n * (side == 0 ? smallEps : -smallEps);
                            auto lA = isInsideMesh(pt, rep.n, *bvh1.ptr, bvh1, rep.n, isUnion);
                            auto lB = isInsideMesh(pt, rep.n, *bvh2.ptr, bvh2, rep.n, isUnion);
                            if (lA.loc == MeshLocation::BOUNDARY || lB.loc == MeshLocation::BOUNDARY)
                            {
                                uncertainSmall = true;
                                break;
                            }
                            const bool inA = lA.loc == MeshLocation::INSIDE;
                            const bool inB = lB.loc == MeshLocation::INSIDE;
                            inMatSmall[side] = isUnion ? (inA || inB) : (inA && !inB);
                        }
                        kill = !uncertainSmall && inMatSmall[0] == inMatSmall[1];
                    }
                }
                if (kill)
                {
                    for (uint32_t fi : region)
                    {
                        if (state[fi] == ST_KEEP || state[fi] == ST_KEEP_FLIP)
                        {
                            const auto& v = fvid[fi];
                            for (int e = 0; e < 3; e++)
                            {
                                if (v[e] == v[(e + 1) % 3]) continue;
                                auto it = keptEdgeCount.find(globalEdgeKey(v[e], v[(e + 1) % 3]));
                                if (it != keptEdgeCount.end() && it->second > 0) it->second--;
                            }
                            state[fi] = ST_PRE_DROP;
                        }
                    }
                }
            }
        }

        // NOTE: a "thin-strip rescue" (dropped small region adjacent to kept same-plane
        // faces adopts the neighbours' decision) was implemented and measured here in two
        // variants (unconditional and BOUNDARY-ambiguity-gated): both were net negative
        // (suite 26 -> 25, interior sheets up, test53b/53d sharply worse) because the
        // rescue cannot distinguish rim strips that belong to the surface from strips
        // that were correctly cut away. The strip-vote problem needs provenance-aware
        // voting (which rim pair the strip lies between), not adjacency heuristics.
        // Run labels fable5_rescue / fable5_rescue2 in csg_test/v5.

        // ---- same-plane overlap resolution ----------------------------------
        // The normalized mesh is not a perfect partition: verbatim re-added
        // faces and the CDTs of near-coplanar planes can overlap within the
        // same geometric plane. The old per-fragment classifier dropped such
        // fragments by ray luck; region classification keeps them
        // consistently, which turns each overlap into non-manifold edges
        // (observed on test44: micro CDT fragments overlapping kept faces).
        // Resolve deterministically: within one geometric plane, when two kept
        // fragments of the same effective orientation overlap (the smaller
        // fragment's interior sample lies inside the larger one), drop the
        // smaller fragment.
        {
            auto planeKeyOf = [&](const FaceData& f, bool flipped) -> uint64_t {
                Vec n = flipped ? -f.n : f.n;
                double ax = std::fabs(n.x), ay = std::fabs(n.y), az = std::fabs(n.z);
                double dom = (az >= ax && az >= ay) ? n.z : (ay >= ax ? n.y : n.x);
                int sign = 1;
                if (dom < 0.0) { n = -n; sign = -1; }
                double d = glm::dot(n, f.a);
                // quantize at normalize tolerance; orientation sign in bit 0
                int64_t qx = (int64_t)std::llround(n.x * 1000.0);
                int64_t qy = (int64_t)std::llround(n.y * 1000.0);
                int64_t qz = (int64_t)std::llround(n.z * 1000.0);
                int64_t qd = (int64_t)std::llround(d / cutTol);
                uint64_t h = (uint64_t)(qx * 73856093LL ^ qy * 19349663LL ^ qz * 83492791LL ^ qd * 2654435761LL);
                return (h << 1) | (uint64_t)(sign > 0 ? 1 : 0);
            };
            std::unordered_map<uint64_t, std::vector<uint32_t>> keptByPlane;
            for (uint32_t i = 0; i < nClassify; i++)
            {
                if (state[i] != ST_KEEP && state[i] != ST_KEEP_FLIP) continue;
                keptByPlane[planeKeyOf(fd[i], state[i] == ST_KEEP_FLIP)].push_back(i);
            }
            auto pointInTri2D = [](double px, double py, double ax, double ay, double bx, double by, double cx, double cy) -> bool {
                double d1 = (px - bx) * (ay - by) - (ax - bx) * (py - by);
                double d2 = (px - cx) * (by - cy) - (bx - cx) * (py - cy);
                double d3 = (px - ax) * (cy - ay) - (cx - ax) * (py - ay);
                bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
                bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);
                return !(hasNeg && hasPos);
            };
            for (auto& [key, faces] : keptByPlane)
            {
                if (faces.size() < 2) continue;
                if ((iteration++ & 63ULL) == 0)
                {
                    budget.CheckDeadline("clip overlap resolve");
                }
                // project on dominant axis of the group normal
                const Vec n0 = fd[faces[0]].n;
                double ax = std::fabs(n0.x), ay = std::fabs(n0.y), az = std::fabs(n0.z);
                int dropAxis = (az >= ax && az >= ay) ? 2 : (ay >= ax ? 1 : 0);
                auto proj = [&](const Vec& p, double& u, double& v) {
                    switch (dropAxis)
                    {
                    case 0: u = p.y; v = p.z; break;
                    case 1: u = p.x; v = p.z; break;
                    default: u = p.x; v = p.y; break;
                    }
                };
                std::sort(faces.begin(), faces.end(), [&](uint32_t x, uint32_t y) { return fd[x].area > fd[y].area; });
                for (size_t s = 1; s < faces.size(); s++)
                {
                    if (state[faces[s]] == ST_PRE_DROP) continue;
                    const FaceData& frag = fd[faces[s]];
                    // A fragment is overlap-redundant only when FULLY covered
                    // by a larger kept fragment: all three vertices and the
                    // centroid inside (boundary counts as inside). Dropping
                    // partial overlaps removes real surface and manufactures
                    // holes (observed on test53d with a centroid-only rule).
                    double pu, pv, p0u, p0v, p1u, p1v, p2u, p2v;
                    proj((frag.a + frag.b + frag.c) / 3.0, pu, pv);
                    proj(frag.a, p0u, p0v);
                    proj(frag.b, p1u, p1v);
                    proj(frag.c, p2u, p2v);
                    for (size_t l = 0; l < s; l++)
                    {
                        if (state[faces[l]] == ST_PRE_DROP) continue;
                        const FaceData& big = fd[faces[l]];
                        double au, av, bu, bv, cu, cv;
                        proj(big.a, au, av);
                        proj(big.b, bu, bv);
                        proj(big.c, cu, cv);
                        if (pointInTri2D(pu, pv, au, av, bu, bv, cu, cv)
                            && pointInTri2D(p0u, p0v, au, av, bu, bv, cu, cv)
                            && pointInTri2D(p1u, p1v, au, av, bu, bv, cu, cv)
                            && pointInTri2D(p2u, p2v, au, av, bu, bv, cu, cv))
                        {
                            // Which of the two is garbage? A fragment that was
                            // emitted FLIPPED only got there through winding
                            // inversion -- weaker evidence than a fragment kept
                            // in its native orientation. When a flipped and a
                            // native fragment overlap on the same oriented
                            // plane, the flipped one is the misclassification
                            // (observed on test43##: a large flipped wedge
                            // overlapping the true top surface), regardless of
                            // size. Between same-kind fragments the smaller
                            // loses.
                            const bool bigFlipped = state[faces[l]] == ST_KEEP_FLIP;
                            const bool fragFlipped = state[faces[s]] == ST_KEEP_FLIP;
                            if (bigFlipped && !fragFlipped)
                            {
                                state[faces[l]] = ST_PRE_DROP;
                                // the big one is gone; re-test this fragment
                                // against the remaining larger fragments
                                continue;
                            }
                            state[faces[s]] = ST_PRE_DROP;
                            break;
                        }
                    }

                    // Union coverage (second chance): a strip fragment overlapping SEVERAL
                    // larger kept fragments is contained in none of them individually, so the
                    // single-coverer test above keeps it and the overlap surfaces as
                    // non-manifold triples along the strip's edges (test53d: bottom-edge
                    // strips, nm=3 edge census). Sample the fragment (vertices, edge
                    // midpoints, centroid): when every sample lies in SOME larger kept
                    // fragment of the same oriented plane, the strip is redundant surface and
                    // the smaller fragment loses. Any uncovered sample keeps it.
                    static const bool noUnionCov = std::getenv("CX_CSG_NO_UNIONCOV") != nullptr;
                    if (!noUnionCov && state[faces[s]] != ST_PRE_DROP)
                    {
                        const double sampleU[7] = { pu, p0u, p1u, p2u, (p0u + p1u) / 2.0, (p1u + p2u) / 2.0, (p2u + p0u) / 2.0 };
                        const double sampleV[7] = { pv, p0v, p1v, p2v, (p0v + p1v) / 2.0, (p1v + p2v) / 2.0, (p2v + p0v) / 2.0 };
                        bool allCovered = true;
                        for (int sp = 0; sp < 7 && allCovered; sp++)
                        {
                            bool covered = false;
                            for (size_t l = 0; l < s; l++)
                            {
                                if (state[faces[l]] == ST_PRE_DROP) continue;
                                const FaceData& big = fd[faces[l]];
                                double au, av, bu, bv, cu, cv;
                                proj(big.a, au, av);
                                proj(big.b, bu, bv);
                                proj(big.c, cu, cv);
                                if (pointInTri2D(sampleU[sp], sampleV[sp], au, av, bu, bv, cu, cv))
                                {
                                    covered = true;
                                    break;
                                }
                            }
                            allCovered = covered;
                        }
                        if (allCovered)
                        {
                            state[faces[s]] = ST_PRE_DROP;
                        }
                    }
                }
            }
        }


        // ---- emission --------------------------------------------------------
        // Track every emitted classified triangle as an UNDIRECTED canonical
        // triple. The verbatim tail (irrelevant re-adds) can carry the original
        // copy of a face whose re-triangulated twin was already emitted by
        // classification with corrected winding; emitting both creates a
        // reversed-duplicate pair that downstream cleanup rightly removes IN
        // FULL -- which deleted every cavity cap of test3 (54 caps x 2 copies)
        // and left 162 open edges. The classified copy wins; the tail copy is
        // skipped.
        struct UTriKey {
            uint64_t k1, k2;
            bool operator==(const UTriKey& o) const { return k1 == o.k1 && k2 == o.k2; }
        };
        struct UTriKeyHash { size_t operator()(const UTriKey& t) const { return std::hash<uint64_t>()(t.k1) ^ (std::hash<uint64_t>()(t.k2) * 0x9e3779b97f4a7c15ULL); } };
        std::unordered_set<UTriKey, UTriKeyHash> emittedTriples;
        emittedTriples.reserve(nClassify);
        auto makeUTri = [&](uint32_t va, uint32_t vb, uint32_t vc) -> UTriKey {
            uint32_t s0 = va, s1 = vb, s2 = vc;
            if (s0 > s1) std::swap(s0, s1);
            if (s1 > s2) std::swap(s1, s2);
            if (s0 > s1) std::swap(s0, s1);
            return UTriKey{ ((uint64_t)s0 << 32) | s1, (uint64_t)s2 };
        };

        for (uint32_t i = 0; i < nClassify; i++)
        {
            const FaceData& f = fd[i];
            if (state[i] == ST_KEEP || state[i] == ST_PRE_KEEP)
            {
                result.AddFace(f.a, f.b, f.c, f.pId);
            }
            else if (state[i] == ST_KEEP_FLIP)
            {
                result.AddFace(f.b, f.a, f.c, f.pId);
            }
            else
            {
                continue;
            }
            // canonicalize via the weld grid (pre-keep faces have no fvid)
            emittedTriples.insert(makeUTri(findOrAdd(f.a), findOrAdd(f.b), findOrAdd(f.c)));
        }

        for (uint32_t i = mesh.data; i < mesh.numFaces; i++)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("clip tail loop");
                budget.CheckFaceCount(result.numFaces, "clip tail loop");
            }
            Face tri = mesh.GetFace(i);
            const Vec a = mesh.GetPoint(tri.i0);
            const Vec b = mesh.GetPoint(tri.i1);
            const Vec c = mesh.GetPoint(tri.i2);
            if (!emittedTriples.empty())
            {
                const UTriKey key = makeUTri(findOrAdd(a), findOrAdd(b), findOrAdd(c));
                if (emittedTriples.count(key))
                {
                    continue; // already emitted by classification (correct winding)
                }
            }
            result.AddFace(a, b, c, tri.pId);
        }

        budget.CheckFaceCount(result.numFaces, isUnion ? "clipJoin complete" : "clipSubtract complete");
    }

    // ------------------------------------------------------------------
    // T-junction stitch (post-clip).
    //
    // Normalize re-triangulates every shared plane independently, so two
    // faces meeting across a shared boundary edge can disagree about the
    // subdivision points placed on that edge. The clipped result then
    // LOOKS closed but carries hairline open-edge pairs wherever a
    // subdivision vertex of one side lies in the interior of the other
    // side's edge -- and because boolean chains feed each result back
    // into Normalize, the damage compounds with every operation
    // (measured on test42: 42 open edges after step 2, 506 after step
    // 6). Repair at the source: for every open edge, find welded result
    // vertices lying on its interior and fan-split the owning triangle
    // so both sides agree about the subdivision again.
    // ------------------------------------------------------------------
    // Cheap open/non-manifold edge census on welded vertices, used to verify
    // that a repair pass actually improved the mesh before keeping it.
    static std::pair<uint32_t, uint32_t> CountOpenAndNonManifoldEdges(const Geometry& g, double weldTol)
    {
        std::pair<uint32_t, uint32_t> res{ 0, 0 };
        if (g.numFaces == 0) return res;
        const double weldTolSq = weldTol * weldTol;
        using GridKey = std::tuple<int64_t, int64_t, int64_t>;
        struct GridKeyHash {
            size_t operator()(const GridKey& k) const {
                size_t h = std::hash<int64_t>()(std::get<0>(k));
                h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9u + (h << 6) + (h >> 2);
                h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9u + (h << 6) + (h >> 2);
                return h;
            }
        };
        std::unordered_map<GridKey, std::vector<std::pair<uint32_t, Vec>>, GridKeyHash> grid;
        uint32_t nextId = 0;
        auto findOrAdd = [&](const Vec& p) -> uint32_t {
            GridKey c{ (int64_t)std::floor(p.x / weldTol), (int64_t)std::floor(p.y / weldTol), (int64_t)std::floor(p.z / weldTol) };
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        GridKey nk{ std::get<0>(c) + dx, std::get<1>(c) + dy, std::get<2>(c) + dz };
                        auto it = grid.find(nk);
                        if (it == grid.end()) continue;
                        for (auto& [id, pos] : it->second)
                        {
                            Vec d = p - pos;
                            if (glm::dot(d, d) < weldTolSq) return id;
                        }
                    }
            uint32_t id = nextId++;
            grid[c].emplace_back(id, p);
            return id;
        };
        std::unordered_map<uint64_t, uint32_t> edgeCount;
        edgeCount.reserve((size_t)g.numFaces * 3);
        for (uint32_t i = 0; i < g.numFaces; i++)
        {
            Face tri = g.GetFace(i);
            uint32_t v[3] = { findOrAdd(g.GetPoint(tri.i0)), findOrAdd(g.GetPoint(tri.i1)), findOrAdd(g.GetPoint(tri.i2)) };
            for (int e = 0; e < 3; e++)
            {
                uint32_t a = v[e], b = v[(e + 1) % 3];
                if (a == b) continue;
                if (a > b) std::swap(a, b);
                ++edgeCount[((uint64_t)a << 32) | b];
            }
        }
        for (auto& kv : edgeCount)
        {
            if (kv.second == 1) ++res.first;
            else if (kv.second > 2) ++res.second;
        }
        return res;
    }

    // Both-sides-outside membrane census (the artifact the viewer shows as
    // sheets across openings). Rays are cast along +/- face normal with a
    // slight jitter -- a ray with dir . n > 0 from the +n side can never
    // re-cross the face plane, so excluding the tested face (required for
    // true membranes, which are not volume boundaries) cannot delete a
    // genuine crossing. Same algorithm as the fixed export pass 6.
    static uint32_t CountMembraneArtifactFaces(Geometry& g)
    {
        if (g.numFaces == 0 || g.numFaces > 20000) return 0;
        const AABB box = g.GetAABB();
        const double diag = glm::length(box.max - box.min);
        if (!(diag > 0.0)) return 0;
        const double off = std::max(diag * 1e-5, 1e-6);
        BVH bvh = MakeBVH(g);
        auto sideOutside = [&](const Vec& c, const Vec& n, double s, uint32_t skip) -> bool {
            Vec axis = std::abs(n.x) < 0.9 ? Vec(1, 0, 0) : Vec(0, 1, 0);
            Vec u = glm::normalize(glm::cross(n, axis));
            Vec w = glm::cross(n, u);
            const Vec base = n * s;
            const Vec dirs[2] = {
                glm::normalize(base + u * 0.11 + w * 0.07),
                glm::normalize(base - u * 0.13 + w * 0.17)
            };
            const Vec p = c + base * off;
            int oddCount = 0;
            for (const Vec& dir : dirs)
            {
                int crossings = 0;
                std::unordered_set<uint32_t> seen;
                bvh.IntersectRay(p, dir, [&](uint32_t i) -> bool {
                    if (i == skip) return false;
                    if (!seen.insert(i).second) return false;
                    Face f2 = g.GetFace(i);
                    Vec a2 = g.GetPoint(f2.i0);
                    Vec b2 = g.GetPoint(f2.i1);
                    Vec c2 = g.GetPoint(f2.i2);
                    Vec hit;
                    double dist, dPlane;
                    if (intersect_ray_triangle(p, p + dir, a2, b2, c2, hit, dist, dPlane, true))
                    {
                        // Ignore crossings inside the offset near-field: a
                        // coincident double-cover twin of the tested face
                        // (stitch fan overlap, kernel double-emission) sits at
                        // distance ~off and would flip the parity of BOTH
                        // layers, making legitimate geometry score as
                        // membranes. Real second surfaces are feature-scale
                        // away. Export pass 4b removes the twins themselves.
                        if (dist > off * 2.5)
                        {
                            ++crossings;
                        }
                    }
                    return false;
                });
                if ((crossings & 1) != 0) ++oddCount;
            }
            // 0 = both rays even (outside), 2 = both odd (inside), 1 = unknown
            return oddCount;
        };
        uint32_t cnt = 0;
        for (uint32_t i = 0; i < g.numFaces; i++)
        {
            Face f = g.GetFace(i);
            Vec a = g.GetPoint(f.i0);
            Vec b = g.GetPoint(f.i1);
            Vec c = g.GetPoint(f.i2);
            Vec cr = glm::cross(b - a, c - a);
            double len = glm::length(cr);
            if (len < 1e-15) continue;
            Vec n = cr / len;
            Vec ctr = (a + b + c) / 3.0;
            const int sPlus = sideOutside(ctr, n, 1.0, i);
            const int sMinus = sideOutside(ctr, n, -1.0, i);
            // A boundary face has material on exactly one side. Both sides
            // outside = floating membrane; both sides inside = buried interior
            // sheet (the test107 regression class). Mixed parity = unknown,
            // benefit of the doubt.
            if ((sPlus == 0 && sMinus == 0) || (sPlus == 2 && sMinus == 2)) ++cnt;
        }
        return cnt;
    }

    // Combined mesh quality score for choosing between candidate boolean
    // results: open edges weigh 1, non-manifold edges 3 (harder to repair),
    // membranes 10 (the user-visible artifact class the edge census cannot
    // see). Lower is better.
    static uint64_t MeshQualityScore(Geometry& g)
    {
        const AABB box = g.GetAABB();
        const double diag = glm::length(box.max - box.min);
        const double weldTol = std::max(1.0E-07, diag * 1e-9);
        const auto census = CountOpenAndNonManifoldEdges(g, weldTol);
        const uint32_t membranes = CountMembraneArtifactFaces(g);
        return (uint64_t)census.first + 3ULL * census.second + 10ULL * membranes;
    }

    static void StitchTJunctions(Geometry& g, const BooleanBudget& budget)
    {
        static const bool stitchOff = std::getenv("CSG_NOSTITCH") != nullptr;
        if (stitchOff) return;
        if (g.numFaces == 0 || g.numFaces > 20000) return;
        const AABB box = g.GetAABB();
        const double diag = glm::length(box.max - box.min);
        if (!(diag > 0.0)) return;
        const double weldTol = std::max(1.0E-07, diag * 1e-9);
        // Snapshot + verify: the fan splits are only kept when they strictly
        // improved the census (fewer open edges, no new non-manifold edges).
        // On meshes whose opens are NOT T-junctions the splits can land on
        // coincident sheets and mint non-manifold edges (measured: test58,
        // test60, test61); reverting per-mesh keeps the stitch strictly
        // no-harm.
        const auto censusBefore = CountOpenAndNonManifoldEdges(g, weldTol);
        if (censusBefore.first == 0) return;
        const Geometry snapshot = g;
        const double weldTolSq = weldTol * weldTol;
        const double stitchTol = std::max(_TOLERANCE_PLANE_DEVIATION, diag * 2e-6);
        const double stitchTolSq = stitchTol * stitchTol;

        using GridKey = std::tuple<int64_t, int64_t, int64_t>;
        struct GridKeyHash {
            size_t operator()(const GridKey& k) const {
                size_t h = std::hash<int64_t>()(std::get<0>(k));
                h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9u + (h << 6) + (h >> 2);
                h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9u + (h << 6) + (h >> 2);
                return h;
            }
        };
        struct EdgeKey {
            uint32_t v0, v1;
            bool operator==(const EdgeKey& o) const { return v0 == o.v0 && v1 == o.v1; }
        };
        struct EdgeKeyHash {
            size_t operator()(const EdgeKey& k) const {
                size_t h = std::hash<uint32_t>()(k.v0);
                h ^= std::hash<uint32_t>()(k.v1) + 0x9e3779b9u + (h << 6) + (h >> 2);
                return h;
            }
        };

        for (int iter = 0; iter < 6; ++iter)
        {
            budget.CheckDeadline("stitch iter");
            const uint32_t nF = g.numFaces;

            std::unordered_map<GridKey, std::vector<std::pair<uint32_t, Vec>>, GridKeyHash> grid;
            std::vector<Vec> repPos;
            repPos.reserve(nF * 3);
            auto cellOf = [&](const Vec& p) -> GridKey {
                return { (int64_t)std::floor(p.x / weldTol),
                        (int64_t)std::floor(p.y / weldTol),
                        (int64_t)std::floor(p.z / weldTol) };
            };
            auto findOrAdd = [&](const Vec& p) -> uint32_t {
                auto c = cellOf(p);
                for (int dx = -1; dx <= 1; ++dx)
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dz = -1; dz <= 1; ++dz)
                        {
                            GridKey nk{ std::get<0>(c) + dx, std::get<1>(c) + dy, std::get<2>(c) + dz };
                            auto it = grid.find(nk);
                            if (it == grid.end()) continue;
                            for (auto& [id, pos] : it->second)
                            {
                                Vec d = p - pos;
                                if (glm::dot(d, d) < weldTolSq) return id;
                            }
                        }
                uint32_t id = (uint32_t)repPos.size();
                repPos.push_back(p);
                grid[c].emplace_back(id, p);
                return id;
            };

            std::vector<std::array<uint32_t, 3>> fv(nF);
            std::vector<std::array<Vec, 3>> fpos(nF);
            std::vector<uint32_t> fpid(nF);
            for (uint32_t i = 0; i < nF; i++)
            {
                Face tri = g.GetFace(i);
                Vec a = g.GetPoint(tri.i0);
                Vec b = g.GetPoint(tri.i1);
                Vec c = g.GetPoint(tri.i2);
                fv[i] = { findOrAdd(a), findOrAdd(b), findOrAdd(c) };
                fpos[i] = { a, b, c };
                fpid[i] = tri.pId;
            }

            // Undirected edge census; remember the single owner of each edge
            // seen exactly once.
            struct EdgeUse { uint32_t count = 0; uint32_t face = 0; uint8_t slot = 0; };
            std::unordered_map<EdgeKey, EdgeUse, EdgeKeyHash> edges;
            edges.reserve(nF * 3);
            for (uint32_t i = 0; i < nF; i++)
            {
                for (int e = 0; e < 3; e++)
                {
                    uint32_t a = fv[i][e], b = fv[i][(e + 1) % 3];
                    if (a == b) continue;
                    EdgeKey k{ std::min(a, b), std::max(a, b) };
                    auto& u = edges[k];
                    if (u.count == 0) { u.face = i; u.slot = (uint8_t)e; }
                    ++u.count;
                }
            }

            std::vector<std::pair<uint32_t, uint8_t>> openOwners;
            for (auto& kv : edges)
            {
                if (kv.second.count == 1) openOwners.emplace_back(kv.second.face, kv.second.slot);
            }
            if (openOwners.empty()) break;
            if ((uint64_t)openOwners.size() * repPos.size() > 80000000ULL) break;

            // One split per face per iteration; the fan introduces fresh edges
            // that the next iteration re-examines.
            std::vector<uint8_t> faceDirty(nF, 0);
            struct Fan { uint32_t face; uint8_t slot; std::vector<std::pair<double, uint32_t>> mids; };
            std::vector<Fan> fans;
            for (auto& [fi, slot] : openOwners)
            {
                if (faceDirty[fi]) continue;
                const Vec A = fpos[fi][slot];
                const Vec B = fpos[fi][(slot + 1) % 3];
                const Vec AB = B - A;
                const double len2 = glm::dot(AB, AB);
                if (len2 < weldTolSq * 4.0) continue;
                const double len = std::sqrt(len2);
                const double tEps = std::max(2.0 * weldTol / len, 1e-9);
                const uint32_t ia = fv[fi][slot], ib = fv[fi][(slot + 1) % 3];

                std::vector<std::pair<double, uint32_t>> mids;
                for (uint32_t v = 0; v < (uint32_t)repPos.size(); ++v)
                {
                    if (v == ia || v == ib) continue;
                    const Vec& p = repPos[v];
                    const double t = glm::dot(p - A, AB) / len2;
                    if (t <= tEps || t >= 1.0 - tEps) continue;
                    const Vec foot = A + AB * t;
                    const Vec d = p - foot;
                    if (glm::dot(d, d) >= stitchTolSq) continue;
                    mids.emplace_back(t, v);
                }
                if (mids.empty()) continue;
                std::sort(mids.begin(), mids.end());
                faceDirty[fi] = 1;
                fans.push_back(Fan{ fi, slot, std::move(mids) });
            }
            if (fans.empty()) break;

            Geometry ng;
            ng.planes = g.planes;
            ng.hasPlanes = g.hasPlanes;
            ng.mBoolOpCount = g.mBoolOpCount;
            for (uint32_t i = 0; i < nF; i++)
            {
                if (faceDirty[i]) continue;
                ng.AddFace(fpos[i][0], fpos[i][1], fpos[i][2], fpid[i]);
            }
            for (auto& fan : fans)
            {
                const Vec A = fpos[fan.face][fan.slot];
                const Vec B = fpos[fan.face][(fan.slot + 1) % 3];
                const Vec C = fpos[fan.face][(fan.slot + 2) % 3];
                Vec prev = A;
                for (auto& [t, v] : fan.mids)
                {
                    ng.AddFace(prev, repPos[v], C, fpid[fan.face]);
                    prev = repPos[v];
                }
                ng.AddFace(prev, B, C, fpid[fan.face]);
            }
            g = std::move(ng);
        }

        // Keep the stitched mesh only if the full quality score (open edges,
        // non-manifold edges, membranes) improved. Single-signal criteria were
        // measured and rejected: census-only guards either starve chains whose
        // intermediates need closing (test53d) or admit membrane debris the
        // edge census cannot see (test63/64).
        Geometry stitched = std::move(g);
        g = snapshot;
        const uint64_t scoreBefore = MeshQualityScore(g);
        const uint64_t scoreAfter = MeshQualityScore(stitched);
        if (scoreAfter < scoreBefore)
        {
            g = std::move(stitched);
        }
    }

    static Geometry clipJoin(Geometry& mesh, BVH bvh1, BVH bvh2, const BooleanBudget& budget, uint32_t planesA = UINT32_MAX)
    {
        Geometry resultingMesh;
        clipMeshRegions(mesh, bvh1, bvh2, resultingMesh, budget, true, planesA);
        StitchTJunctions(resultingMesh, budget);
        return resultingMesh;
    }

    static Geometry clipJoin(Geometry& mesh, BVH bvh1, BVH bvh2)
    {
        auto budget = BooleanBudget::Unlimited(mesh.numFaces);
        return clipJoin(mesh, bvh1, bvh2, budget);
    }

    static Geometry clipSubtract(Geometry& mesh, BVH bvh1, BVH bvh2, const BooleanBudget& budget, uint32_t planesA = UINT32_MAX)
    {
        Geometry resultingMesh;
        clipMeshRegions(mesh, bvh1, bvh2, resultingMesh, budget, false, planesA);
        StitchTJunctions(resultingMesh, budget);
        return resultingMesh;
    }

    static Geometry clipSubtract(Geometry& mesh, BVH bvh1, BVH bvh2)
    {
        auto budget = BooleanBudget::Unlimited(mesh.numFaces);
        return clipSubtract(mesh, bvh1, bvh2, budget);
    }
}
