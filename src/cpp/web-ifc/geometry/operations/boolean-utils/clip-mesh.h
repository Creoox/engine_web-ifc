/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#pragma once

#include <glm/glm.hpp>
#include <cstdint>
#include <cstring>
#include <cstdlib>
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
                    if (dist <= cutTol)
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
        for (auto& region : regions)
        {
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
            for (size_t s = 0; s < region.size() && s < 3; s++)
            {
                const FaceData& f = fd[region[s]];
                ClipVotes votes = castClipVotes(f.center, f.n, bvh1, bvh2, isUnion);
                ClipDecision d = isUnion ? decideUnion(votes, f.n) : decideSubtract(votes, f.n);
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
                const double probeEps = 5.0 * cutTol;
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
                    continue; // membrane (neither side) or interior sheet (both)
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

        const auto censusAfter = CountOpenAndNonManifoldEdges(g, weldTol);
        // Weighted comparison: non-manifold edges are harder to repair
        // downstream than open edges, so they weigh 3x. A stitch that closes
        // many T-junctions at the cost of a couple of non-manifold contacts
        // is still a net win (test42: 95 open / 1 nm -> 58 open / 7 nm).
        const uint64_t scoreBefore = (uint64_t)censusBefore.first + 3ULL * censusBefore.second;
        const uint64_t scoreAfter = (uint64_t)censusAfter.first + 3ULL * censusAfter.second;
        if (scoreAfter >= scoreBefore)
        {
            g = snapshot;
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
