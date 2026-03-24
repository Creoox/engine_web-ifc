#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"

#include <queue>
#include <array>
#include <tuple>
#include <functional>
#include <algorithm>

namespace fuzzybools
{
	inline void SetEpsilons(double TOLERANCE_PLANE_INTERSECTION, double TOLERANCE_PLANE_DEVIATION, double TOLERANCE_BACK_DEVIATION_DISTANCE, double TOLERANCE_INSIDE_OUTSIDE_PERIMETER, double TOLERANCE_BOUNDING_BOX, double BOOLSTATUS)
	{
		_TOLERANCE_PLANE_INTERSECTION = TOLERANCE_PLANE_INTERSECTION;
		_TOLERANCE_PLANE_DEVIATION = TOLERANCE_PLANE_DEVIATION;
		_TOLERANCE_BACK_DEVIATION_DISTANCE = TOLERANCE_BACK_DEVIATION_DISTANCE;
		_TOLERANCE_INSIDE_OUTSIDE_PERIMETER = TOLERANCE_INSIDE_OUTSIDE_PERIMETER;
		_TOLERANCE_BOUNDING_BOX = TOLERANCE_BOUNDING_BOX;
		_BOOLSTATUS = BOOLSTATUS;
	}

	// ---------------------------------------------------------------
	// Post-boolean cleanup.
	//
	// Phase A -- strip near-degenerate sliver triangles (area < 1e-9).
	//
	// Phase B -- remove duplicate faces at junction edges.
	//   For each edge shared by 3+ faces, find coplanar face pairs
	//   and remove the one that is a membrane (identified by a
	//   ray-based backing test: a solid face has another face behind
	//   it within the wall thickness, a membrane does not).
	//
	// Phase C -- remove disconnected non-manifold open-shell
	//   components with near-zero volume/area (orphan fragments).
	// ---------------------------------------------------------------
	inline void CleanNonManifoldShells(Geometry &result)
	{
		// -- Phase A: strip sliver triangles -----------------------------
		{
			constexpr double SLIVER_AREA_THRESHOLD = 1e-9;
			const uint32_t n = result.numFaces;
			uint32_t sliverCount = 0;
			for (uint32_t i = 0; i < n; i++)
			{
				Face f = result.GetFace(i);
				if (areaOfTriangle(result.GetPoint(f.i0),
				                   result.GetPoint(f.i1),
				                   result.GetPoint(f.i2)) < SLIVER_AREA_THRESHOLD)
					sliverCount++;
			}
			if (sliverCount > 0)
			{
				Geometry tmp;
				tmp.planes = result.planes;
				tmp.hasPlanes = result.hasPlanes;
				tmp.data = result.data;
				for (uint32_t i = 0; i < n; i++)
				{
					Face f = result.GetFace(i);
					Vec a = result.GetPoint(f.i0);
					Vec b = result.GetPoint(f.i1);
					Vec c = result.GetPoint(f.i2);
					if (areaOfTriangle(a, b, c) >= SLIVER_AREA_THRESHOLD)
						tmp.AddFace(a, b, c, f.pId);
				}
				result = tmp;
			}
		}

		// -- Shared setup for Phases B & C --------------------------------
		const uint32_t nFaces = result.numFaces;
		if (nFaces < 4) return;

		// -- Step 1: cache face info -------------------------------------
		struct FV {
			Vec a, b, c;
			uint32_t pId;
			Vec center, normal;
			double area, maxEdge;
		};
		std::vector<FV> fv(nFaces);
		for (uint32_t i = 0; i < nFaces; i++)
		{
			Face f = result.GetFace(i);
			Vec a = result.GetPoint(f.i0);
			Vec b = result.GetPoint(f.i1);
			Vec c = result.GetPoint(f.i2);
			Vec crossP = glm::cross(b - a, c - a);
			double crossLen = glm::length(crossP);
			double e0 = glm::length(b - a);
			double e1 = glm::length(c - b);
			double e2 = glm::length(a - c);
			fv[i] = {a, b, c,
			         static_cast<uint32_t>(f.pId),
			         (a + b + c) / 3.0,
			         crossLen > 1e-15 ? crossP / crossLen : Vec(0),
			         crossLen * 0.5,
			         std::max(e0, std::max(e1, e2))};
		}

		// -- Step 2: spatial-hash vertex deduplication -------------------
		const double cellSize = toleranceVectorEquality * 1.5;
		const double cellSizeSq = cellSize * cellSize;

		using GridKey = std::tuple<int64_t, int64_t, int64_t>;
		struct GridKeyHash {
			size_t operator()(const GridKey &k) const {
				size_t h = std::hash<int64_t>()(std::get<0>(k));
				h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
				h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
				return h;
			}
		};

		std::unordered_map<GridKey, std::vector<std::pair<uint32_t, Vec>>, GridKeyHash> vtxGrid;
		std::vector<std::array<uint32_t, 3>> fvid(nFaces);
		uint32_t nextVid = 0;

		auto getCell = [&](const Vec &p) -> GridKey {
			return {static_cast<int64_t>(std::floor(p.x / cellSize)),
			        static_cast<int64_t>(std::floor(p.y / cellSize)),
			        static_cast<int64_t>(std::floor(p.z / cellSize))};
		};

		auto findOrAdd = [&](const Vec &p) -> uint32_t {
			auto center = getCell(p);
			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz)
					{
						GridKey nk = {std::get<0>(center) + dx,
						              std::get<1>(center) + dy,
						              std::get<2>(center) + dz};
						auto it = vtxGrid.find(nk);
						if (it != vtxGrid.end())
						{
							for (auto &[id, pos] : it->second)
							{
								Vec d = p - pos;
								if (glm::dot(d, d) < cellSizeSq)
									return id;
							}
						}
					}
			uint32_t id = nextVid++;
			vtxGrid[center].emplace_back(id, p);
			return id;
		};

		for (uint32_t i = 0; i < nFaces; i++)
		{
			fvid[i][0] = findOrAdd(fv[i].a);
			fvid[i][1] = findOrAdd(fv[i].b);
			fvid[i][2] = findOrAdd(fv[i].c);
		}

		// -- Step 3: edge -> face adjacency -----------------------------
		struct EKey {
			uint32_t v0, v1;
			bool operator==(const EKey &o) const { return v0 == o.v0 && v1 == o.v1; }
		};
		struct EKeyHash {
			size_t operator()(const EKey &k) const {
				size_t h = std::hash<uint32_t>()(k.v0);
				h ^= std::hash<uint32_t>()(k.v1) + 0x9e3779b9 + (h << 6) + (h >> 2);
				return h;
			}
		};
		std::unordered_map<EKey, std::vector<uint32_t>, EKeyHash> edgeFaces;
		for (uint32_t i = 0; i < nFaces; i++)
			for (int e = 0; e < 3; e++)
			{
				uint32_t va = fvid[i][e];
				uint32_t vb = fvid[i][(e + 1) % 3];
				EKey ek = {std::min(va, vb), std::max(va, vb)};
				edgeFaces[ek].push_back(i);
			}

		// -- Step 4: face adjacency -------------------------------------
		std::vector<std::vector<uint32_t>> faceAdj(nFaces);
		for (auto &[ek, fl] : edgeFaces)
			for (size_t a = 0; a < fl.size(); a++)
				for (size_t b = a + 1; b < fl.size(); b++)
				{
					faceAdj[fl[a]].push_back(fl[b]);
					faceAdj[fl[b]].push_back(fl[a]);
				}

		// ================================================================
		// Phase B: thin membrane detection & removal
		// ================================================================
		constexpr double THIN_THRESHOLD = 1e-3;      // 1 mm
		constexpr double MEMBRANE_THICKNESS = 5e-3;   // 5 mm
		std::vector<bool> thinMarked(nFaces, false);
		BVH resultBVH = MakeBVH(result);

		// B.1: double-layer detection (within 1 mm)
		for (uint32_t i = 0; i < nFaces; i++)
		{
			if (thinMarked[i]) continue;
			if (glm::length(fv[i].normal) < 0.5) continue;
			for (int sign = -1; sign <= 1; sign += 2)
			{
				if (thinMarked[i]) break;
				Vec rayDir = fv[i].normal * static_cast<double>(sign);
				Vec rayEnd = fv[i].center + rayDir * THIN_THRESHOLD;
				resultBVH.IntersectRay(fv[i].center, rayDir,
					[&](uint32_t fj) -> bool
				{
					if (fj == i) return false;
					if (glm::dot(fv[i].normal, fv[fj].normal) > -0.7)
						return false;
					Vec hitPos; double t, dp;
					if (intersect_ray_triangle(fv[i].center, rayEnd,
						fv[fj].a, fv[fj].b, fv[fj].c, hitPos, t, dp, false))
					{
						double dist = glm::length(hitPos - fv[i].center);
						if (dist > 1e-6 && dist < THIN_THRESHOLD)
						{ thinMarked[i] = true; thinMarked[fj] = true; return true; }
					}
					return false;
				});
			}
		}

		// B.2: Manifold-edge component analysis.
		//   Build adjacency using ONLY manifold edges (count == 2).
		//   Junction edges (3+) separate the solid from membranes.
		//   Components with near-zero thickness are membranes.
		//   Structural edges (perpendicular to junction neighbor) are kept.
		{
			std::vector<std::vector<uint32_t>> manifoldAdj(nFaces);
			for (auto &[ek2, fl2] : edgeFaces)
			{
				if (fl2.size() != 2) continue;
				manifoldAdj[fl2[0]].push_back(fl2[1]);
				manifoldAdj[fl2[1]].push_back(fl2[0]);
			}
			std::vector<int> mCompId(nFaces, -1);
			int mNumComp = 0;
			for (uint32_t i = 0; i < nFaces; i++)
			{
				if (mCompId[i] >= 0) continue;
				int cid = mNumComp++;
				std::queue<uint32_t> q;
				q.push(i); mCompId[i] = cid;
				while (!q.empty())
				{
					uint32_t cur = q.front(); q.pop();
					for (uint32_t nb : manifoldAdj[cur])
						if (mCompId[nb] < 0) { mCompId[nb] = cid; q.push(nb); }
				}
			}
			if (mNumComp > 1)
			{
				struct MCI { uint32_t count = 0; Vec centroid{0}; double area = 0; };
				std::vector<MCI> mci(mNumComp);
				for (uint32_t i = 0; i < nFaces; i++)
				{
					mci[mCompId[i]].count++;
					mci[mCompId[i]].centroid += fv[i].center;
					mci[mCompId[i]].area += fv[i].area;
				}
				for (int c = 0; c < mNumComp; c++)
					if (mci[c].count > 0) mci[c].centroid /= (double)mci[c].count;
				std::vector<double> mVol(mNumComp, 0.0);
				for (uint32_t i = 0; i < nFaces; i++)
				{
					int ci = mCompId[i];
					Vec va = fv[i].a - mci[ci].centroid;
					Vec vb = fv[i].b - mci[ci].centroid;
					Vec vc = fv[i].c - mci[ci].centroid;
					mVol[ci] += glm::dot(va, glm::cross(vb, vc)) / 6.0;
				}
				int largestComp = 0;
				for (int c = 1; c < mNumComp; c++)
					if (std::abs(mVol[c]) > std::abs(mVol[largestComp]))
						largestComp = c;
				// Junction-normal check: structural edges are perpendicular
				std::vector<double> compMaxJuncDot(mNumComp, -1.0);
				for (auto &[ek2, fl2] : edgeFaces)
				{
					if (fl2.size() < 3) continue;
					for (size_t ja = 0; ja < fl2.size(); ja++)
						for (size_t jb = ja + 1; jb < fl2.size(); jb++)
						{
							int ca = mCompId[fl2[ja]], cb = mCompId[fl2[jb]];
							if (ca == cb) continue;
							double d = std::abs(glm::dot(
								fv[fl2[ja]].normal, fv[fl2[jb]].normal));
							compMaxJuncDot[ca] = std::max(compMaxJuncDot[ca], d);
							compMaxJuncDot[cb] = std::max(compMaxJuncDot[cb], d);
						}
				}
				for (uint32_t i = 0; i < nFaces; i++)
				{
					int ci = mCompId[i];
					if (ci == largestComp) continue;
					double thickness = mci[ci].area > 1e-12
						? std::abs(mVol[ci]) / mci[ci].area : 0;
					if (thickness < MEMBRANE_THICKNESS && compMaxJuncDot[ci] >= 0.5)
						thinMarked[i] = true;
				}
			}
		}

		// B.3: Junction-edge duplicate removal.
		//   For each junction edge with 3+ faces, find coplanar pairs.
		//   The face with fewer manifold-edge connections is the membrane.
		{
			std::vector<uint8_t> manifoldEdgeCnt(nFaces, 0);
			for (auto &[ek2, fl2] : edgeFaces)
				if (fl2.size() == 2)
					for (uint32_t fi : fl2)
						manifoldEdgeCnt[fi]++;
			for (auto &[ek2, fl2] : edgeFaces)
			{
				if (fl2.size() < 3) continue;
				for (size_t a = 0; a < fl2.size(); a++)
				{
					if (thinMarked[fl2[a]]) continue;
					for (size_t b = a + 1; b < fl2.size(); b++)
					{
						if (thinMarked[fl2[b]]) continue;
						double dot = glm::dot(fv[fl2[a]].normal, fv[fl2[b]].normal);
						if (std::abs(dot) < 0.7) continue;
						uint32_t fa = fl2[a], fb = fl2[b];
						if (manifoldEdgeCnt[fa] > manifoldEdgeCnt[fb])
							thinMarked[fb] = true;
						else if (manifoldEdgeCnt[fb] > manifoldEdgeCnt[fa])
							thinMarked[fa] = true;
					}
				}
			}
		}

		// Safety: if > 75% marked, disable phase B
		uint32_t thinCount = 0;
		for (uint32_t i = 0; i < nFaces; i++)
			if (thinMarked[i]) thinCount++;
		bool thinEnabled = (thinCount > 0 && thinCount < nFaces * 3 / 4);

		// ================================================================
		// Phase C: non-manifold zero-volume component removal
		//   Uses centroid-relative volume / area (thickness) to avoid
		//   the origin-relative volume bug.
		// ================================================================
		std::vector<int> compId(nFaces, -1);
		int numComp = 0;
		for (uint32_t i = 0; i < nFaces; i++)
		{
			if (compId[i] >= 0) continue;
			int cid = numComp++;
			std::queue<uint32_t> q;
			q.push(i);
			compId[i] = cid;
			while (!q.empty())
			{
				uint32_t cur = q.front(); q.pop();
				for (uint32_t nb : faceAdj[cur])
					if (compId[nb] < 0)
					{
						compId[nb] = cid;
						q.push(nb);
					}
			}
		}

		struct CompInfo {
			uint32_t faceCount = 0;
			uint32_t boundaryEdges = 0;
			Vec centroid{0};
			double area = 0;
			double volume = 0;
		};
		std::vector<CompInfo> comps(numComp);

		for (uint32_t i = 0; i < nFaces; i++)
		{
			int cid = compId[i];
			comps[cid].faceCount++;
			comps[cid].centroid += fv[i].center;
			comps[cid].area += fv[i].area;
		}
		for (int cid = 0; cid < numComp; cid++)
			if (comps[cid].faceCount > 0)
				comps[cid].centroid /= static_cast<double>(comps[cid].faceCount);

		// Centroid-relative signed volume
		for (uint32_t i = 0; i < nFaces; i++)
		{
			int cid = compId[i];
			Vec va = fv[i].a - comps[cid].centroid;
			Vec vb = fv[i].b - comps[cid].centroid;
			Vec vc = fv[i].c - comps[cid].centroid;
			comps[cid].volume += glm::dot(va, glm::cross(vb, vc)) / 6.0;
		}
		for (auto &[ek, fl] : edgeFaces)
			if (fl.size() != 2)
				comps[compId[fl[0]]].boundaryEdges++;

		std::vector<bool> badComp(numComp, false);
		uint32_t removedByComp = 0;
		if (numComp > 1)
		{
			// Find largest-volume component (always kept)
			int largestComp = 0;
			for (int cid = 1; cid < numComp; cid++)
				if (std::abs(comps[cid].volume) > std::abs(comps[largestComp].volume))
					largestComp = cid;

			for (int cid = 0; cid < numComp; cid++)
			{
				if (cid == largestComp) continue;
				bool isManifold = (comps[cid].boundaryEdges == 0);
				double thickness = comps[cid].area > 1e-12
					? std::abs(comps[cid].volume) / comps[cid].area : 0;
				if (!isManifold && thickness < MEMBRANE_THICKNESS)
				{
					badComp[cid] = true;
					removedByComp += comps[cid].faceCount;
				}
			}
		}

		// ================================================================
		// Rebuild
		// ================================================================
		uint32_t totalRemoved = (thinEnabled ? thinCount : 0) + removedByComp;
		if (totalRemoved == 0 || totalRemoved >= nFaces) return;

		Geometry cleaned;
		cleaned.planes = result.planes;
		cleaned.hasPlanes = result.hasPlanes;
		for (uint32_t i = 0; i < nFaces; i++)
		{
			if (thinEnabled && thinMarked[i]) continue;
			if (badComp[compId[i]]) continue;
			cleaned.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
		}
		cleaned.data = result.data;
		result = cleaned;
	}

	// Returns true when the AABBs of A and B overlap by at least
	// minPenetration in every axis.  Operands that merely touch at a
	// surface (overlap ~ 0) return false -- the boolean would produce
	// no volume change but can create coplanar-face artifacts.
	inline bool HasVolumeOverlap(const Geometry &A, const Geometry &B,
	                             double minPenetration = 1e-3)
	{
		AABB a = A.GetAABB();
		AABB b = B.GetAABB();
		double ox = std::min(a.max.x, b.max.x) - std::max(a.min.x, b.min.x);
		double oy = std::min(a.max.y, b.max.y) - std::max(a.min.y, b.min.y);
		double oz = std::min(a.max.z, b.max.z) - std::max(a.min.z, b.min.z);
		return ox > minPenetration && oy > minPenetration && oz > minPenetration;
	}

	// Count boundary (non-manifold) edges in a geometry for diagnostics.
	inline uint32_t CountBoundaryEdges(const Geometry &g)
	{
		const double cs = toleranceVectorEquality * 1.5;
		const double cs2 = cs * cs;
		using GK = std::tuple<int64_t, int64_t, int64_t>;
		struct GKH { size_t operator()(const GK &k) const {
			size_t h = std::hash<int64_t>()(std::get<0>(k));
			h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h<<6) + (h>>2);
			h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h<<6) + (h>>2);
			return h; }};
		std::unordered_map<GK, std::vector<std::pair<uint32_t, Vec>>, GKH> grid;
		uint32_t nv = 0;
		auto cell = [&](const Vec &p) -> GK {
			return {(int64_t)std::floor(p.x/cs),(int64_t)std::floor(p.y/cs),(int64_t)std::floor(p.z/cs)};};
		auto foa = [&](const Vec &p) -> uint32_t {
			auto c0 = cell(p);
			for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++) for (int dz=-1;dz<=1;dz++){
				GK nk={std::get<0>(c0)+dx,std::get<1>(c0)+dy,std::get<2>(c0)+dz};
				auto it=grid.find(nk); if (it!=grid.end())
					for (auto &[id,pos]:it->second){ Vec d=p-pos; if (glm::dot(d,d)<cs2) return id; }}
			uint32_t id=nv++; grid[c0].emplace_back(id,p); return id; };
		struct EK2 { uint32_t v0,v1; bool operator==(const EK2&o)const{return v0==o.v0&&v1==o.v1;}};
		struct EKH2 { size_t operator()(const EK2&k)const{
			size_t h=std::hash<uint32_t>()(k.v0);
			h^=std::hash<uint32_t>()(k.v1)+0x9e3779b9+(h<<6)+(h>>2); return h;}};
		std::unordered_map<EK2,uint32_t,EKH2> ec;
		for (uint32_t i=0;i<g.numFaces;i++){
			Face f=g.GetFace(i);
			uint32_t ids[3]={foa(g.GetPoint(f.i0)),foa(g.GetPoint(f.i1)),foa(g.GetPoint(f.i2))};
			for (int e=0;e<3;e++){
				uint32_t a=ids[e],b=ids[(e+1)%3];
				EK2 ek={std::min(a,b),std::max(a,b)}; ec[ek]++;}}
		uint32_t boundary=0;
		for (auto &[ek,cnt]:ec) if (cnt!=2) boundary++;
		return boundary;
	}

	// ---------------------------------------------------------------
	// RepairBoundaryEdges -- after clipSubtract, some faces may have
	// been incorrectly rejected, creating boundary (open) edges.
	// This function identifies boundary edges in `result`, looks
	// for faces in the normalized mesh `normalized` (up to index
	// `nDataLimit`) that were NOT included in `result` but share
	// a boundary edge, and adds them back if doing so reduces the
	// total number of boundary edges.  Iterates until no improvement.
	// ---------------------------------------------------------------
	inline void RepairBoundaryEdges(Geometry &result,
	                                const Geometry &normalized,
	                                uint32_t nDataLimit)
	{
		const double cs  = toleranceVectorEquality * 1.5;
		const double cs2 = cs * cs;

		// -- vertex dedup helpers (reusable) ----
		using GK = std::tuple<int64_t, int64_t, int64_t>;
		struct GKH { size_t operator()(const GK &k) const {
			size_t h = std::hash<int64_t>()(std::get<0>(k));
			h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h<<6) + (h>>2);
			h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h<<6) + (h>>2);
			return h; }};
		struct EK { uint32_t v0, v1;
			bool operator==(const EK &o) const { return v0==o.v0 && v1==o.v1; }};
		struct EKH { size_t operator()(const EK &k) const {
			size_t h = std::hash<uint32_t>()(k.v0);
			h ^= std::hash<uint32_t>()(k.v1) + 0x9e3779b9 + (h<<6) + (h>>2);
			return h; }};

		// Cache the normalized mesh faces + deduplicated vertex IDs
		struct NF { Vec a, b, c; uint32_t pId; uint32_t vid[3]; };
		std::vector<NF> nf(nDataLimit);
		{
			std::unordered_map<GK, std::vector<std::pair<uint32_t,Vec>>, GKH> grid;
			uint32_t nv = 0;
			auto cell = [&](const Vec &p) -> GK {
				return {(int64_t)std::floor(p.x/cs),
				        (int64_t)std::floor(p.y/cs),
				        (int64_t)std::floor(p.z/cs)};};
			auto foa = [&](const Vec &p) -> uint32_t {
				auto c0 = cell(p);
				for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++)
				for (int dz=-1;dz<=1;dz++){
					GK nk={std::get<0>(c0)+dx,std::get<1>(c0)+dy,
					        std::get<2>(c0)+dz};
					auto it=grid.find(nk); if (it!=grid.end())
						for (auto &[id,pos]:it->second){
							Vec d=p-pos;
							if (glm::dot(d,d)<cs2) return id; }}
				uint32_t id=nv++; grid[c0].emplace_back(id,p);
				return id; };
			for (uint32_t i = 0; i < nDataLimit; i++)
			{
				Face f = normalized.GetFace(i);
				Vec a = normalized.GetPoint(f.i0);
				Vec b = normalized.GetPoint(f.i1);
				Vec c = normalized.GetPoint(f.i2);
				nf[i] = {a, b, c, static_cast<uint32_t>(f.pId),
				         {foa(a), foa(b), foa(c)}};
			}
		}

		// Build a set of "result face signatures" so we can quickly
		// check which normalized faces are already in the result.
		// Signature: sorted triple of deduplicated vertex IDs.
		struct TriSig {
			uint32_t v[3];
			bool operator==(const TriSig &o) const {
				return v[0]==o.v[0] && v[1]==o.v[1] && v[2]==o.v[2]; }
		};
		struct TriSigH { size_t operator()(const TriSig &s) const {
			size_t h = std::hash<uint32_t>()(s.v[0]);
			h ^= std::hash<uint32_t>()(s.v[1]) + 0x9e3779b9 + (h<<6) + (h>>2);
			h ^= std::hash<uint32_t>()(s.v[2]) + 0x9e3779b9 + (h<<6) + (h>>2);
			return h; }};

		auto makeTriSig = [](uint32_t a, uint32_t b, uint32_t c) -> TriSig {
			uint32_t v[3] = {a, b, c};
			if (v[0]>v[1]) std::swap(v[0],v[1]);
			if (v[1]>v[2]) std::swap(v[1],v[2]);
			if (v[0]>v[1]) std::swap(v[0],v[1]);
			return {v[0],v[1],v[2]};
		};

		// Build dedup grid that covers BOTH result and normalized verts
		// (re-do from scratch so IDs are consistent).
		std::unordered_map<GK, std::vector<std::pair<uint32_t,Vec>>, GKH> vGrid;
		uint32_t vNext = 0;
		auto vCell = [&](const Vec &p) -> GK {
			return {(int64_t)std::floor(p.x/cs),
			        (int64_t)std::floor(p.y/cs),
			        (int64_t)std::floor(p.z/cs)};};
		auto vFoa = [&](const Vec &p) -> uint32_t {
			auto c0 = vCell(p);
			for (int dx=-1;dx<=1;dx++) for (int dy=-1;dy<=1;dy++)
			for (int dz=-1;dz<=1;dz++){
				GK nk={std::get<0>(c0)+dx,std::get<1>(c0)+dy,
				        std::get<2>(c0)+dz};
				auto it=vGrid.find(nk); if (it!=vGrid.end())
					for (auto &[id,pos]:it->second){
						Vec d=p-pos;
						if (glm::dot(d,d)<cs2) return id; }}
			uint32_t id=vNext++; vGrid[c0].emplace_back(id,p);
			return id; };

		// Build edge counts and result face sigs from current result.
		// These are updated INCREMENTALLY as faces are added.
		vGrid.clear(); vNext = 0;
		std::unordered_map<EK, uint32_t, EKH> edgeCnt;
		std::unordered_set<TriSig, TriSigH> resultSigs;

		uint32_t rn = result.numFaces;
		for (uint32_t i = 0; i < rn; i++)
		{
			Face f = result.GetFace(i);
			uint32_t ids[3] = {vFoa(result.GetPoint(f.i0)),
			                   vFoa(result.GetPoint(f.i1)),
			                   vFoa(result.GetPoint(f.i2))};
			for (int e = 0; e < 3; e++)
			{
				EK ek = {std::min(ids[e], ids[(e+1)%3]),
				         std::max(ids[e], ids[(e+1)%3])};
				edgeCnt[ek]++;
			}
			resultSigs.insert(makeTriSig(ids[0], ids[1], ids[2]));
		}

		for (int iter = 0; iter < 100; iter++)
		{
			bool added = false;
			for (uint32_t i = 0; i < nDataLimit; i++)
			{
				uint32_t ids[3] = {vFoa(nf[i].a), vFoa(nf[i].b), vFoa(nf[i].c)};
				if (ids[0]==ids[1]||ids[1]==ids[2]||ids[0]==ids[2]) continue;
				TriSig ts = makeTriSig(ids[0], ids[1], ids[2]);
				if (resultSigs.count(ts)) continue;

				int delta = 0;
				bool anyBdry = false;
				EK faceEdges[3];
				for (int e = 0; e < 3; e++)
				{
					faceEdges[e] = {std::min(ids[e], ids[(e+1)%3]),
					                std::max(ids[e], ids[(e+1)%3])};
					auto it = edgeCnt.find(faceEdges[e]);
					if (it == edgeCnt.end())
						delta++;
					else if (it->second == 1)
						{ delta--; anyBdry = true; }
					else if (it->second == 2)
						delta++;
				}
				if (!anyBdry || delta >= 0) continue;

				result.AddFace(nf[i].a, nf[i].b, nf[i].c, nf[i].pId);
				resultSigs.insert(ts);
				for (int e = 0; e < 3; e++) edgeCnt[faceEdges[e]]++;
				added = true;
			}
			if (!added) break;
		}
	}

	inline Geometry Subtract(const Geometry &A, const Geometry &B)
	{
		if (!HasVolumeOverlap(A, B))
			return A;

		uint32_t bdryA = CountBoundaryEdges(A);

		fuzzybools::SharedPosition sp;
		sp.Construct(A, B, false);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(A, B, sp, false);

		auto result = fuzzybools::clipSubtract(geom, bvh1, bvh2);

		uint32_t bdryClip = CountBoundaryEdges(result);
		if (bdryClip > bdryA)
			RepairBoundaryEdges(result, geom, geom.data);

		CleanNonManifoldShells(result);
		return result;
	}

	inline Geometry Union(const Geometry &A, const Geometry &B)
	{
		if (!HasVolumeOverlap(A, B))
		{
			// No real overlap -- just merge both geometries.
			Geometry merged = A;
			for (uint32_t i = 0; i < B.numFaces; i++)
			{
				Face f = B.GetFace(i);
				merged.AddFace(B.GetPoint(f.i0), B.GetPoint(f.i1),
				               B.GetPoint(f.i2), f.pId);
			}
			return merged;
		}

		fuzzybools::SharedPosition sp;
		sp.Construct(A, B, true);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(A, B, sp, true);

		auto result = fuzzybools::clipJoin(geom, bvh1, bvh2);
		CleanNonManifoldShells(result);
		return result;
	}
}
