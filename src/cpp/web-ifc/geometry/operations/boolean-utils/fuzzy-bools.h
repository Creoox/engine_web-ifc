#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"

#include <queue>
#include <array>
#include <tuple>
#include <functional>

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
	// Post-boolean cleanup: remove non-manifold open-shell connected
	// components whose signed volume is near zero.
	//
	// After nested boolean operations, floating-point drift can leave
	// orphaned face fragments that are not part of any closed solid.
	// A valid CSG result must be a closed manifold; isolated open
	// shells with negligible volume are artifacts.
	//
	// Uses the same spatial-hash vertex deduplication as
	// SharedPosition::AddPoint (cell size = toleranceVectorEquality,
	// 27-cell neighbourhood, exact distance check) so that vertices
	// that were considered identical during the boolean operation are
	// still merged here.
	// ---------------------------------------------------------------
	inline void CleanNonManifoldShells(Geometry &result)
	{
		const uint32_t nFaces = result.numFaces;
		if (nFaces < 4) return; // need >= 4 faces for any closed solid

		// -- Step 1: cache face vertices --------------------------------
		struct FV { Vec a, b, c; uint32_t pId; };
		std::vector<FV> fv(nFaces);
		for (uint32_t i = 0; i < nFaces; i++)
		{
			Face f = result.GetFace(i);
			fv[i] = {result.GetPoint(f.i0),
			         result.GetPoint(f.i1),
			         result.GetPoint(f.i2),
			         static_cast<uint32_t>(f.pId)};
		}

		// -- Step 2: spatial-hash vertex deduplication -------------------
		// Mirror SharedPosition::AddPoint: grid cell = floor(coord/cellSize),
		// 27-cell neighbourhood search, exact Euclidean distance check.
		const double cellSize = toleranceVectorEquality;          // 1e-4
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

		// canonical vertex id per position
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
			// search 27 neighbours
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
			// not found — insert
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

		// -- Step 3: edge → face adjacency ------------------------------
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

		// -- Step 4: BFS connected components ---------------------------
		std::vector<std::vector<uint32_t>> faceAdj(nFaces);
		for (auto &[ek, fl] : edgeFaces)
			for (size_t a = 0; a < fl.size(); a++)
				for (size_t b = a + 1; b < fl.size(); b++)
				{
					faceAdj[fl[a]].push_back(fl[b]);
					faceAdj[fl[b]].push_back(fl[a]);
				}

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

		if (numComp <= 1) return; // single component — nothing to clean

		// -- Step 5: per-component manifold check + signed volume -------
		struct CompInfo {
			uint32_t faceCount = 0;
			uint32_t boundaryEdges = 0; // edges shared by != 2 faces
			double   volume = 0;
		};
		std::vector<CompInfo> comps(numComp);

		for (uint32_t i = 0; i < nFaces; i++)
		{
			int cid = compId[i];
			comps[cid].faceCount++;
			// Signed volume via tetrahedron-with-origin formula
			comps[cid].volume += glm::dot(fv[i].a, glm::cross(fv[i].b, fv[i].c)) / 6.0;
		}
		for (auto &[ek, fl] : edgeFaces)
			if (fl.size() != 2)
				comps[compId[fl[0]]].boundaryEdges++;

		// -- Step 6: discard non-manifold, near-zero-volume shells ------
		// 1e-6 m³ = 1 mm³ — well below any real architectural feature
		// but safely above numerical noise from flat face fragments.
		constexpr double VOLUME_THRESHOLD = 1e-6;
		uint32_t removedFaces = 0;
		for (int cid = 0; cid < numComp; cid++)
		{
			bool isManifold = (comps[cid].boundaryEdges == 0);
			bool hasVolume  = (std::abs(comps[cid].volume) >= VOLUME_THRESHOLD);
			if (!isManifold && !hasVolume)
				removedFaces += comps[cid].faceCount;
		}

		if (removedFaces == 0 || removedFaces >= nFaces) return;

		// Rebuild geometry keeping only valid components
		Geometry cleaned;
		cleaned.planes = result.planes;
		cleaned.hasPlanes = result.hasPlanes;
		for (uint32_t i = 0; i < nFaces; i++)
		{
			int cid = compId[i];
			bool isManifold = (comps[cid].boundaryEdges == 0);
			bool hasVolume  = (std::abs(comps[cid].volume) >= VOLUME_THRESHOLD);
			if (isManifold || hasVolume)
				cleaned.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
		}
		cleaned.data = result.data; // preserve face-count metadata
		result = cleaned;
	}

	inline Geometry Subtract(const Geometry &A, const Geometry &B)
	{
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B, false);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(A, B, sp, false);

#ifdef CSG_DEBUG_OUTPUT
//	DumpGeometry(geom, L"Post-normalize.obj");
#endif

		auto result = fuzzybools::clipSubtract(geom, bvh1, bvh2);
		CleanNonManifoldShells(result);
		return result;
	}

	inline Geometry Union(const Geometry &A, const Geometry &B)
	{
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