#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"

#include <queue>
#include <array>
#include <tuple>
#include <functional>
#include <algorithm>
#include <unordered_set>

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
	// Shared mesh-edge data produced by BuildEdgeMap.
	// Holds spatial-hash vertex deduplication results, per-face
	// canonical vertex IDs, the edge-to-face map, and the count of
	// boundary (open) edges.  Passed by reference so the caller can
	// reuse the data without rebuilding.
	// ---------------------------------------------------------------
	struct EdgeMapData {
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

		std::unordered_map<uint32_t, Vec> vidPos;                              // canonical vertex ID -> position
		std::vector<std::array<uint32_t, 3>> fvid;                             // per-face canonical vertex IDs
		std::vector<uint32_t> fpid;                                            // per-face plane ID
		std::unordered_map<EKey, std::vector<uint32_t>, EKeyHash> edgeFaces;   // edge -> face indices
		uint32_t openEdgeCount = 0;
	};

	// ---------------------------------------------------------------
	// Build (or rebuild) the edge map for a Geometry.
	// Uses spatial-hash vertex deduplication identical to
	// SharedPosition::AddPoint (cell = toleranceVectorEquality * 1.5,
	// 27-cell neighbourhood, exact Euclidean distance check).
	// ---------------------------------------------------------------
	inline void BuildEdgeMap(const Geometry &geo, EdgeMapData &out)
	{
		out = EdgeMapData{};
		const uint32_t nf = geo.numFaces;
		if (nf < 3) return;

		const double cs = toleranceVectorEquality * 1.5;
		const double csSq = cs * cs;

		using GK = std::tuple<int64_t, int64_t, int64_t>;
		struct GKH {
			size_t operator()(const GK &k) const {
				size_t h = std::hash<int64_t>()(std::get<0>(k));
				h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
				h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
				return h;
			}
		};

		std::unordered_map<GK, std::vector<std::pair<uint32_t, Vec>>, GKH> grid;
		uint32_t nextId = 0;

		auto cell = [&](const Vec &p) -> GK {
			return {static_cast<int64_t>(std::floor(p.x / cs)),
			        static_cast<int64_t>(std::floor(p.y / cs)),
			        static_cast<int64_t>(std::floor(p.z / cs))};
		};

		auto findOrAdd = [&](const Vec &p) -> uint32_t {
			auto c = cell(p);
			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz)
					{
						GK nk = {std::get<0>(c) + dx,
						          std::get<1>(c) + dy,
						          std::get<2>(c) + dz};
						auto it = grid.find(nk);
						if (it != grid.end())
							for (auto &[id, pos] : it->second)
							{
								Vec d = p - pos;
								if (glm::dot(d, d) < csSq) return id;
							}
					}
			uint32_t id = nextId++;
			grid[c].emplace_back(id, p);
			out.vidPos[id] = p;
			return id;
		};

		out.fvid.resize(nf);
		out.fpid.resize(nf);
		for (uint32_t i = 0; i < nf; i++)
		{
			Face f = geo.GetFace(i);
			out.fvid[i][0] = findOrAdd(geo.GetPoint(f.i0));
			out.fvid[i][1] = findOrAdd(geo.GetPoint(f.i1));
			out.fvid[i][2] = findOrAdd(geo.GetPoint(f.i2));
			out.fpid[i] = static_cast<uint32_t>(f.pId);
			for (int e = 0; e < 3; e++)
			{
				uint32_t a = out.fvid[i][e], b = out.fvid[i][(e + 1) % 3];
				EdgeMapData::EKey ek = {std::min(a, b), std::max(a, b)};
				out.edgeFaces[ek].push_back(i);
			}
		}

		for (auto &[ek, fl] : out.edgeFaces)
			if (fl.size() == 1) out.openEdgeCount++;
	}

	inline bool meshSanityCheck(Geometry& result) {
		bool ok = true;

		// indexData must hold exactly 3 indices per face
		if (result.indexData.size() != result.numFaces * 3) {
			printf("CleanNonManifoldShells: indexData.size()=%zu != numFaces*3=%u\n",
				result.indexData.size(), result.numFaces * 3);
			ok = false;
		}

		// planeData must hold exactly 1 entry per face
		if (result.planeData.size() != result.numFaces) {
			printf("CleanNonManifoldShells: planeData.size()=%zu != numFaces=%u\n",
				result.planeData.size(), result.numFaces);
			ok = false;
		}

		// vertexData must hold exactly 6 doubles per point
		if (result.vertexData.size() != result.numPoints * VERTEX_FORMAT_SIZE_FLOATS) {
			printf("CleanNonManifoldShells: vertexData.size()=%zu != numPoints*6=%u\n",
				result.vertexData.size(),
				result.numPoints * VERTEX_FORMAT_SIZE_FLOATS);
			ok = false;
		}

		// Every index must reference a valid vertex
		for (size_t i = 0; i < result.indexData.size(); i++) {
			if (result.indexData[i] >= result.numPoints) {
				printf("CleanNonManifoldShells: indexData[%zu]=%u >= numPoints=%u\n",
					i, result.indexData[i], result.numPoints);
				ok = false;
				break;
			}
		}

		// If planes are used, the planes vector must cover all referenced IDs
		if (result.hasPlanes && !result.planes.empty()) {
			for (size_t i = 0; i < result.numFaces; i++) {
				Face f = result.GetFace(i);
				uint32_t planeIdx = f.pId;
				
				if (planeIdx >= result.planes.size()) {
					printf("CleanNonManifoldShells: face %zu has invalid plane ID\n",i);
					ok = false;
					break;
				}
			}
		}

		// Check for NaN in vertex positions
		for (uint32_t i = 0; i < result.numPoints; i++) {
			size_t base = i * VERTEX_FORMAT_SIZE_FLOATS;
			if (std::isnan(result.vertexData[base + 0]) ||
				std::isnan(result.vertexData[base + 1]) ||
				std::isnan(result.vertexData[base + 2])) {
				printf("CleanNonManifoldShells: NaN in vertex %u position\n", i);
				ok = false;
				break;
			}
		}

		if (!ok) {
			printf("CleanNonManifoldShells: consistency check failed\n");
			// no roll-back in debug only! That could lead to a crash happening only in release!
		}
		return ok;
	}

	// ---------------------------------------------------------------
	// Post-boolean cleanup, four phases:
	//
	// Phase A — strip near-degenerate sliver triangles (area < 1e-9 m²)
	//   that accumulate at intersection boundaries across multiple
	//   boolean iterations.
	//
	// Phase B — detect and remove thin membrane regions.
	//   B.1   Double-layer detection: ray-cast from each face centre
	//         along ±normal within THIN_THRESHOLD (1 mm). Opposing
	//         faces with anti-parallel normals -> both removed.
	//   B.2   Manifold-edge component analysis: build connectivity
	//         using only manifold edges (count == 2).  Junction edges
	//         (count 3+) are excluded, naturally separating the solid
	//         body from membrane artifacts.  Signed volume is computed
	//         relative to each component's centroid (flat open
	//         surfaces -> V ≈ 0; closed solids -> V = enclosed vol).
	//         Near-zero-volume components are marked for removal.
	//   B.3   Erosion: faces where ≥ half of edge-neighbours are
	//         already thin-marked are also removed.
	//   B.4   Safety: if > 75 % of faces are thin-marked, disable.
	//   Uses a BVH for B.1 ray queries.
	//
	// Phase C — remove non-manifold open-shell connected components
	//   whose signed volume is near zero (fallback for fragments
	//   fully disconnected from the solid body).
	//
	// Vertex deduplication mirrors SharedPosition::AddPoint
	// (cell size = toleranceVectorEquality, 27-cell neighbourhood,
	// exact distance check).
	// ---------------------------------------------------------------
	inline void CleanNonManifoldShells(fuzzybools::Geometry &result) {
		result.hasPlanes = false;  // rebuild planes, since they could be invalid after fuzzybools::Subtract
		result.planes.clear();
		result.planeData.clear();
		result.buildPlanes();
		Geometry backup = result;
		EdgeMapData emd;
		BuildEdgeMap(result, emd);
		uint32_t openBefore = emd.openEdgeCount;
#if defined( CSG_DEBUG_OUTPUT ) || defined(_DEBUG)
		if (openBefore > 0) {
			DumpGeometry(result, L"CleanNonManifoldShells-entry.obj");
		}
		bool meshValidOnEntry = meshSanityCheck(result);
#endif

		// -- Phase A: strip sliver triangles -----------------------------
		// Threshold 1e-9 m² sits well above the toleranceAddFace filter
		// (~5e-11 m²) so it removes only absolute dregs, while staying
		// far below the smallest real feature in any meter-scale model.
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
					if (areaOfTriangle(a, b, c) >= SLIVER_AREA_THRESHOLD) {
						tmp.AddFace(a, b, c, f.pId);
#ifdef _DEBUG
						if (f.pId == UINT32_MAX) {
							printf("Warning: face %u has no plane assigned\n", i);
						}
#endif
					}
				}
				result = tmp;
				BuildEdgeMap(result, emd); // refresh after Phase A
			}
		}

		// -- Phase A.5: close obvious open-edge loops ----------------------
		// Detect boundary edges (shared by exactly 1 face), build an
		// undirected boundary adjacency graph, trace simple loops
		// (every vertex has degree 2), and fill them with ear-clipping
		// triangulation.  Winding is determined by comparing against
		// the adjacent face's edge direction.  Reuses the EdgeMapData
		// built by BuildEdgeMap to avoid a redundant spatial-hash pass.
		{
			const uint32_t nf = result.numFaces;
			bool loopsClosed = false;

			if (nf >= 3)
			{
				auto &vidPos    = emd.vidPos;
				auto &fids      = emd.fvid;
				auto &fpids     = emd.fpid;
				auto &edgeFaces = emd.edgeFaces;

				// Collect undirected boundary adjacency.
				// An open edge (edgeFaces.size()==1) connects two boundary
				// vertices.  Using undirected edges is more robust than
				// directed ones because face winding around the boundary
				// may be inconsistent after boolean operations.
				std::unordered_map<uint32_t, std::vector<uint32_t>> bdryAdj;
				std::unordered_map<uint32_t, uint32_t> bdryPid;

				for (auto &[ek, fl] : edgeFaces)
				{
					if (fl.size() != 1) continue;
					bdryAdj[ek.v0].push_back(ek.v1);
					bdryAdj[ek.v1].push_back(ek.v0);
					bdryPid[ek.v0] = fpids[fl[0]];
					bdryPid[ek.v1] = fpids[fl[0]];
				}

				// Trace simple loops (every boundary vertex has degree 2)
				std::unordered_set<uint32_t> visited;
				for (auto &[start, adj] : bdryAdj)
				{
					if (visited.count(start)) continue;
					if (adj.size() != 2) continue;

					std::vector<uint32_t> loop;
					loop.push_back(start);
					visited.insert(start);
					uint32_t prev = start;
					uint32_t cur = adj[0];
					bool valid = true;

					while (cur != start)
					{
						if (visited.count(cur)) { valid = false; break; }
						auto it = bdryAdj.find(cur);
						if (it == bdryAdj.end() || it->second.size() != 2)
						{
							valid = false;
							break;
						}
						loop.push_back(cur);
						visited.insert(cur);
						uint32_t next = (it->second[0] == prev)
						              ? it->second[1] : it->second[0];
						prev = cur;
						cur = next;
					}

					if (!valid || loop.size() < 3) continue;

					// Determine correct winding: the fill face for an
					// open edge should use the OPPOSITE edge direction
					// from the existing face.  Check the first edge of
					// the loop against its adjacent face and reverse
					// the loop if it goes the same way.
					{
						EdgeMapData::EKey ek0 = {
							std::min(loop[0], loop[1]),
							std::max(loop[0], loop[1])};
						auto eit = edgeFaces.find(ek0);
						if (eit != edgeFaces.end() && !eit->second.empty())
						{
							uint32_t fi = eit->second[0];
							for (int e = 0; e < 3; e++)
							{
								uint32_t fa = fids[fi][e];
								uint32_t fb = fids[fi][(e + 1) % 3];
								if (fa == loop[0] && fb == loop[1])
								{
									// Loop traverses this edge in the
									// same direction as the face -- the
									// fill face needs the reverse.
									std::reverse(loop.begin(), loop.end());
									break;
								}
								if (fa == loop[1] && fb == loop[0])
									break; // already correct
							}
						}
					}

					// -- Ear-clipping triangulation of the loop ------

					// Average normal from boundary edges
					Vec centroid(0);
					for (uint32_t vid : loop) centroid += vidPos[vid];
					centroid /= static_cast<double>(loop.size());

					Vec avgNormal(0);
					for (size_t i = 0; i < loop.size(); i++)
					{
						const Vec &pa = vidPos[loop[i]];
						const Vec &pb = vidPos[loop[(i + 1) % loop.size()]];
						avgNormal += glm::cross(pa - centroid, pb - centroid);
					}
					double nLen = glm::length(avgNormal);
					if (nLen < 1e-15) continue;
					avgNormal /= nLen;

					// Orthonormal basis for 2D projection
					Vec ref = (std::abs(avgNormal.x) < 0.9)
					        ? Vec(1, 0, 0) : Vec(0, 1, 0);
					Vec uAxis = glm::normalize(glm::cross(avgNormal, ref));
					Vec vAxis = glm::cross(avgNormal, uAxis);

					size_t n = loop.size();
					std::vector<glm::dvec2> pts2d(n);
					for (size_t i = 0; i < n; i++)
						pts2d[i] = glm::dvec2( glm::dot(vidPos[loop[i]], uAxis), glm::dot(vidPos[loop[i]], vAxis));

					// Signed area -> winding direction in 2D
					double sa = 0;
					for (size_t i = 0; i < n; i++)
					{
						size_t j = (i + 1) % n;
						sa += pts2d[i].x * pts2d[j].y - pts2d[j].x * pts2d[i].y;
					}

					// Build index list; reverse for CCW if input is CW
					std::vector<uint32_t> idx(n);
					for (size_t i = 0; i < n; i++)
						idx[i] = static_cast<uint32_t>(i);
					bool reversed = (sa < 0);
					if (reversed)
						std::reverse(idx.begin(), idx.end());

					uint32_t pId = bdryPid.count(loop[0]) ? bdryPid[loop[0]] : 0;

					// Clip ears until only 2 vertices remain
					while (idx.size() > 2)
					{
						bool found = false;
						size_t sz = idx.size();
						for (size_t i = 0; i < sz; i++)
						{
							uint32_t pi = idx[(i + sz - 1) % sz];
							uint32_t ci = idx[i];
							uint32_t ni = idx[(i + 1) % sz];

							glm::dvec2 pa = pts2d[pi];
							glm::dvec2 pb = pts2d[ci];
							glm::dvec2 pc = pts2d[ni];
							double cross = (pb.x - pa.x) * (pc.y - pa.y) -
								(pb.y - pa.y) * (pc.x - pa.x);
							if (cross <= 1e-15) continue;

							// Reject if any other vertex is inside
							bool inside = false;
							for (size_t j = 0; j < sz; j++)
							{
								uint32_t ti = idx[j];
								if (ti == pi || ti == ci || ti == ni)
									continue;
								glm::dvec2 pt = pts2d[ti];
								double d1 = (pt.x - pa.x) * (pb.y - pa.y)
								          - (pt.y - pa.y) * (pb.x - pa.x);
								double d2 = (pt.x - pb.x) * (pc.y - pb.y)
								          - (pt.y - pb.y) * (pc.x - pb.x);
								double d3 = (pt.x - pc.x) * (pa.y - pc.y)
								          - (pt.y - pc.y) * (pa.x - pc.x);
								if (d1 >= 0 && d2 >= 0 && d3 >= 0)
								{
									inside = true;
									break;
								}
							}
							if (inside) continue;

							// Emit triangle with correct 3D winding
							Vec ta = vidPos[loop[pi]];
							Vec tb = vidPos[loop[ci]];
							Vec tc = vidPos[loop[ni]];
							if (reversed)
								result.AddFace(ta, tc, tb, pId);
							else
								result.AddFace(ta, tb, tc, pId);

#ifdef _DEBUG
							if (pId == UINT32_MAX) {
								printf("Warning: no plane assigned\n");
							}
#endif

							idx.erase(idx.begin() + i);
							found = true;
							loopsClosed = true;
							break;
						}
						if (!found) break;
					}
				}
			}

			// Refresh emd if the mesh was modified
			if (loopsClosed) {
				BuildEdgeMap(result, emd);

#if defined( CSG_DEBUG_OUTPUT ) || defined(_DEBUG)
				DumpGeometry(result, L"CleanNonManifoldShells-loopsClosed.obj");
#endif
			}
		}

		// -- Shared setup for Phases B & C --------------------------------
		const uint32_t nFaces = result.numFaces;
		if (nFaces < 4) // need >= 4 faces for any closed solid
		{
			BuildEdgeMap(result, emd);
			if (emd.openEdgeCount > openBefore) result = backup;
			return;
		}

		// -- Step 1: cache face vertices + per-face geometric info --------
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

		// -- Steps 2 & 3: reuse emd (vertex dedup + edge map) -----------
		auto &fvid      = emd.fvid;
		auto &edgeFaces = emd.edgeFaces;

		// -- Step 4: face adjacency (shared by Phase B & C) -------------
		std::vector<std::vector<uint32_t>> faceAdj(nFaces);
		for (auto &[ek, fl] : edgeFaces)
			for (size_t a = 0; a < fl.size(); a++)
				for (size_t b = a + 1; b < fl.size(); b++)
				{
					faceAdj[fl[a]].push_back(fl[b]);
					faceAdj[fl[b]].push_back(fl[a]);
				}

		// ================================================================
		// thin membrane detection
		// ================================================================

		constexpr double THIN_THRESHOLD = 1e-3; // 1 mm
		constexpr double VOLUME_THRESHOLD = 1e-6; // 1 mm³
		std::vector<bool> thinMarked(nFaces, false);

		// Build BVH of the result mesh (used by B.1)
		BVH resultBVH = MakeBVH(result);

		// B.1: double-layer detection — probe along ±normal for nearby
		//      opposing faces (within THIN_THRESHOLD, anti-parallel normals).
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
					[&](uint32_t faceIdx) -> bool
				{
					if (faceIdx == i) return false;
					if (glm::dot(fv[i].normal, fv[faceIdx].normal) > -0.7)
						return false;

					Vec hitPos;
					double t, d_plane;
					if (intersect_ray_triangle(
						fv[i].center, rayEnd,
						fv[faceIdx].a, fv[faceIdx].b, fv[faceIdx].c,
						hitPos, t, d_plane, false))
					{
						double dist = glm::length(hitPos - fv[i].center);
						if (dist > 1e-6 && dist < THIN_THRESHOLD)
						{
							thinMarked[i] = true;
							thinMarked[faceIdx] = true;
							return true;
						}
					}
					return false;
				});
			}
		}

		// B.2: manifold-edge component analysis
		// Build connectivity using ONLY manifold edges (count == 2).
		// Junction edges (count 3+) are excluded, so the solid body
		// and membrane artifacts fall into separate components.
		// Signed volume relative to each component's centroid is used
		// to detect membranes (flat -> V ≈ 0) vs solids (V >> 0).
		{
			std::vector<std::vector<uint32_t>> manifoldAdj(nFaces);
			for (auto &[ek, fl] : edgeFaces)
			{
				if (fl.size() != 2) continue;
				manifoldAdj[fl[0]].push_back(fl[1]);
				manifoldAdj[fl[1]].push_back(fl[0]);
			}

			std::vector<int> mCompId(nFaces, -1);
			int mNumComp = 0;
			for (uint32_t i = 0; i < nFaces; i++)
			{
				if (mCompId[i] >= 0) continue;
				int cid = mNumComp++;
				std::queue<uint32_t> q;
				q.push(i);
				mCompId[i] = cid;
				while (!q.empty())
				{
					uint32_t cur = q.front(); q.pop();
					for (uint32_t nb : manifoldAdj[cur])
						if (mCompId[nb] < 0)
						{
							mCompId[nb] = cid;
							q.push(nb);
						}
				}
			}

			if (mNumComp > 1)
			{
				// Centroid per component (average of face centres)
				struct MCInfo { uint32_t count = 0; Vec centroid{0}; };
				std::vector<MCInfo> mci(mNumComp);
				for (uint32_t i = 0; i < nFaces; i++)
				{
					mci[mCompId[i]].count++;
					mci[mCompId[i]].centroid += fv[i].center;
				}
				for (int c = 0; c < mNumComp; c++)
					if (mci[c].count > 0)
						mci[c].centroid /= static_cast<double>(mci[c].count);

				// Signed volume relative to centroid:
				//   closed surface -> actual enclosed volume (origin-invariant)
				//   flat open surface -> V ≈ 0 (points coplanar with centroid)
				std::vector<double> mVol(mNumComp, 0.0);
				for (uint32_t i = 0; i < nFaces; i++)
				{
					int ci = mCompId[i];
					Vec va = fv[i].a - mci[ci].centroid;
					Vec vb = fv[i].b - mci[ci].centroid;
					Vec vc = fv[i].c - mci[ci].centroid;
					mVol[ci] += glm::dot(va, glm::cross(vb, vc)) / 6.0;
				}

				for (uint32_t i = 0; i < nFaces; i++)
					if (std::abs(mVol[mCompId[i]]) < VOLUME_THRESHOLD)
						thinMarked[i] = true;
			}
		}

		// B.3: iterative erosion — mark narrow bridge faces that
		//       connect the membrane to the solid body.
		//       A face is eroded if:
		//         - its minimum altitude < THIN_THRESHOLD (it is narrow)
		//         - at least half its edge-neighbours are thin-marked
		//       Max 20 iterations.
		for (int iter = 0; iter < 20; iter++)
		{
			bool changed = false;
			for (uint32_t i = 0; i < nFaces; i++)
			{
				if (thinMarked[i]) continue;
				double minH = fv[i].maxEdge > 0
				            ? 2.0 * fv[i].area / fv[i].maxEdge : 0;
				if (minH >= THIN_THRESHOLD) continue;
				int thinNb = 0, totalNb = 0;
				for (uint32_t nb : faceAdj[i])
				{
					totalNb++;
					if (thinMarked[nb]) thinNb++;
				}
				if (totalNb > 0 && thinNb * 2 >= totalNb)
				{
					thinMarked[i] = true;
					changed = true;
				}
			}
			if (!changed) break;
		}

		// B.4: safety — if > 75 % of faces are thin, something went
		//       wrong (e.g. very thin but valid geometry); disable.
		uint32_t thinCount = 0;
		for (uint32_t i = 0; i < nFaces; i++)
			if (thinMarked[i]) thinCount++;
		bool thinEnabled = (thinCount > 0 && thinCount < nFaces * 3 / 4);

		// ================================================================
		// Phase C: non-manifold zero-volume component removal
		// ================================================================
		// BFS connected components on ALL faces (independent of Phase B).
		// Components that are non-manifold AND have near-zero signed
		// volume are discarded as orphaned shell fragments.

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

		// Per-component manifold check + signed volume
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

		// Identify bad components (non-manifold AND near-zero volume)
		std::vector<bool> badComp(numComp, false);
		uint32_t removedByComp = 0;
		if (numComp > 1)
		{
			for (int cid = 0; cid < numComp; cid++)
			{
				bool isManifold = (comps[cid].boundaryEdges == 0);
				bool hasVolume  = (std::abs(comps[cid].volume) >= VOLUME_THRESHOLD);
				if (!isManifold && !hasVolume)
				{
					badComp[cid] = true;
					removedByComp += comps[cid].faceCount;
				}
			}
		}

		// ================================================================
		// Rebuild: keep faces not marked thin AND not in bad components
		// ================================================================
		uint32_t totalRemoved = (thinEnabled ? thinCount : 0) + removedByComp;
		if (totalRemoved == 0 || totalRemoved >= nFaces)
		{
			BuildEdgeMap(result, emd);
			if (emd.openEdgeCount > openBefore) result = backup;
			return;
		}

		Geometry cleaned;
		cleaned.planes = result.planes;
		cleaned.hasPlanes = result.hasPlanes;
		for (uint32_t i = 0; i < nFaces; i++)
		{
			if (thinEnabled && thinMarked[i]) continue;
			if (badComp[compId[i]]) continue;
			cleaned.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);

#ifdef _DEBUG
			if (fv[i].pId == UINT32_MAX) {
				printf("Warning: face %u has no plane assigned\n", i);
			}
#endif
		}
		cleaned.data = result.data; // preserve face-count metadata
		result = cleaned;

		// -- Safety: revert if open-edge count increased ------------------
		BuildEdgeMap(result, emd);
		if (emd.openEdgeCount > openBefore)
			result = backup;

#if defined( CSG_DEBUG_OUTPUT ) || defined(_DEBUG)
		// Sanity-check the result geometry for internal consistency.
		bool meshValidOnExit = meshSanityCheck(result);
		if (!meshValidOnExit) {
			printf("CleanNonManifoldShells: mesh sanity check failed after cleanup!\n");
		}
#endif
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
