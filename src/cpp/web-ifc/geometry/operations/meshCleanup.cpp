#include <queue>
#include <array>
#include <tuple>
#include <functional>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

#include <mapbox/earcut.hpp>

#include <web-ifc/geometry/operations/bim-geometry/booleanUtils.h>
#include <web-ifc/geometry/operations/bim-geometry/geometry.h>
#include <web-ifc/geometry/operations/boolean-utils/clip-mesh.h>
#include <web-ifc/geometry/operations/boolean-utils/shared-position.h>
#include <web-ifc/geometry/operations/meshCleanup.h>
#include <web-ifc/geometry/IfcGeometryProcessor.h>

#if defined(_DEBUG)
#define DUMP_CSG_MESHES
#include "../../test/io_helpers.h"
#include "../../test/dumpToThree.h"
#endif

using namespace fuzzybools;

namespace meshCleanup {
	static int SelectProjectionAxis(const Vec& normal) {
		double ax = std::abs(normal.x);
		double ay = std::abs(normal.y);
		double az = std::abs(normal.z);
		if (az >= ax && az >= ay) return 2;
		if (ay >= ax) return 1;
		return 0;
	}

	static std::array<double, 2> ProjectTo2D(const Vec& p, int dropAxis) {
		switch (dropAxis) {
		case 0: return { p.y, p.z };
		case 1: return { p.x, p.z };
		default: return { p.x, p.y };
		}
	}

	static double SignedArea2D(const std::vector<std::array<double, 2>>& ring) {
		if (ring.size() < 3) return 0.0;

		double area = 0.0;
		for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
			area += ring[j][0] * ring[i][1] - ring[i][0] * ring[j][1];
		}
		return area * 0.5;
	}

	static bool PointInPolygon2D(const std::array<double, 2>& point, const std::vector<std::array<double, 2>>& ring) {
		bool inside = false;
		for (size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
			const auto& a = ring[j];
			const auto& b = ring[i];

			bool crossesScanline = (a[1] > point[1]) != (b[1] > point[1]);
			if (!crossesScanline) continue;

			double xAtY = (b[0] - a[0]) * (point[1] - a[1]) / (b[1] - a[1]) + a[0];
			if (point[0] < xAtY) {
				inside = !inside;
			}
		}
		return inside;
	}

	MeshWatertightInfo isMeshWatertight(const fuzzybools::Geometry& geom) {
		MeshWatertightInfo info;
		info.numFaces = geom.numFaces;
		if (geom.numFaces == 0) {
			return info;
		}

		// Vertices are not shared between triangles, so edges must be compared
		// by position. Snap coordinates to a grid to handle floating-point noise,
		// then count how many triangles share each edge by position. A watertight mesh has every edge shared by exactly 2 triangles.

		constexpr double SNAP = 1e4;  // snap to millimeter precision
		constexpr int STRIDE = fuzzybools::VERTEX_FORMAT_SIZE_FLOATS;

		auto snapCoord = [](double v) -> int64_t {
			return static_cast<int64_t>(std::round(v * SNAP));
			};

		struct Vec3Key {
			int64_t x, y, z;
			bool operator==(const Vec3Key& o) const { return x == o.x && y == o.y && z == o.z; }
		};

		struct Vec3Hash {
			size_t operator()(const Vec3Key& k) const {
				size_t h = std::hash<int64_t>{}(k.x);
				h ^= std::hash<int64_t>{}(k.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
				h ^= std::hash<int64_t>{}(k.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
				return h;
			}
		};

		// Map each unique snapped position to a canonical index
		std::unordered_map<Vec3Key, uint32_t, Vec3Hash> posToIndex;
		std::vector<uint32_t> canonicalIndex(geom.numPoints);

		uint32_t nextIndex = 0;
		for (uint32_t i = 0; i < geom.numPoints; ++i) {
			Vec3Key key{
				snapCoord(geom.vertexData[i * STRIDE + 0]),
				snapCoord(geom.vertexData[i * STRIDE + 1]),
				snapCoord(geom.vertexData[i * STRIDE + 2])
			};
			auto [it, inserted] = posToIndex.emplace(key, nextIndex);
			if (inserted) {
				++nextIndex;
			}
			canonicalIndex[i] = it->second;
		}

		// Count edges using canonical indices
		auto edgeKey = [](uint32_t a, uint32_t b) -> uint64_t {
			if (a > b) std::swap(a, b);
			return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
			};

		std::unordered_map<uint64_t, int> edgeCount;
		edgeCount.reserve(geom.numFaces * 3);

		for (uint32_t f = 0; f < geom.numFaces; ++f) {
			uint32_t i0 = canonicalIndex[geom.indexData[f * 3 + 0]];
			uint32_t i1 = canonicalIndex[geom.indexData[f * 3 + 1]];
			uint32_t i2 = canonicalIndex[geom.indexData[f * 3 + 2]];

			edgeCount[edgeKey(i0, i1)]++;
			edgeCount[edgeKey(i1, i2)]++;
			edgeCount[edgeKey(i2, i0)]++;
		}

		info.numUniqueVertices = nextIndex;
		info.numTotalEdges = static_cast<uint32_t>(edgeCount.size());

		for (const auto& [key, count] : edgeCount) {
			if (count != 2) {
				++info.numOpenEdges;
				if (count == 1) {
					++info.numBoundaryEdges;
				}
			}
		}

		info.watertight = (info.numOpenEdges == 0);
		return info;
	}

	// ---------------------------------------------------------------
	// Patch coplanar holes: find boundary-edge loops and fill them with earcut triangulation when all loop vertices are coplanar.
	// ---------------------------------------------------------------
	static void PatchCoplanarHoles(Geometry& geom) {
		const uint32_t nFaces = geom.numFaces;
		if (nFaces == 0) return;

		// -- vertex deduplication (same spatial-hash as PostBooleanOperationMeshCleanup) --
		const double cellSize = toleranceVectorEquality;
		const double cellSizeSq = cellSize * cellSize;

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
		std::vector<Vec> canonPos;   // canonical id -> position
		uint32_t nextVid = 0;

		auto getCell = [&](const Vec& p) -> GridKey {
			return { static_cast<int64_t>(std::floor(p.x / cellSize)),
					static_cast<int64_t>(std::floor(p.y / cellSize)),
					static_cast<int64_t>(std::floor(p.z / cellSize)) };
		};

		auto findOrAdd = [&](const Vec& p) -> uint32_t {
			auto center = getCell(p);
			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz) {
						GridKey nk = { std::get<0>(center) + dx,
									  std::get<1>(center) + dy,
									  std::get<2>(center) + dz };
						auto it = vtxGrid.find(nk);
						if (it != vtxGrid.end()) {
							for (auto& [id, pos] : it->second) {
								Vec d = p - pos;
								if (glm::dot(d, d) < cellSizeSq)
									return id;
							}
						}
					}
			uint32_t id = nextVid++;
			vtxGrid[center].emplace_back(id, p);
			canonPos.push_back(p);
			return id;
		};

		// Build per-face canonical vertex ids and cache face normals/planeIds
		struct FaceInfo {
			std::array<uint32_t, 3> vid;
			Vec normal;
			uint32_t pId;
		};
		std::vector<FaceInfo> faces(nFaces);
		for (uint32_t i = 0; i < nFaces; i++) {
			Face f = geom.GetFace(i);
			Vec a = geom.GetPoint(f.i0);
			Vec b = geom.GetPoint(f.i1);
			Vec c = geom.GetPoint(f.i2);
			Vec cr = glm::cross(b - a, c - a);
			double len = glm::length(cr);
			faces[i].vid = { findOrAdd(a), findOrAdd(b), findOrAdd(c) };
			faces[i].normal = len > 1e-15 ? cr / len : Vec(0);
			faces[i].pId = static_cast<uint32_t>(f.pId);
		}

		// -- Build undirected edge map and find boundary edges --
		struct EKey {
			uint32_t v0, v1;
			bool operator==(const EKey& o) const { return v0 == o.v0 && v1 == o.v1; }
		};
		struct EKeyHash {
			size_t operator()(const EKey& k) const {
				size_t h = std::hash<uint32_t>()(k.v0);
				h ^= std::hash<uint32_t>()(k.v1) + 0x9e3779b9 + (h << 6) + (h >> 2);
				return h;
			}
		};
		struct EdgeInfo {
			uint32_t count;
			uint32_t faceIdx; // one adjacent face (for normal/planeId lookup)
		};
		std::unordered_map<EKey, EdgeInfo, EKeyHash> edgeMap;
		for (uint32_t i = 0; i < nFaces; i++) {
			auto& v = faces[i].vid;
			for (int e = 0; e < 3; e++) {
				uint32_t va = v[e], vb = v[(e + 1) % 3];
				EKey ek = { std::min(va, vb), std::max(va, vb) };
				auto it = edgeMap.find(ek);
				if (it != edgeMap.end())
					it->second.count++;
				else
					edgeMap[ek] = { 1, i };
			}
		}

		// Build boundary vertex adjacency (undirected)
		std::unordered_map<uint32_t, std::vector<uint32_t>> boundaryAdj;
		std::unordered_map<uint32_t, uint32_t> vertexAdjFace; // vertex -> one adjacent boundary face
		for (auto& [ek, info] : edgeMap) {
			if (info.count != 1) continue;
			boundaryAdj[ek.v0].push_back(ek.v1);
			boundaryAdj[ek.v1].push_back(ek.v0);
			vertexAdjFace[ek.v0] = info.faceIdx;
			vertexAdjFace[ek.v1] = info.faceIdx;
		}

		if (boundaryAdj.empty()) return; // mesh is already closed

		// -- Trace boundary loops via adjacency --
		std::unordered_set<uint32_t> visited;
		std::vector<std::vector<uint32_t>> loops;
		std::vector<uint32_t> loopAdjacentFace; // one adjacent face per loop

		for (auto& [start, _] : boundaryAdj) {
			if (visited.count(start)) continue;
			// Skip vertices with valence != 2 (non-manifold junctions)
			if (boundaryAdj[start].size() != 2) continue;

			std::vector<uint32_t> loop;
			uint32_t prev = UINT32_MAX;
			uint32_t cur = start;
			bool valid = true;
			while (true) {
				visited.insert(cur);
				loop.push_back(cur);

				auto& neighbors = boundaryAdj[cur];
				if (neighbors.size() != 2) {
					valid = false;
					break;
				}

				uint32_t next = (neighbors[0] == prev) ? neighbors[1] : neighbors[0];
				prev = cur;
				cur = next;

				if (cur == start) break;

				if (visited.count(cur)) {
					valid = false;
					break;
				}

				if (loop.size() > canonPos.size()) {
					valid = false;
					break;
				}
			}

			if (valid && loop.size() >= 3) {
				loops.push_back(std::move(loop));
				loopAdjacentFace.push_back(vertexAdjFace[start]);
			}
		}

		if (loops.empty()) return;

		const double PLANE_EPS = 1e-5;
		struct LoopInfo {
			std::vector<uint32_t> vid;
			Vec planeNormal;
			Vec planePoint;
			Vec adjNormal;
			uint32_t adjPlaneId;
		};

		std::unordered_map<uint32_t, Vec> rawBoundaryPos;
		rawBoundaryPos.reserve(boundaryAdj.size());

		auto getRawBoundaryPos = [&](uint32_t vid) -> Vec {
			auto rawIt = rawBoundaryPos.find(vid);
			if (rawIt != rawBoundaryPos.end()) {
				return rawIt->second;
			}

			Vec raw = canonPos[vid];
			auto it = vertexAdjFace.find(vid);
			if (it != vertexAdjFace.end()) {
				uint32_t fIdx = it->second;
				Face f = geom.GetFace(fIdx);
				size_t rawIdx[3] = { (size_t)f.i0, (size_t)f.i1, (size_t)f.i2 };
				for (int j = 0; j < 3; j++) {
					if (faces[fIdx].vid[j] == vid) {
						raw = geom.GetPoint(rawIdx[j]);
						break;
					}
				}
			}

			rawBoundaryPos[vid] = raw;
			return raw;
		};

		std::vector<LoopInfo> validLoops;
		validLoops.reserve(loops.size());
		for (size_t li = 0; li < loops.size(); li++) {
			auto& loop = loops[li];
			if (loop.size() < 3) continue;

			// Find plane normal from first non-collinear triple
			Vec planeNormal(0);
			Vec planePoint = canonPos[loop[0]];
			bool foundPlane = false;
			for (size_t i = 1; i + 1 < loop.size(); i++) {
				Vec e1 = canonPos[loop[i]] - planePoint;
				Vec e2 = canonPos[loop[i + 1]] - planePoint;
				Vec cr = glm::cross(e1, e2);
				double len = glm::length(cr);
				if (len > 1e-12) {
					planeNormal = cr / len;
					foundPlane = true;
					break;
				}
			}
			if (!foundPlane) continue;

			// Check all vertices lie on the plane
			bool coplanar = true;
			for (auto vid : loop) {
				double dist = std::abs(glm::dot(canonPos[vid] - planePoint, planeNormal));
				if (dist > PLANE_EPS) {
					coplanar = false;
					break;
				}
			}
			if (!coplanar) continue;

			validLoops.push_back({ loop, planeNormal, planePoint, faces[loopAdjacentFace[li]].normal,
				faces[loopAdjacentFace[li]].pId });
		}

		if (validLoops.empty()) return;

		auto samePlane = [&](const LoopInfo& a, const LoopInfo& b) -> bool {
			if (std::abs(glm::dot(a.planeNormal, b.planeNormal)) < 1.0 - 1e-6) {
				return false;
			}
			double dist = std::abs(glm::dot(b.planePoint - a.planePoint, a.planeNormal));
			return dist <= PLANE_EPS;
		};

		std::vector<bool> grouped(validLoops.size(), false);
		for (size_t i = 0; i < validLoops.size(); i++) {
			if (grouped[i]) continue;

			std::vector<size_t> planeGroup = { i };
			grouped[i] = true;
			for (size_t j = i + 1; j < validLoops.size(); j++) {
				if (grouped[j]) continue;
				if (!samePlane(validLoops[i], validLoops[j])) continue;
				grouped[j] = true;
				planeGroup.push_back(j);
			}

			struct ProjectedLoop {
				size_t loopIndex;
				std::vector<std::array<double, 2>> projected;
				double area = 0.0;
				int parent = -1;
				int depth = 0;
			};

			const int dropAxis = SelectProjectionAxis(validLoops[planeGroup[0]].planeNormal);
			std::vector<ProjectedLoop> projectedLoops;
			projectedLoops.reserve(planeGroup.size());

			for (size_t loopIndex : planeGroup) {
				ProjectedLoop projectedLoop;
				projectedLoop.loopIndex = loopIndex;
				projectedLoop.projected.reserve(validLoops[loopIndex].vid.size());
				for (uint32_t vid : validLoops[loopIndex].vid) {
					projectedLoop.projected.push_back(ProjectTo2D(canonPos[vid], dropAxis));
				}
				projectedLoop.area = std::abs(SignedArea2D(projectedLoop.projected));
				if (projectedLoop.area > 1e-12) {
					projectedLoops.push_back(std::move(projectedLoop));
				}
			}

			for (size_t loopIdx = 0; loopIdx < projectedLoops.size(); loopIdx++) {
				const auto samplePoint = projectedLoops[loopIdx].projected[0];
				double bestParentArea = std::numeric_limits<double>::max();

				for (size_t candidateIdx = 0; candidateIdx < projectedLoops.size(); candidateIdx++) {
					if (candidateIdx == loopIdx) continue;
					if (projectedLoops[candidateIdx].area <= projectedLoops[loopIdx].area + 1e-12) continue;
					if (!PointInPolygon2D(samplePoint, projectedLoops[candidateIdx].projected)) continue;

					if (projectedLoops[candidateIdx].area < bestParentArea) {
						bestParentArea = projectedLoops[candidateIdx].area;
						projectedLoops[loopIdx].parent = static_cast<int>(candidateIdx);
					}
				}
			}

			for (auto& projectedLoop : projectedLoops) {
				for (int parent = projectedLoop.parent; parent >= 0; parent = projectedLoops[parent].parent) {
					projectedLoop.depth++;
				}
			}

			for (size_t outerIdx = 0; outerIdx < projectedLoops.size(); outerIdx++) {
				if ((projectedLoops[outerIdx].depth & 1) != 0) continue;

				const LoopInfo& outerLoop = validLoops[projectedLoops[outerIdx].loopIndex];

				std::vector<std::vector<std::array<double, 2>>> polygon;
				std::vector<uint32_t> polygonVertexIds;
				polygon.reserve(projectedLoops.size());
				polygonVertexIds.reserve(outerLoop.vid.size());

				polygon.push_back(projectedLoops[outerIdx].projected);
				polygonVertexIds.insert(polygonVertexIds.end(), outerLoop.vid.begin(), outerLoop.vid.end());

				for (size_t holeIdx = 0; holeIdx < projectedLoops.size(); holeIdx++) {
					if (projectedLoops[holeIdx].parent != static_cast<int>(outerIdx)) continue;
					polygon.push_back(projectedLoops[holeIdx].projected);

					const auto& holeLoop = validLoops[projectedLoops[holeIdx].loopIndex];
					polygonVertexIds.insert(polygonVertexIds.end(), holeLoop.vid.begin(), holeLoop.vid.end());
				}

				std::vector<uint32_t> indices = mapbox::earcut<uint32_t>(polygon);
				if (indices.size() < 3) continue;

				Vec fillNormal(0);
				for (size_t tri = 0; tri < indices.size(); tri += 3) {
					Vec ta = getRawBoundaryPos(polygonVertexIds[indices[tri]]);
					Vec tb = getRawBoundaryPos(polygonVertexIds[indices[tri + 1]]);
					Vec tc = getRawBoundaryPos(polygonVertexIds[indices[tri + 2]]);
					fillNormal += glm::cross(tb - ta, tc - ta);
				}

				bool flipWinding = glm::dot(fillNormal, outerLoop.adjNormal) < 0;
				Vec n = flipWinding ? -outerLoop.planeNormal : outerLoop.planeNormal;

				for (size_t tri = 0; tri < indices.size(); tri += 3) {
					Vec va = getRawBoundaryPos(polygonVertexIds[indices[tri]]);
					Vec vb = getRawBoundaryPos(polygonVertexIds[indices[tri + 1]]);
					Vec vc = getRawBoundaryPos(polygonVertexIds[indices[tri + 2]]);

					// Use index-based AddFace to avoid rejection of near-degenerate triangles by computeSafeNormal.
					// These triangles may have near-zero area but are topologically necessary to close the hole.
					if (flipWinding) {
						geom.AddPoint(va, n);
						geom.AddPoint(vc, n);
						geom.AddPoint(vb, n);
					}
					else {
						geom.AddPoint(va, n);
						geom.AddPoint(vb, n);
						geom.AddPoint(vc, n);
					}
					geom.AddFace(geom.numPoints - 3, geom.numPoints - 2, geom.numPoints - 1,
						outerLoop.adjPlaneId);
				}
			}
		}
		if (geom.numFaces > nFaces) {
			// faces have been added -> rebuild planes
			geom.hasPlanes = false;

#ifdef DUMP_CSG_MESHES
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(geom);
			webifc::io::DumpIfcGeometry(inputWebIfc, "PostBooleanOperationMeshCleanup-patched.obj");
#endif
		}
	}

	// ---------------------------------------------------------------
	// Post-boolean cleanup, five phases: 
	// Phase A - strip near-degenerate sliver triangles (area < 1e-9 m^2)
	//   that accumulate at intersection boundaries across multiple
	// boolean iterations. 
	// Phase B - detect and remove thin membrane regions.
	//   B.1   Double-layer detection: ray-cast from each face centre along ±normal within THIN_THRESHOLD (1 mm). 
	//		   Opposing faces with anti-parallel normals -> both removed.
	//   B.2   Manifold-edge component analysis: build connectivity using only manifold edges (count == 2). 
	//         Junction edges (count 3+) are excluded, naturally separating the solid
	//         body from membrane artifacts.  Signed volume is computed relative to each component's 
	//         centroid (flat open surfaces -> V ≈ 0; closed solids -> V = enclosed vol).
	//         Near-zero-volume components are marked for removal.
	//   B.3   Erosion: faces where ≥ half of edge-neighbours are
	//         already thin-marked are also removed.
	//   B.4   Safety: if > 75 % of faces are thin-marked, disable.
	// Uses a BVH for B.1 ray queries. 
	// Phase C - remove non-manifold open-shell connected components whose signed volume is near zero
	// (fallback for fragments fully disconnected from the solid body). 
	// Phase D - patch coplanar holes: trace boundary-edge loops and fill in-plane holes with earcut triangulation. 
	// Vertex deduplication mirrors SharedPosition::AddPoint (cell size = toleranceVectorEquality, 27-cell 
	// neighbourhood, exact distance check).
	// ---------------------------------------------------------------
	void PostBooleanOperationMeshCleanup(fuzzybools::Geometry& input) {
		if (input.numFaces > 8000) {
#ifdef _DEBUG
			std::cout << "PostBooleanOperationMeshCleanup: skipping mesh with " << input.numFaces << " faces" << std::endl;
#endif
			return;
		}
		fuzzybools::Geometry workingMesh = input;
		auto meshInfoOnEntry = meshCleanup::isMeshWatertight(input);

#ifdef DUMP_CSG_MESHES
		if (!meshInfoOnEntry.watertight) {
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(input);
			webifc::io::DumpIfcGeometry(inputWebIfc, "PostBooleanOperationMeshCleanup-input.obj");
		}
#endif

		// -- Phase A: strip sliver triangles -----------------------------
		// Threshold 1e-9 m^2 sits well above the toleranceAddFace filter (~5e-11 m^2) so it removes only absolute dregs, while staying
		// far below the smallest real feature in any meter-scale model.
		{
			constexpr double SLIVER_AREA_THRESHOLD = 1e-9;
			const uint32_t n = workingMesh.numFaces;
			uint32_t sliverCount = 0;
			for (uint32_t i = 0; i < n; i++) {
				Face f = workingMesh.GetFace(i);
				double area = areaOfTriangle(workingMesh.GetPoint(f.i0), workingMesh.GetPoint(f.i1), workingMesh.GetPoint(f.i2));
				if (area < SLIVER_AREA_THRESHOLD)
					sliverCount++;
			}

			if (sliverCount > 0) {
				Geometry tmp;
				tmp.planes = workingMesh.planes;
				tmp.hasPlanes = workingMesh.hasPlanes;
				tmp.data = workingMesh.data;
				for (uint32_t i = 0; i < n; i++) {
					Face f = workingMesh.GetFace(i);
					Vec a = workingMesh.GetPoint(f.i0);
					Vec b = workingMesh.GetPoint(f.i1);
					Vec c = workingMesh.GetPoint(f.i2);
					if (areaOfTriangle(a, b, c) >= SLIVER_AREA_THRESHOLD)
						tmp.AddFace(a, b, c, f.pId);
				}
				workingMesh = tmp;
			}
		}

		// -- Shared setup for Phases B & C 
		const uint32_t nFaces = workingMesh.numFaces;

		// -- Step 1: cache face vertices + per-face geometric info
		struct FV {
			Vec a, b, c;
			uint32_t pId;
			Vec center, normal;
			double area, maxEdge;
		};
		std::vector<FV> fv(nFaces);
		for (uint32_t i = 0; i < nFaces; i++) {
			Face f = workingMesh.GetFace(i);
			Vec a = workingMesh.GetPoint(f.i0);
			Vec b = workingMesh.GetPoint(f.i1);
			Vec c = workingMesh.GetPoint(f.i2);
			Vec crossP = glm::cross(b - a, c - a);
			double crossLen = glm::length(crossP);
			double e0 = glm::length(b - a);
			double e1 = glm::length(c - b);
			double e2 = glm::length(a - c);
			fv[i] = { a, b, c,
					 static_cast<uint32_t>(f.pId),
					 (a + b + c) / 3.0,
					 crossLen > 1e-15 ? crossP / crossLen : Vec(0),
					 crossLen * 0.5,
					 std::max(e0, std::max(e1, e2)) };
		}

		// -- Step 2: spatial-hash vertex deduplication
		// Mirror SharedPosition::AddPoint: grid cell = floor(coord/cellSize), 27-cell neighbourhood search, exact Euclidean distance check.
		const double cellSize = toleranceVectorEquality * 1.0;          // 1e-4
		const double cellSizeSq = cellSize * cellSize;

		using GridKey = std::tuple<int64_t, int64_t, int64_t>;
		struct GridKeyHash {
			size_t operator()(const GridKey& k) const {
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

		auto getCell = [&](const Vec& p) -> GridKey {
			return { static_cast<int64_t>(std::floor(p.x / cellSize)),
					static_cast<int64_t>(std::floor(p.y / cellSize)),
					static_cast<int64_t>(std::floor(p.z / cellSize)) };
			};

		auto findOrAdd = [&](const Vec& p) -> uint32_t {
			auto center = getCell(p);
			// search 27 neighbours
			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz) {
						GridKey nk = { std::get<0>(center) + dx,
									  std::get<1>(center) + dy,
									  std::get<2>(center) + dz };
						auto it = vtxGrid.find(nk);
						if (it != vtxGrid.end()) {
							for (auto& [id, pos] : it->second) {
								Vec d = p - pos;
								if (glm::dot(d, d) < cellSizeSq)
									return id;
							}
						}
					}
			// not found - insert
			uint32_t id = nextVid++;
			vtxGrid[center].emplace_back(id, p);
			return id;
			};

		for (uint32_t i = 0; i < nFaces; i++) {
			fvid[i][0] = findOrAdd(fv[i].a);
			fvid[i][1] = findOrAdd(fv[i].b);
			fvid[i][2] = findOrAdd(fv[i].c);
		}

		// -- Step 3: edge -> face adjacency ------------------------------
		struct EKey {
			uint32_t v0, v1;
			bool operator==(const EKey& o) const { return v0 == o.v0 && v1 == o.v1; }
		};
		struct EKeyHash {
			size_t operator()(const EKey& k) const {
				size_t h = std::hash<uint32_t>()(k.v0);
				h ^= std::hash<uint32_t>()(k.v1) + 0x9e3779b9 + (h << 6) + (h >> 2);
				return h;
			}
		};
		std::unordered_map<EKey, std::vector<uint32_t>, EKeyHash> edgeFaces;
		for (uint32_t i = 0; i < nFaces; i++)
			for (int e = 0; e < 3; e++) {
				uint32_t va = fvid[i][e];
				uint32_t vb = fvid[i][(e + 1) % 3];
				EKey ek = { std::min(va, vb), std::max(va, vb) };
				edgeFaces[ek].push_back(i);
			}

		// -- Step 4: face adjacency (shared by Phase B & C) -------------
		std::vector<std::vector<uint32_t>> faceAdj(nFaces);
		for (auto& [ek, fl] : edgeFaces)
			for (size_t a = 0; a < fl.size(); a++)
				for (size_t b = a + 1; b < fl.size(); b++) {
					faceAdj[fl[a]].push_back(fl[b]);
					faceAdj[fl[b]].push_back(fl[a]);
				}

		// Phase B: thin membrane detection
		constexpr double THIN_THRESHOLD = 1e-3; // 1 mm
		constexpr double VOLUME_THRESHOLD = 1e-6; // 1 mm³
		std::vector<bool> thinMarked(nFaces, false);

		// Build BVH (Bounding Volume Hierarchy) of the result mesh (used by B.1)
		BVH resultBVH = MakeBVH(workingMesh);

		// B.1: double-layer detection - probe along ±normal for nearby
		//      opposing faces (within THIN_THRESHOLD, anti-parallel normals).
		for (uint32_t i = 0; i < nFaces; i++) {
			if (thinMarked[i]) continue;
			if (glm::length(fv[i].normal) < 0.5) continue;

			for (int sign = -1; sign <= 1; sign += 2) {
				if (thinMarked[i]) break;

				Vec rayDir = fv[i].normal * static_cast<double>(sign);
				Vec rayEnd = fv[i].center + rayDir * THIN_THRESHOLD;

				resultBVH.IntersectRay(fv[i].center, rayDir, [&](uint32_t faceIdx) -> bool {
						if (faceIdx == i) return false;
						if (glm::dot(fv[i].normal, fv[faceIdx].normal) > -0.7)
							return false;

						Vec hitPos;
						double t, d_plane;
						if (intersect_ray_triangle( fv[i].center, rayEnd, fv[faceIdx].a, fv[faceIdx].b,
							fv[faceIdx].c, hitPos, t, d_plane, false)) {
							double dist = glm::length(hitPos - fv[i].center);
							if (dist > 1e-6 && dist < THIN_THRESHOLD) {
								thinMarked[i] = true;
								thinMarked[faceIdx] = true;
								return true;
							}
						}
						return false;
					});
			}
		}

		// B.2: manifold-edge component analysis Build connectivity using ONLY manifold edges (count == 2).
		// Junction edges (count 3+) are excluded, so the solid body and membrane artifacts fall into separate components.
		// Signed volume relative to each component's centroid is used to detect membranes (flat -> V ≈ 0) vs solids (V >> 0).
		{
			std::vector<std::vector<uint32_t>> manifoldAdj(nFaces);
			for (auto& [ek, fl] : edgeFaces) {
				if (fl.size() != 2) continue;
				manifoldAdj[fl[0]].push_back(fl[1]);
				manifoldAdj[fl[1]].push_back(fl[0]);
			}

			std::vector<int> mCompId(nFaces, -1);
			int mNumComp = 0;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (mCompId[i] >= 0) continue;
				int cid = mNumComp++;
				std::queue<uint32_t> q;
				q.push(i);
				mCompId[i] = cid;
				while (!q.empty()) {
					uint32_t cur = q.front(); q.pop();
					for (uint32_t nb : manifoldAdj[cur])
						if (mCompId[nb] < 0) {
							mCompId[nb] = cid;
							q.push(nb);
						}
				}
			}

			if (mNumComp > 1) {
				// Centroid per component (average of face centres)
				struct MCInfo { uint32_t count = 0; Vec centroid{ 0 }; };
				std::vector<MCInfo> mci(mNumComp);
				for (uint32_t i = 0; i < nFaces; i++) {
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
				for (uint32_t i = 0; i < nFaces; i++) {
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

		// B.3: iterative erosion - mark narrow bridge faces that connect the membrane to the solid body.
		//       A face is eroded if:
		//         - its minimum altitude < THIN_THRESHOLD (it is narrow)
		//         - at least half its edge-neighbours are thin-marked
		//       Max 20 iterations.
		for (int iter = 0; iter < 20; iter++) {
			bool changed = false;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (thinMarked[i]) continue;
				double minH = fv[i].maxEdge > 0 ? 2.0 * fv[i].area / fv[i].maxEdge : 0;
				if (minH >= THIN_THRESHOLD) continue;
				int thinNb = 0, totalNb = 0;
				for (uint32_t nb : faceAdj[i]) {
					totalNb++;
					if (thinMarked[nb]) thinNb++;
				}
				if (totalNb > 0 && thinNb * 2 >= totalNb) {
					thinMarked[i] = true;
					changed = true;
				}
			}
			if (!changed) break;
		}

		// B.4: safety - if > 75 % of faces are thin, something went wrong (e.g. very thin but valid geometry); disable.
		uint32_t thinCount = 0;
		for (uint32_t i = 0; i < nFaces; i++)
			if (thinMarked[i]) thinCount++;
		bool thinEnabled = (thinCount > 0 && thinCount < nFaces * 3 / 4);

		// Phase C: non-manifold zero-volume component removal
		// BFS connected components on ALL faces (independent of Phase B). Components that are non-manifold AND have 
		// near-zero signed volume are discarded as orphaned shell fragments.
		std::vector<int> compId(nFaces, -1);
		int numComp = 0;
		for (uint32_t i = 0; i < nFaces; i++) {
			if (compId[i] >= 0) continue;
			int cid = numComp++;
			std::queue<uint32_t> q;
			q.push(i);
			compId[i] = cid;
			while (!q.empty()) {
				uint32_t cur = q.front(); q.pop();
				for (uint32_t nb : faceAdj[cur])
					if (compId[nb] < 0) {
						compId[nb] = cid;
						q.push(nb);
					}
			}
		}

		// Per-component manifold check + signed volume
		struct CompInfo {
			uint32_t faceCount = 0;
			uint32_t boundaryEdges = 0; // edges shared by != 2 faces
			Vec      centroid{ 0 };
			double   volume = 0;
		};
		std::vector<CompInfo> comps(numComp);

		for (uint32_t i = 0; i < nFaces; i++) {
			int cid = compId[i];
			comps[cid].faceCount++;
			comps[cid].centroid += fv[i].center;
		}
		for (auto& comp : comps)
			if (comp.faceCount > 0)
				comp.centroid /= static_cast<double>(comp.faceCount);
		for (uint32_t i = 0; i < nFaces; i++) {
			int cid = compId[i];
			Vec va = fv[i].a - comps[cid].centroid;
			Vec vb = fv[i].b - comps[cid].centroid;
			Vec vc = fv[i].c - comps[cid].centroid;
			comps[cid].volume += glm::dot(va, glm::cross(vb, vc)) / 6.0;
		}
		for (auto& [ek, fl] : edgeFaces)
			if (fl.size() != 2)
				comps[compId[fl[0]]].boundaryEdges++;

		// Identify bad components (non-manifold AND near-zero volume)
		std::vector<bool> badComp(numComp, false);
		uint32_t removedByComp = 0;
		if (numComp > 1) {
			for (int cid = 0; cid < numComp; cid++) {
				bool isManifold = (comps[cid].boundaryEdges == 0);
				bool hasVolume = (std::abs(comps[cid].volume) >= VOLUME_THRESHOLD);
				if (!isManifold && !hasVolume) {
					badComp[cid] = true;
					removedByComp += comps[cid].faceCount;
				}
			}
		}

		// Rebuild: keep faces not marked thin AND not in bad components
		uint32_t totalRemoved = (thinEnabled ? thinCount : 0) + removedByComp;
		if (totalRemoved == 0 || totalRemoved >= nFaces) {
			if (meshInfoOnEntry.watertight) {
				return;
			}
#ifdef _DEBUG
			//std::cout << "Not able to fix open mesh by removing degenerated faces. Continue with hole patching" << std::endl;
#endif
		}

		Geometry cleaned;
		cleaned.planes = workingMesh.planes;
		cleaned.hasPlanes = workingMesh.hasPlanes;
		for (uint32_t i = 0; i < nFaces; i++) {
			if (thinEnabled && thinMarked[i]) continue;
			if (badComp[compId[i]]) continue;
			cleaned.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
		}
		cleaned.data = workingMesh.data; // preserve face-count metadata

		// -- Phase D: patch coplanar holes --
		PatchCoplanarHoles(cleaned);

		auto meshInfoOnExit = meshCleanup::isMeshWatertight(cleaned);
		bool isWatertightOnExit = meshInfoOnExit.watertight;
#ifdef DUMP_CSG_MESHES
		if (!isWatertightOnExit) {
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(cleaned);
			webifc::io::DumpIfcGeometry(inputWebIfc, "PostBooleanOperationMeshCleanup-cleaned.obj");
		}
#endif

		//if (isWatertightOnExit) {
		if(meshInfoOnEntry.numOpenEdges > meshInfoOnExit.numOpenEdges){
			input = std::move(cleaned);
		}
	}
}
