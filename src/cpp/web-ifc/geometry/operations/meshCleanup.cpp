#include <queue>
#include <array>
#include <tuple>
#include <filesystem>
#include <functional>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>

#include <mapbox/earcut.hpp>
#include <CDT.h>

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
		// by position. Use spatial-hash vertex deduplication with 27-cell
		// neighbourhood search and exact distance check (same as the rest of
		// this file) to avoid grid-boundary mismatches that a simple snap can
		// cause. Then count how many triangles share each edge by position.

		constexpr int STRIDE = fuzzybools::VERTEX_FORMAT_SIZE_FLOATS;
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
		std::vector<uint32_t> canonicalIndex(geom.numPoints);
		uint32_t nextIndex = 0;

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
			uint32_t id = nextIndex++;
			vtxGrid[center].emplace_back(id, p);
			return id;
		};

		for (uint32_t i = 0; i < geom.numPoints; ++i) {
			Vec p(geom.vertexData[i * STRIDE + 0],
				  geom.vertexData[i * STRIDE + 1],
				  geom.vertexData[i * STRIDE + 2]);
			canonicalIndex[i] = findOrAdd(p);
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
			if (count == 1) {
				++info.numOpenEdges;
			}
		}

		info.watertight = (info.numOpenEdges == 0);
		return info;
	}

	// ---------------------------------------------------------------
	// Remove disconnected non-manifold zero-volume components.
	// Safe to run at any point -- purely topological, no heuristics.
	// ---------------------------------------------------------------
	static void RemoveDisconnectedFragments(Geometry& workingMesh, std::string step,
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
		const uint32_t nFaces = workingMesh.numFaces;
		if (nFaces == 0) { meshInfoResult = meshInfoInput; return; }

		constexpr double VOLUME_THRESHOLD = 1e-6;

		// Cache face vertices
		struct FV { Vec a, b, c; uint32_t pId; Vec center; };
		std::vector<FV> fv(nFaces);
		for (uint32_t i = 0; i < nFaces; i++) {
			Face f = workingMesh.GetFace(i);
			Vec a = workingMesh.GetPoint(f.i0);
			Vec b = workingMesh.GetPoint(f.i1);
			Vec c = workingMesh.GetPoint(f.i2);
			fv[i] = { a, b, c, static_cast<uint32_t>(f.pId), (a + b + c) / 3.0 };
		}

		// Vertex dedup
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
		std::vector<std::array<uint32_t, 3>> fvid(nFaces);
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
								if (glm::dot(d, d) < cellSizeSq) return id;
							}
						}
					}
			uint32_t id = nextVid++;
			vtxGrid[center].emplace_back(id, p);
			return id;
		};
		for (uint32_t i = 0; i < nFaces; i++) {
			fvid[i][0] = findOrAdd(fv[i].a);
			fvid[i][1] = findOrAdd(fv[i].b);
			fvid[i][2] = findOrAdd(fv[i].c);
		}

		// Edge -> face adjacency
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
				uint32_t va = fvid[i][e], vb = fvid[i][(e + 1) % 3];
				edgeFaces[{std::min(va, vb), std::max(va, vb)}].push_back(i);
			}

		// Face adjacency + BFS components
		std::vector<std::vector<uint32_t>> faceAdj(nFaces);
		for (auto& [ek, fl] : edgeFaces)
			for (size_t a = 0; a < fl.size(); a++)
				for (size_t b = a + 1; b < fl.size(); b++) {
					faceAdj[fl[a]].push_back(fl[b]);
					faceAdj[fl[b]].push_back(fl[a]);
				}

		std::vector<int> compId(nFaces, -1);
		int numComp = 0;
		for (uint32_t i = 0; i < nFaces; i++) {
			if (compId[i] >= 0) continue;
			int cid = numComp++;
			std::queue<uint32_t> q;
			q.push(i); compId[i] = cid;
			while (!q.empty()) {
				uint32_t cur = q.front(); q.pop();
				for (uint32_t nb : faceAdj[cur])
					if (compId[nb] < 0) { compId[nb] = cid; q.push(nb); }
			}
		}

		if (numComp <= 1) { meshInfoResult = meshInfoInput; return; }

		// Per-component analysis
		struct CompInfo { uint32_t faceCount = 0; uint32_t boundaryEdges = 0; Vec centroid{0}; double volume = 0; };
		std::vector<CompInfo> comps(numComp);
		for (uint32_t i = 0; i < nFaces; i++) {
			comps[compId[i]].faceCount++;
			comps[compId[i]].centroid += fv[i].center;
		}
		for (auto& comp : comps)
			if (comp.faceCount > 0) comp.centroid /= static_cast<double>(comp.faceCount);
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

		std::vector<bool> badComp(numComp, false);
		uint32_t removedByComp = 0;
		for (int cid = 0; cid < numComp; cid++) {
			bool isManifold = (comps[cid].boundaryEdges == 0);
			bool hasVolume = (std::abs(comps[cid].volume) >= VOLUME_THRESHOLD);
			if (!isManifold && !hasVolume) {
				badComp[cid] = true;
				removedByComp += comps[cid].faceCount;
			}
		}

		if (removedByComp == 0 || removedByComp >= nFaces) {
			meshInfoResult = meshInfoInput;
			return;
		}

		Geometry cleaned;
		cleaned.planes = workingMesh.planes;
		cleaned.hasPlanes = workingMesh.hasPlanes;
		for (uint32_t i = 0; i < nFaces; i++) {
			if (badComp[compId[i]]) continue;
			cleaned.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
		}
		cleaned.data = workingMesh.data;

		auto infoCleaned = meshCleanup::isMeshWatertight(cleaned);
		if (infoCleaned.numOpenEdges < meshInfoInput.numOpenEdges) {
			workingMesh = std::move(cleaned);
			meshInfoResult = infoCleaned;
		}
		else {
			meshInfoResult = meshInfoInput;
		}
	}

	// -- Phase A: strip sliver triangles -----------------------------
		// Threshold 1e-9 m^2 sits well above the toleranceAddFace filter (~5e-11 m^2) so it removes only absolute dregs, while staying
		// far below the smallest real feature in any meter-scale model.
		// Additionally, detect altitude-based slivers: triangles with very long edges but nearly collinear vertices
		// (minAltitude < toleranceVectorEquality). For these, snap the tip vertex onto the opposite edge
		// so that adjacent faces can later be split at the snap point (T-junction resolution).
	void removeDegeneratedTriangles(Geometry& workingMesh, std::string step,
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
		constexpr double SLIVER_AREA_THRESHOLD = 1e-9;
		const double SLIVER_ALTITUDE_THRESHOLD = toleranceVectorEquality; // 1e-4 m
		const uint32_t n = workingMesh.numFaces;

		// Detect slivers and build a snap map: original tip position -> projected position on opposite edge.
		// Uses spatial hash to match vertices across faces.
		const double snapCellSize = toleranceVectorEquality;
		const double snapCellSizeSq = snapCellSize * snapCellSize;

		using SnapGridKey = std::tuple<int64_t, int64_t, int64_t>;
		struct SnapGridKeyHash {
			size_t operator()(const SnapGridKey& k) const {
				size_t h = std::hash<int64_t>()(std::get<0>(k));
				h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
				h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
				return h;
			}
		};

		// Snap map entries: source position -> projected position
		struct SnapEntry {
			Vec source;
			Vec target;
		};
		std::vector<SnapEntry> snapEntries;

		auto getSnapCell = [&](const Vec& p) -> SnapGridKey {
			return { static_cast<int64_t>(std::floor(p.x / snapCellSize)),
					static_cast<int64_t>(std::floor(p.y / snapCellSize)),
					static_cast<int64_t>(std::floor(p.z / snapCellSize)) };
			};

		// Build snap entries from altitude-based slivers
		for (uint32_t i = 0; i < n; i++) {
			Face f = workingMesh.GetFace(i);
			Vec va = workingMesh.GetPoint(f.i0);
			Vec vb = workingMesh.GetPoint(f.i1);
			Vec vc = workingMesh.GetPoint(f.i2);

			double area = areaOfTriangle(va, vb, vc);
			if (area < SLIVER_AREA_THRESHOLD) continue; // already caught by area threshold

			// Compute edge lengths and find longest edge
			double e0 = glm::length(vb - va); // edge opposite vc
			double e1 = glm::length(vc - vb); // edge opposite va
			double e2 = glm::length(va - vc); // edge opposite vb
			double maxEdge = std::max(e0, std::max(e1, e2));
			if (maxEdge < 1e-15) continue;

			double minAlt = 2.0 * area / maxEdge;
			if (minAlt >= SLIVER_ALTITUDE_THRESHOLD) continue;

			// This is an altitude-based sliver. Find tip vertex and project onto opposite edge.
			Vec tip, edgeA, edgeB;
			if (maxEdge == e0) { tip = vc; edgeA = va; edgeB = vb; }
			else if (maxEdge == e1) { tip = va; edgeA = vb; edgeB = vc; }
			else { tip = vb; edgeA = vc; edgeB = va; }

			Vec edgeDir = edgeB - edgeA;
			double edgeLenSq = glm::dot(edgeDir, edgeDir);
			double t = glm::dot(tip - edgeA, edgeDir) / edgeLenSq;
			if (t < 0.0) t = 0.0;
			else if (t > 1.0) t = 1.0;
			Vec projected = edgeA + t * edgeDir;

			snapEntries.push_back({ tip, projected });
		}

		// Apply snap map + area filter during rebuild
		auto applySnap = [&](Vec& v) {
			for (const auto& entry : snapEntries) {
				Vec d = v - entry.source;
				if (glm::dot(d, d) < snapCellSizeSq) {
					v = entry.target;
					return;
				}
			}
			};

		bool hasAreaSlivers = false;
		if (snapEntries.empty()) {
			// No altitude-based slivers -- check for area-based slivers only
			for (uint32_t i = 0; i < n; i++) {
				Face f = workingMesh.GetFace(i);
				double area = areaOfTriangle(workingMesh.GetPoint(f.i0), workingMesh.GetPoint(f.i1), workingMesh.GetPoint(f.i2));
				if (area < SLIVER_AREA_THRESHOLD) {
					hasAreaSlivers = true;
					break;
				}
			}
		}

		if (snapEntries.empty() && !hasAreaSlivers) {
			meshInfoResult = meshInfoInput;
			return;
		}
		Geometry tmp;
		tmp.planes = workingMesh.planes;
		tmp.hasPlanes = workingMesh.hasPlanes;
		tmp.data = workingMesh.data;
		for (uint32_t i = 0; i < n; i++) {
			Face f = workingMesh.GetFace(i);
			Vec a = workingMesh.GetPoint(f.i0);
			Vec b = workingMesh.GetPoint(f.i1);
			Vec c = workingMesh.GetPoint(f.i2);

			// Apply snap map to vertices
			if (!snapEntries.empty()) {
				applySnap(a);
				applySnap(b);
				applySnap(c);
			}

			// Skip degenerate faces (area check covers both original slivers and collapsed altitude-slivers)
			if (areaOfTriangle(a, b, c) >= SLIVER_AREA_THRESHOLD) {
				tmp.AddFace(a, b, c, f.pId);
				tmp.hasPlanes = false; // enforce re-compute of planes
			}
		}

		auto meshInfoRemovedSlivers = meshCleanup::isMeshWatertight(tmp);

		if (meshInfoRemovedSlivers.numOpenEdges < meshInfoInput.numOpenEdges) {
			// improvement found. TODO: check if this condition is sufficient, or if other conditions make it safer
			workingMesh = tmp;
			meshInfoResult = meshInfoRemovedSlivers;
		}
		else {
			meshInfoResult = meshInfoInput;
#ifdef DUMP_CSG_MESHES
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
			webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup"+step+"-fail.obj");
#endif
		}
	}

	// ---------------------------------------------------------------
	// Resolve T-junctions: find boundary vertices that lie on the surface
	// of another face and re-triangulate that face to incorporate them.
	// This closes open edges caused by boolean operations where one face's
	// edge lies on another face's surface without topological connection.
	// ---------------------------------------------------------------
	static void ResolveTJunctions(Geometry& geom, std::string step, 
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
		const uint32_t nFaces = geom.numFaces;
		if (nFaces == 0) {
			meshInfoResult = meshInfoInput;
			return;
		}

		if (meshInfoInput.numOpenEdges == 0) {
			meshInfoResult = meshInfoInput;
			return;
		}

		// -- vertex deduplication (same spatial-hash pattern) --
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
		std::vector<Vec> canonPos;
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

		// Build per-face canonical vertex ids and cache face info
		struct FaceInfo {
			std::array<uint32_t, 3> vid;
			Vec normal;
			Vec verts[3]; // raw 3D positions
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
			faces[i].verts[0] = a;
			faces[i].verts[1] = b;
			faces[i].verts[2] = c;
			faces[i].pId = static_cast<uint32_t>(f.pId);
		}

		// -- Build edge map and find boundary edges --
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
		std::unordered_map<EKey, uint32_t, EKeyHash> edgeCount;
		for (uint32_t i = 0; i < nFaces; i++) {
			auto& v = faces[i].vid;
			for (int e = 0; e < 3; e++) {
				uint32_t va = v[e], vb = v[(e + 1) % 3];
				EKey ek = { std::min(va, vb), std::max(va, vb) };
				edgeCount[ek]++;
			}
		}

		// Collect boundary vertex IDs (vertices on edges with count == 1)
		std::unordered_set<uint32_t> boundaryVids;
		for (auto& [ek, cnt] : edgeCount) {
			if (cnt == 1) {
				boundaryVids.insert(ek.v0);
				boundaryVids.insert(ek.v1);
			}
		}

		if (boundaryVids.empty()) return;

		// -- For each boundary vertex, find faces it lies on --
		// Map: faceIdx -> list of {canonical vid, 3D position} to insert
		struct InsertionPoint {
			uint32_t vid;
			Vec pos;
		};
		std::unordered_map<uint32_t, std::vector<InsertionPoint>> faceInsertions;

		// Helper: squared distance from 2D point to segment, returns parameter t
		auto pointToSegDistSq2D = [](const std::array<double, 2>& p,
									 const std::array<double, 2>& a,
									 const std::array<double, 2>& b,
									 double& t) -> double {
			double dx = b[0] - a[0], dy = b[1] - a[1];
			double lenSq = dx * dx + dy * dy;
			if (lenSq < 1e-30) {
				t = 0.0;
				double cx = a[0] - p[0], cy = a[1] - p[1];
				return cx * cx + cy * cy;
			}
			t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / lenSq;
			if (t < 0.0) t = 0.0;
			else if (t > 1.0) t = 1.0;
			double cx = a[0] + t * dx - p[0];
			double cy = a[1] + t * dy - p[1];
			return cx * cx + cy * cy;
		};

		// Helper: check if 2D point is inside triangle using barycentric coords
		auto pointInTriangle2D = [](const std::array<double, 2>& p,
									const std::array<double, 2>& a,
									const std::array<double, 2>& b,
									const std::array<double, 2>& c) -> bool {
			double v0x = c[0] - a[0], v0y = c[1] - a[1];
			double v1x = b[0] - a[0], v1y = b[1] - a[1];
			double v2x = p[0] - a[0], v2y = p[1] - a[1];
			double dot00 = v0x * v0x + v0y * v0y;
			double dot01 = v0x * v1x + v0y * v1y;
			double dot02 = v0x * v2x + v0y * v2y;
			double dot11 = v1x * v1x + v1y * v1y;
			double dot12 = v1x * v2x + v1y * v2y;
			double inv = dot00 * dot11 - dot01 * dot01;
			if (std::abs(inv) < 1e-30) return false;
			inv = 1.0 / inv;
			double u = (dot11 * dot02 - dot01 * dot12) * inv;
			double v = (dot00 * dot12 - dot01 * dot02) * inv;
			return (u >= -1e-6) && (v >= -1e-6) && (u + v <= 1.0 + 1e-6);
		};

		const double edgeTolSq = cellSize * cellSize;

		for (uint32_t bvid : boundaryVids) {
			Vec bpos = canonPos[bvid];

			for (uint32_t fi = 0; fi < nFaces; fi++) {
				// Skip if this vertex is already a vertex of the face
				if (faces[fi].vid[0] == bvid || faces[fi].vid[1] == bvid || faces[fi].vid[2] == bvid)
					continue;

				// Skip degenerate faces
				if (glm::length(faces[fi].normal) < 0.5) continue;

				// Plane distance check
				double planeDist = std::abs(glm::dot(bpos - faces[fi].verts[0], faces[fi].normal));
				if (planeDist > cellSize) continue;

				// Project to 2D
				int dropAxis = SelectProjectionAxis(faces[fi].normal);
				auto p2d = ProjectTo2D(bpos, dropAxis);
				std::array<std::array<double, 2>, 3> tri2d = {
					ProjectTo2D(faces[fi].verts[0], dropAxis),
					ProjectTo2D(faces[fi].verts[1], dropAxis),
					ProjectTo2D(faces[fi].verts[2], dropAxis)
				};

				// Check if point is on any edge of the triangle
				bool onFace = false;
				for (int e = 0; e < 3; e++) {
					double t;
					double distSq = pointToSegDistSq2D(p2d, tri2d[e], tri2d[(e + 1) % 3], t);
					if (distSq < edgeTolSq && t > 1e-4 && t < (1.0 - 1e-4)) {
						onFace = true;
						break;
					}
				}

				// If not on an edge, check if inside the triangle
				if (!onFace) {
					onFace = pointInTriangle2D(p2d, tri2d[0], tri2d[1], tri2d[2]);
				}

				if (onFace) {
					// Verify this vertex is not too close to any existing face vertex
					bool tooClose = false;
					for (int v = 0; v < 3; v++) {
						Vec d = bpos - faces[fi].verts[v];
						if (glm::dot(d, d) < cellSizeSq) {
							tooClose = true;
							break;
						}
					}
					if (!tooClose) {
						faceInsertions[fi].push_back({ bvid, bpos });
					}
				}
			}
		}

		if (faceInsertions.empty()) return;

		// Safety: if too many insertions, skip
		uint32_t totalInsertions = 0;
		for (auto& [fi, pts] : faceInsertions) totalInsertions += static_cast<uint32_t>(pts.size());
		if (totalInsertions > nFaces) return;

		// -- Re-triangulate affected faces and rebuild geometry --
		Geometry rebuilt;
		rebuilt.planes = geom.planes;
		rebuilt.hasPlanes = geom.hasPlanes;
		rebuilt.data = geom.data;

		for (uint32_t fi = 0; fi < nFaces; fi++) {
			auto insertIt = faceInsertions.find(fi);
			if (insertIt == faceInsertions.end()) {
				// No insertions -- copy face as-is
				rebuilt.AddFace(faces[fi].verts[0], faces[fi].verts[1], faces[fi].verts[2], faces[fi].pId);
				continue;
			}

			const auto& insertions = insertIt->second;
			const Vec& faceNormal = faces[fi].normal;
			int dropAxis = SelectProjectionAxis(faceNormal);

			// Build CDT vertices: triangle corners (0,1,2) + insertion points (3,...)
			std::vector<CDT::V2d<double>> cdtVerts;
			std::vector<Vec> verts3D;
			cdtVerts.reserve(3 + insertions.size());
			verts3D.reserve(3 + insertions.size());

			for (int v = 0; v < 3; v++) {
				auto p2d = ProjectTo2D(faces[fi].verts[v], dropAxis);
				cdtVerts.push_back({ p2d[0], p2d[1] });
				verts3D.push_back(faces[fi].verts[v]);
			}
			for (auto& ins : insertions) {
				auto p2d = ProjectTo2D(ins.pos, dropAxis);
				cdtVerts.push_back({ p2d[0], p2d[1] });
				verts3D.push_back(ins.pos);
			}

			// Constraint edges: triangle boundary
			std::vector<CDT::Edge> cdtEdges;
			cdtEdges.push_back(CDT::Edge(0, 1));
			cdtEdges.push_back(CDT::Edge(1, 2));
			cdtEdges.push_back(CDT::Edge(2, 0));

			try {
				CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::AsProvided);
				cdt.insertVertices(cdtVerts);
				cdt.insertEdges(cdtEdges);
				cdt.eraseSuperTriangle();

				// Project original triangle corners to 2D for containment test
				std::vector<std::array<double, 2>> triRing = {
					{ cdtVerts[0].x, cdtVerts[0].y },
					{ cdtVerts[1].x, cdtVerts[1].y },
					{ cdtVerts[2].x, cdtVerts[2].y }
				};

				bool anyAdded = false;
				for (const auto& tri : cdt.triangles) {
					// Get centroid of CDT triangle
					double cx = 0, cy = 0;
					for (int v = 0; v < 3; v++) {
						cx += cdt.vertices[tri.vertices[v]].x;
						cy += cdt.vertices[tri.vertices[v]].y;
					}
					cx /= 3.0;
					cy /= 3.0;

					// Keep only triangles whose centroid is inside the original triangle
					std::array<double, 2> centroid = { cx, cy };
					if (!PointInPolygon2D(centroid, triRing)) continue;

					// Map CDT vertex indices back to 3D
					Vec va = verts3D[tri.vertices[0]];
					Vec vb = verts3D[tri.vertices[1]];
					Vec vc = verts3D[tri.vertices[2]];

					// Check winding consistency with original face normal
					Vec cr = glm::cross(vb - va, vc - va);
					if (glm::dot(cr, faceNormal) < 0) {
						std::swap(vb, vc);
					}

					Vec n = faceNormal;
					rebuilt.AddPoint(va, n);
					rebuilt.AddPoint(vb, n);
					rebuilt.AddPoint(vc, n);
					rebuilt.AddFace(rebuilt.numPoints - 3, rebuilt.numPoints - 2, rebuilt.numPoints - 1,
						faces[fi].pId);
					anyAdded = true;
				}

				if (!anyAdded) {
					// CDT produced no valid triangles -- keep original face
					rebuilt.AddFace(faces[fi].verts[0], faces[fi].verts[1], faces[fi].verts[2], faces[fi].pId);
				}
			}
			catch (...) {
				// CDT failed -- keep original face unmodified
				rebuilt.AddFace(faces[fi].verts[0], faces[fi].verts[1], faces[fi].verts[2], faces[fi].pId);
			}
		}

		MeshWatertightInfo infoRebuilt = meshCleanup::isMeshWatertight(rebuilt);
		if (infoRebuilt.numOpenEdges < meshInfoInput.numOpenEdges) {
			geom = std::move(rebuilt);
			meshInfoResult = infoRebuilt;
		}
		else {
			meshInfoResult = meshInfoInput;
#ifdef _DEBUG
			webifc::geometry::IfcGeometry webifcGeom = webifc::geometry::booleanManager::convertToWebIfc(geom);
			webifc::geometry::IfcGeometry geomFail = webifc::geometry::booleanManager::convertToWebIfc(rebuilt);
			webifc::io::DumpIfcGeometry(webifcGeom, "meshCleanup" + step + "-input.obj");
			webifc::io::DumpIfcGeometry(geomFail, "meshCleanup" + step + "-fail.obj");
#endif
		}
	}

	// ---------------------------------------------------------------
	// Patch coplanar holes: find boundary-edge loops and fill them with earcut triangulation when all loop vertices are coplanar.
	// ---------------------------------------------------------------
	static void PatchCoplanarHoles(Geometry& geom, std::string step, 
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
		if (meshInfoInput.numOpenEdges == 0) {
			meshInfoResult = meshInfoInput;
			return;
		}
		const uint32_t nFaces = geom.numFaces;
		if (nFaces == 0) return;

		// Save backup for revert if patching doesn't improve the mesh
		Geometry backup = geom;

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

		// -- Build edge -> face adjacency --
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
		auto makeEKey = [](uint32_t a, uint32_t b) -> EKey {
			return { std::min(a, b), std::max(a, b) };
		};
		std::unordered_map<EKey, std::vector<uint32_t>, EKeyHash> edgeFaces;
		for (uint32_t i = 0; i < nFaces; i++) {
			auto& v = faces[i].vid;
			for (int e = 0; e < 3; e++) {
				edgeFaces[makeEKey(v[e], v[(e + 1) % 3])].push_back(i);
			}
		}

		// Collect boundary edges and build vertexAdjFace for raw position lookup
		std::unordered_map<uint32_t, uint32_t> vertexAdjFace;
		bool hasBoundaryEdges = false;
		for (auto& [ek, fl] : edgeFaces) {
			if (fl.size() == 1) {
				hasBoundaryEdges = true;
				vertexAdjFace[ek.v0] = fl[0];
				vertexAdjFace[ek.v1] = fl[0];
			}
		}

		if (!hasBoundaryEdges) return; // mesh is already closed

		// -- Face-fan walk: given boundary edge prev->cur on face prevFace,
		//    find the next boundary edge cur->next by walking around cur
		//    through the face fan. Returns {next vertex, face of edge cur->next}.
		auto findNextBoundary = [&](uint32_t cur, uint32_t prev, uint32_t prevFace,
								   uint32_t& outNext, uint32_t& outFace) -> bool {
			// Find the "other" vertex of prevFace (the one that is neither prev nor cur)
			uint32_t other = UINT32_MAX;
			for (int j = 0; j < 3; j++) {
				uint32_t v = faces[prevFace].vid[j];
				if (v != prev && v != cur) { other = v; break; }
			}
			if (other == UINT32_MAX) return false;

			// Walk the face fan around cur: start from edge cur-other
			uint32_t walkEdgeOther = other;
			uint32_t walkFace = prevFace;
			for (uint32_t safety = 0; safety < nFaces; safety++) {
				EKey ek = makeEKey(cur, walkEdgeOther);
				auto it = edgeFaces.find(ek);
				if (it == edgeFaces.end()) return false;

				auto& fl = it->second;
				if (fl.size() == 1) {
					// Boundary edge found -- this is the next edge in the loop
					outNext = walkEdgeOther;
					outFace = fl[0];
					return true;
				}

				// Internal edge (size >= 2): cross to the other face
				uint32_t nextFace = UINT32_MAX;
				for (uint32_t fi : fl) {
					if (fi != walkFace) { nextFace = fi; break; }
				}
				if (nextFace == UINT32_MAX) return false;

				// Find the next vertex in nextFace (the one that is neither cur nor walkEdgeOther)
				uint32_t nextOther = UINT32_MAX;
				for (int j = 0; j < 3; j++) {
					uint32_t v = faces[nextFace].vid[j];
					if (v != cur && v != walkEdgeOther) { nextOther = v; break; }
				}
				if (nextOther == UINT32_MAX) return false;

				walkFace = nextFace;
				walkEdgeOther = nextOther;
			}
			return false; // fan walk exceeded safety limit
		};

		// -- Trace boundary loops using face-fan walk --
		// Track visited directed boundary edges (as prev<<32|cur) to handle
		// junction vertices where two loops share a vertex.
		std::unordered_set<uint64_t> visitedEdges;
		auto dirEdgeKey = [](uint32_t from, uint32_t to) -> uint64_t {
			return (static_cast<uint64_t>(from) << 32) | static_cast<uint64_t>(to);
		};

		std::vector<std::vector<uint32_t>> loops;
		std::vector<uint32_t> loopAdjacentFace;

		// Find all boundary edges to use as potential starting edges
		struct BoundaryEdge { uint32_t v0, v1; uint32_t faceIdx; };
		std::vector<BoundaryEdge> boundaryEdges;
		for (auto& [ek, fl] : edgeFaces) {
			if (fl.size() == 1) {
				boundaryEdges.push_back({ ek.v0, ek.v1, fl[0] });
			}
		}

		for (auto& be : boundaryEdges) {
			// Try starting from this boundary edge in both directions
			for (int dir = 0; dir < 2; dir++) {
				uint32_t startPrev = dir == 0 ? be.v0 : be.v1;
				uint32_t startCur = dir == 0 ? be.v1 : be.v0;

				if (visitedEdges.count(dirEdgeKey(startPrev, startCur))) continue;

				std::vector<uint32_t> loop;
				uint32_t prev = startPrev;
				uint32_t cur = startCur;
				uint32_t curFace = be.faceIdx;
				uint32_t firstFace = be.faceIdx;
				bool valid = true;

				loop.push_back(startPrev);

				while (true) {
					visitedEdges.insert(dirEdgeKey(prev, cur));
					loop.push_back(cur);

					if (cur == startPrev && prev == loop[loop.size() - 2]) {
						// Check: did we return to the start edge?
						// We pushed startPrev at the beginning. If cur == startPrev,
						// we've closed the loop.
						break;
					}

					uint32_t next, nextFace;
					if (!findNextBoundary(cur, prev, curFace, next, nextFace)) {
						valid = false;
						break;
					}

					prev = cur;
					cur = next;
					curFace = nextFace;

					// Check if we've returned to starting vertex to close the loop
					if (cur == startPrev) {
						visitedEdges.insert(dirEdgeKey(prev, cur));
						break;
					}

					if (visitedEdges.count(dirEdgeKey(prev, cur))) {
						valid = false;
						break;
					}

					if (loop.size() > canonPos.size() + 1) {
						valid = false;
						break;
					}
				}

				if (valid && loop.size() >= 3) {
					loops.push_back(std::move(loop));
					loopAdjacentFace.push_back(firstFace);
				}
			}
		}

		// Deduplicate loops: the tracer tries both directions for each starting
		// edge, so each physical loop is traced forward and reverse. Remove
		// duplicates by comparing vertex sets.
		{
			std::vector<bool> isDuplicate(loops.size(), false);
			for (size_t i = 0; i < loops.size(); i++) {
				if (isDuplicate[i]) continue;
				std::unordered_set<uint32_t> setI(loops[i].begin(), loops[i].end());
				for (size_t j = i + 1; j < loops.size(); j++) {
					if (isDuplicate[j]) continue;
					if (loops[j].size() != loops[i].size()) continue;
					std::unordered_set<uint32_t> setJ(loops[j].begin(), loops[j].end());
					if (setI == setJ) {
						isDuplicate[j] = true;
					}
				}
			}
			std::vector<std::vector<uint32_t>> uniqueLoops;
			std::vector<uint32_t> uniqueLoopFaces;
			for (size_t i = 0; i < loops.size(); i++) {
				if (!isDuplicate[i]) {
					uniqueLoops.push_back(std::move(loops[i]));
					uniqueLoopFaces.push_back(loopAdjacentFace[i]);
				}
			}
			loops = std::move(uniqueLoops);
			loopAdjacentFace = std::move(uniqueLoopFaces);
		}

		if (loops.empty()) {
			meshInfoResult = meshInfoInput;
			return;
		}

		const double PLANE_EPS = 1e-5;
		struct LoopInfo {
			std::vector<uint32_t> vid;
			Vec planeNormal;
			Vec planePoint;
			Vec adjNormal;
			uint32_t adjPlaneId;
		};

		std::unordered_map<uint32_t, Vec> rawBoundaryPos;
		rawBoundaryPos.reserve(vertexAdjFace.size());

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

		if (validLoops.empty()) {
			meshInfoResult = meshInfoInput;
			return;
		}

		// -- Global parent-child containment detection across ALL loops --
		// For each candidate parent, project both loops using the PARENT's normal
		// (the correct axis for that parent's plane). This handles loops with
		// slightly different normals that a single global axis would miss.
		struct ProjectedLoop {
			size_t loopIndex;
			int dropAxis;  // per-loop projection axis
			std::vector<std::array<double, 2>> projected;
			double area = 0.0;
			int parent = -1;
			int depth = 0;
		};

		std::vector<ProjectedLoop> projectedLoops;
		projectedLoops.reserve(validLoops.size());

		for (size_t i = 0; i < validLoops.size(); i++) {
			ProjectedLoop pl;
			pl.loopIndex = i;
			pl.dropAxis = SelectProjectionAxis(validLoops[i].planeNormal);
			pl.projected.reserve(validLoops[i].vid.size());
			for (uint32_t vid : validLoops[i].vid) {
				pl.projected.push_back(ProjectTo2D(canonPos[vid], pl.dropAxis));
			}
			pl.area = std::abs(SignedArea2D(pl.projected));
			if (pl.area > 1e-12) {
				projectedLoops.push_back(std::move(pl));
			}
		}

		// Find parent (smallest enclosing loop) for each loop.
		// For each candidate parent, re-project the child's sample point using
		// the PARENT's axis and test containment against the parent's polygon.
		for (size_t loopIdx = 0; loopIdx < projectedLoops.size(); loopIdx++) {
			double bestParentArea = std::numeric_limits<double>::max();

			for (size_t candidateIdx = 0; candidateIdx < projectedLoops.size(); candidateIdx++) {
				if (candidateIdx == loopIdx) continue;
				if (projectedLoops[candidateIdx].area <= projectedLoops[loopIdx].area + 1e-12) continue;

				// Project the child's first vertex using the CANDIDATE PARENT's axis
				int parentAxis = projectedLoops[candidateIdx].dropAxis;
				uint32_t sampleVid = validLoops[projectedLoops[loopIdx].loopIndex].vid[0];
				auto samplePoint = ProjectTo2D(canonPos[sampleVid], parentAxis);

				if (!PointInPolygon2D(samplePoint, projectedLoops[candidateIdx].projected)) continue;

				if (projectedLoops[candidateIdx].area < bestParentArea) {
					bestParentArea = projectedLoops[candidateIdx].area;
					projectedLoops[loopIdx].parent = static_cast<int>(candidateIdx);
				}
			}
		}

		// Compute nesting depth
		for (auto& pl : projectedLoops) {
			for (int p = pl.parent; p >= 0; p = projectedLoops[p].parent) {
				pl.depth++;
			}
		}

		// Triangulate: even-depth loops are outer boundaries, odd-depth are holes
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

				// Re-project hole using the OUTER loop's axis for consistent earcut input
				const auto& holeLoop = validLoops[projectedLoops[holeIdx].loopIndex];
				int outerAxis = projectedLoops[outerIdx].dropAxis;
				std::vector<std::array<double, 2>> holeProjected;
				holeProjected.reserve(holeLoop.vid.size());
				for (uint32_t vid : holeLoop.vid) {
					holeProjected.push_back(ProjectTo2D(canonPos[vid], outerAxis));
				}
				polygon.push_back(std::move(holeProjected));
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
		if (geom.numFaces > nFaces) {
			// faces have been added -> check if the result actually improved
			geom.hasPlanes = false;
			auto infoPatched = meshCleanup::isMeshWatertight(geom);

			if (infoPatched.numOpenEdges < meshInfoInput.numOpenEdges) {
				meshInfoResult = infoPatched;
#ifdef DUMP_CSG_MESHES
				webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(geom);
				webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup" + step + "-patched.obj");
#endif
			}
			else {
				// patching didn't improve -- revert
				geom = std::move(backup);
				meshInfoResult = meshInfoInput;
			}
		}
		else {
			meshInfoResult = meshInfoInput;
		}
	}

	void RemoveThinMembranes(Geometry& workingMesh, std::string step, 
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
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
					if (intersect_ray_triangle(fv[i].center, rayEnd, fv[faceIdx].a, fv[faceIdx].b,
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
		// Signed volume relative to each component's centroid is used to detect membranes (flat -> V ~= 0) vs solids (V >> 0).
		// B.2 marks are tracked separately (b2Marked) for independent removal pass.

		std::vector<bool> b2Marked(nFaces, false);

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

		uint32_t b2Count = 0;
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
			//   flat open surface -> V ~= 0 (points coplanar with centroid)
			std::vector<double> mVol(mNumComp, 0.0);
			for (uint32_t i = 0; i < nFaces; i++) {
				int ci = mCompId[i];
				Vec va = fv[i].a - mci[ci].centroid;
				Vec vb = fv[i].b - mci[ci].centroid;
				Vec vc = fv[i].c - mci[ci].centroid;
				mVol[ci] += glm::dot(va, glm::cross(vb, vc)) / 6.0;
			}

			for (uint32_t i = 0; i < nFaces; i++) {
				if (std::abs(mVol[mCompId[i]]) < VOLUME_THRESHOLD) {
					b2Marked[i] = true;
					b2Count++;
				}
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

		// B.5: boundary-edge peeling -- detect single-layer membrane flaps.
		// A face with >= 2 boundary edges is a "tip" of a membrane. Iteratively
		// peel tip faces: each removal exposes new boundary edges on adjacent
		// faces, which may become new tips. This catches flat membranes that
		// B.1 misses (not double-layer) and B.2 misses (connected to main body).
		{
			// Count boundary edges per face
			auto countBoundaryEdges = [&](uint32_t fi, const std::vector<bool>& removed) -> int {
				int count = 0;
				for (int e = 0; e < 3; e++) {
					uint32_t va = fvid[fi][e], vb = fvid[fi][(e + 1) % 3];
					EKey ek = { std::min(va, vb), std::max(va, vb) };
					auto it = edgeFaces.find(ek);
					if (it == edgeFaces.end()) continue;
					// Count non-removed faces sharing this edge
					int alive = 0;
					for (uint32_t adj : it->second) {
						if (!removed[adj]) alive++;
					}
					// If only this face uses the edge, it's a boundary edge
					if (alive <= 1) count++;
				}
				return count;
			};

			std::vector<bool> peeled(nFaces, false);
			for (int iter = 0; iter < static_cast<int>(nFaces); iter++) {
				bool changed = false;
				for (uint32_t i = 0; i < nFaces; i++) {
					if (peeled[i] || thinMarked[i]) continue;
					if (countBoundaryEdges(i, peeled) >= 2) {
						peeled[i] = true;
						changed = true;
					}
				}
				if (!changed) break;
			}

			// Safety: only accept if peeled region is < 50% of faces
			uint32_t peelCount = 0;
			for (uint32_t i = 0; i < nFaces; i++)
				if (peeled[i]) peelCount++;

			if (peelCount > 0 && peelCount < nFaces / 2) {
				for (uint32_t i = 0; i < nFaces; i++)
					if (peeled[i]) thinMarked[i] = true;
			}
		}

		// B.6: safety - if > 75 % of faces are thin, something went wrong (e.g. very thin but valid geometry); disable.
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

		// -- Three-pass rebuild: each pass has its own improvement check --
		// Pass 1 (high certainty): remove bad components (Phase C: disconnected non-manifold zero-volume fragments)
		// Pass 2 (medium certainty): also remove B.2 faces (zero-volume manifold-edge components)
		// Pass 3 (lower certainty): also remove B.1+B.3 thin-marked faces (double-layer + erosion)

		MeshWatertightInfo currentInfo = meshInfoInput;

		// Pass 1: remove bad components only (Phase C)
		if (removedByComp > 0 && removedByComp < nFaces) {
			Geometry pass1;
			pass1.planes = workingMesh.planes;
			pass1.hasPlanes = workingMesh.hasPlanes;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (badComp[compId[i]]) continue;
				pass1.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
			}
			pass1.data = workingMesh.data;

			auto infoPass1 = meshCleanup::isMeshWatertight(pass1);
			if (infoPass1.numOpenEdges < currentInfo.numOpenEdges) {
				workingMesh = std::move(pass1);
				currentInfo = infoPass1;
			}
		}

		// Pass 2: also remove B.2 zero-volume manifold-edge components
		// Removing membrane faces may create new boundary edges where the membrane
		// connected to the main body via non-manifold edges. Patch those holes
		// before checking improvement.
		if (b2Count > 0) {
			Geometry pass2;
			pass2.planes = workingMesh.planes;
			pass2.hasPlanes = workingMesh.hasPlanes;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (badComp[compId[i]]) continue;
				if (b2Marked[i]) continue;
				pass2.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
			}
			pass2.data = workingMesh.data;

			MeshWatertightInfo infoPass2Before = meshCleanup::isMeshWatertight(pass2);
			MeshWatertightInfo infoPass2 = infoPass2Before;
			if (infoPass2Before.numOpenEdges > 0) {
				PatchCoplanarHoles(pass2, step + "b", infoPass2Before, infoPass2);
			}
			if (infoPass2.numOpenEdges < currentInfo.numOpenEdges) {
				workingMesh = std::move(pass2);
				currentInfo = infoPass2;
			}
		}

		// Pass 3: also remove B.1+B.3+B.5 thin-marked faces
		if (thinEnabled && thinCount > 0) {
			Geometry pass3;
			pass3.planes = workingMesh.planes;
			pass3.hasPlanes = workingMesh.hasPlanes;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (badComp[compId[i]]) continue;
				if (b2Marked[i]) continue;
				if (thinMarked[i]) continue;
				pass3.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
			}
			pass3.data = workingMesh.data;

			MeshWatertightInfo infoPass3Before = meshCleanup::isMeshWatertight(pass3);
			MeshWatertightInfo infoPass3 = infoPass3Before;
			if (infoPass3Before.numOpenEdges > 0) {
				PatchCoplanarHoles(pass3, step + "c", infoPass3Before, infoPass3);
			}
			if (infoPass3.numOpenEdges < currentInfo.numOpenEdges) {
				workingMesh = std::move(pass3);
				currentInfo = infoPass3;
			}
		}

		meshInfoResult = currentInfo;

#ifdef _DEBUG
		if (currentInfo.numOpenEdges > 0 && currentInfo.numOpenEdges >= meshInfoInput.numOpenEdges) {
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
			webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup" + step + "-fail.obj");
		}
#endif
	}

	void removeTempFiles() {
		const std::filesystem::path dir = R"(E:\work\creoox\cxconverter)";

		for (const auto& entry : std::filesystem::directory_iterator(dir)) {
			if (!entry.is_regular_file())
				continue;

			const auto& path = entry.path();
			const std::string filename = path.filename().string();

			// match: meshCleanup*.obj
			if (filename.rfind("meshCleanup", 0) == 0 && path.extension() == ".obj") {
				try {
					std::filesystem::remove(path);
				}
				catch (const std::exception& e) {
					std::cerr << "Failed to delete: " << path << " (" << e.what() << ")\n";
				}
			}
		}
	}

	// Post-boolean cleanup
	void PostBooleanOperationMeshCleanup(fuzzybools::Geometry& input) {
		if (input.numFaces > 8000) {
#ifdef _DEBUG
			std::cout << "PostBooleanOperationMeshCleanup: skipping mesh with " << input.numFaces << " faces" << std::endl;
#endif
			return;
		}
		fuzzybools::Geometry workingMesh = input;
		auto meshInfoOnEntry = meshCleanup::isMeshWatertight(input);
		if (meshInfoOnEntry.watertight) {
			return;
		}
#ifdef DUMP_CSG_MESHES
		// remove all E:\work\creoox\cxconverter\meshCleanup*.obj
		removeTempFiles();

		if (meshInfoOnEntry.numOpenEdges == 4) {
			int wait = 0;
		}
		
		webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(input);
		webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup1.obj");
#endif

		// 1: remove degenerated triangles
		MeshWatertightInfo meshInfoBeforeRemoveDegeneratedTriangles = meshInfoOnEntry;
		MeshWatertightInfo meshInfoAfterRemoveDegeneratedTriangles = meshInfoOnEntry;
		removeDegeneratedTriangles(workingMesh, "1", meshInfoBeforeRemoveDegeneratedTriangles, meshInfoAfterRemoveDegeneratedTriangles);
		
		// 2: remove disconnected fragments (Phase C only -- safe before hole patching)
		MeshWatertightInfo meshInfoAfterRemoveFragments = meshInfoAfterRemoveDegeneratedTriangles;
		RemoveDisconnectedFragments(workingMesh, "2", meshInfoAfterRemoveDegeneratedTriangles, meshInfoAfterRemoveFragments);

		// 3: patch coplanar holes
		// remove thin membranes later, since a (closable) loop of open edges can cause false membrane detection,
		// membrane removal increases open edges, result gets reverted, nothing gets fixed
		MeshWatertightInfo meshInfoBeforePatchCoplanarHoles = meshInfoAfterRemoveFragments;
		MeshWatertightInfo meshInfoAfterPatchCoplanarHoles = meshInfoAfterRemoveFragments;
		PatchCoplanarHoles(workingMesh, "3", meshInfoBeforePatchCoplanarHoles, meshInfoAfterPatchCoplanarHoles);

		// 4: remove thin membranes
		MeshWatertightInfo meshInfoBeforeRemoveThinMembranes = meshInfoAfterPatchCoplanarHoles;
		MeshWatertightInfo meshInfoAfterRemoveThinMembranes = meshInfoAfterPatchCoplanarHoles;
		RemoveThinMembranes(workingMesh, "4", meshInfoBeforeRemoveThinMembranes, meshInfoAfterRemoveThinMembranes);

		// 5: resolve T-junctions
		MeshWatertightInfo meshInfoBeforeResolveTJunctions = meshInfoAfterRemoveThinMembranes;
		MeshWatertightInfo meshInfoAfterResolveTJunctions = meshInfoAfterRemoveThinMembranes;
		ResolveTJunctions(workingMesh, "5", meshInfoBeforeResolveTJunctions, meshInfoAfterResolveTJunctions);

		// 6: PatchCoplanarHoles re-run
		MeshWatertightInfo meshInfoBeforePatchCoplanarHoles2 = meshInfoAfterResolveTJunctions;
		MeshWatertightInfo meshInfoAfterPatchCoplanarHoles2 = meshInfoAfterResolveTJunctions;
		PatchCoplanarHoles(workingMesh, "6", meshInfoBeforePatchCoplanarHoles2, meshInfoAfterPatchCoplanarHoles2);

		auto meshInfoOnExit = meshInfoAfterPatchCoplanarHoles2;
		if(meshInfoOnEntry.numOpenEdges > meshInfoOnExit.numOpenEdges){
			input = std::move(workingMesh);
			input.hasPlanes = false;

#ifdef DUMP_CSG_MESHES
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(input);
			webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup7-improved.obj");
#endif
		}
	}
}
