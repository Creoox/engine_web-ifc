#include <queue>
#include <array>
#include <tuple>
#include <filesystem>
#include <functional>
#include <algorithm>
#include <limits>
#include <cmath>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include <mapbox/earcut.hpp>
#include <CDT.h>

#include <web-ifc/geometry/IfcGeometryProcessor.h>
#include <web-ifc/geometry/operations/bim-geometry/booleanUtils.h>
#include <web-ifc/geometry/operations/bim-geometry/geometry.h>
#include <web-ifc/geometry/operations/boolean-utils/clip-mesh.h>
#include <web-ifc/geometry/operations/boolean-utils/shared-position.h>
#include "meshCleanup.h"   // keep it like this

#include "cpp/test/io_helpers.h"

//#if defined(_DEBUG)
#define DUMP_CSG_MESHES
//#endif

using namespace fuzzybools;

static std::filesystem::path g_debugDumpDirectory;

// MeshPenaltyScore and MeshPenaltyImprovedNoRegression are now
// inline in meshCleanup.h so they can be used from other translation units.
using meshCleanup::MeshPenaltyScore;
using meshCleanup::MeshPenaltyImprovedNoRegression;

void meshCleanup::SetDebugDumpDirectory(const std::filesystem::path& dir) {
	g_debugDumpDirectory = dir;
}

static std::filesystem::path BuildDebugDumpPath(const std::string& filename) {
	if (g_debugDumpDirectory.empty()) {
		return std::filesystem::path(filename);
	}
	return g_debugDumpDirectory / filename;
}

void meshCleanup::DumpDebugGeometry(const fuzzybools::Geometry& geom, const std::string& filename) {
	if (g_debugDumpDirectory.empty()) return;
	std::error_code ec;
	std::filesystem::create_directories(g_debugDumpDirectory, ec);
	webifc::geometry::IfcGeometry webifcGeom = webifc::geometry::booleanManager::convertToWebIfc(geom);
	std::string filenameWithPath = BuildDebugDumpPath(filename).string();
	webifc::io::DumpIfcGeometry(webifcGeom, filenameWithPath);
}

static void meshCleanup::DumpDebugGeometry(const webifc::geometry::IfcGeometry& geom, const std::string& filename) {
	if (!g_debugDumpDirectory.empty()) {
		std::error_code ec;
		std::filesystem::create_directories(g_debugDumpDirectory, ec);
	}
	std::string filenameWithPath = BuildDebugDumpPath(filename).string();
	webifc::io::DumpIfcGeometry(geom, filenameWithPath);
}

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

meshCleanup::MeshInfo meshCleanup::isMeshWatertight(const fuzzybools::Geometry& geom) {
	MeshInfo info;
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

	// Dedup tolerance only for watertightness accounting. Must be *smaller*
	// than the shortest legitimate edge in typical IFC input; the general
	// geometry tolerance (toleranceVectorEquality = 1e-4) is too loose here,
	// it merges sub-0.1mm distinct vertices in real files (e.g.
	// IfcClosedShell #1165466 in test67.ifc has adjacent polyloop vertices
	// at 2.8e-5 / 3.4e-5 spacing that must NOT be merged or the mesh is
	// wrongly flagged as non-manifold).
	constexpr double watertightDedupTolerance = 1.0E-07;
	const double cellSize = watertightDedupTolerance;
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
		else if (count > 2) {
			++info.numNonManifoldEdges;
		}
	}

	info.watertight = (info.numOpenEdges == 0 && info.numNonManifoldEdges == 0);
	return info;
}

// Remove disconnected non-manifold zero-volume components.
// Safe to run at any point -- purely topological, no heuristics.
uint32_t meshCleanup::RemoveDisconnectedFragments(Geometry& workingMesh, std::string step,
	const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult) {
	const Geometry baseMesh = workingMesh;
	const uint32_t nFaces = baseMesh.numFaces;
	if (nFaces == 0) { meshInfoResult = meshInfoInput; return 0; }

	constexpr double VOLUME_THRESHOLD = 1e-6;

	// Cache face vertices
	struct FV { Vec a, b, c; uint32_t pId; Vec center; };
	std::vector<FV> fv(nFaces);
	for (uint32_t i = 0; i < nFaces; i++) {
		Face f = baseMesh.GetFace(i);
		Vec a = baseMesh.GetPoint(f.i0);
		Vec b = baseMesh.GetPoint(f.i1);
		Vec c = baseMesh.GetPoint(f.i2);
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

	// Face adjacency + BFS components. Only traverse manifold connections so
	// sheets that touch a solid through a non-manifold seam can still be
	// identified and removed as separate zero-volume fragments.
	std::vector<std::vector<uint32_t>> faceAdj(nFaces);
	for (auto& [ek, fl] : edgeFaces) {
		if (fl.size() != 2) continue;
		for (size_t a = 0; a < fl.size(); a++)
			for (size_t b = a + 1; b < fl.size(); b++) {
				faceAdj[fl[a]].push_back(fl[b]);
				faceAdj[fl[b]].push_back(fl[a]);
			}
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

	if (numComp <= 1) { meshInfoResult = meshInfoInput; return 0; }

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
	for (auto& [ek, fl] : edgeFaces) {
		if (fl.size() == 2) continue;
		std::unordered_set<int> touchedComponents;
		for (uint32_t faceIdx : fl) {
			touchedComponents.insert(compId[faceIdx]);
		}
		for (int cid : touchedComponents) {
			comps[cid].boundaryEdges++;
		}
	}

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

	if (removedByComp == 0) {
		meshInfoResult = meshInfoInput;
		return 0;
	}
	if (removedByComp >= nFaces) {
		// All components are non-manifold + zero-volume. Compute total mesh
		// volume to decide: if the whole mesh is zero-volume (degenerate membrane
		// from a boolean subtract that fully consumed the first operand), allow
		// removal. Otherwise preserve the mesh.
		double totalVol = 0.0;
		for (uint32_t i = 0; i < nFaces; i++) {
			int cid = compId[i];
			Vec va = fv[i].a - comps[cid].centroid;
			Vec vb = fv[i].b - comps[cid].centroid;
			Vec vc = fv[i].c - comps[cid].centroid;
			totalVol += glm::dot(va, glm::cross(vb, vc)) / 6.0;
		}
		if (std::abs(totalVol) >= VOLUME_THRESHOLD) {
			meshInfoResult = meshInfoInput;
			return 0;
		}
		// Zero-volume: fall through to rebuild as empty mesh
	}

	Geometry cleaned;
	cleaned.planes = workingMesh.planes;
	cleaned.hasPlanes = workingMesh.hasPlanes;
	for (uint32_t i = 0; i < nFaces; i++) {
		if (badComp[compId[i]]) continue;
		cleaned.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
	}
	cleaned.data = workingMesh.data;

	// Empty mesh from zero-volume membrane removal: skip topology gate.
	if (cleaned.numFaces == 0) {
		workingMesh = std::move(cleaned);
		meshInfoResult = meshCleanup::isMeshWatertight(workingMesh);
		return removedByComp;
	}

	auto infoCleaned = meshCleanup::isMeshWatertight(cleaned);
	if (MeshPenaltyImprovedNoRegression(meshInfoInput, infoCleaned)) {
		workingMesh = std::move(cleaned);
		meshInfoResult = infoCleaned;
		return removedByComp;
	}
	else {
		meshInfoResult = meshInfoInput;
		return 0;
	}
}

// -- Phase A: strip sliver triangles -----------------------------
	// Threshold 1e-9 m^2 sits well above the toleranceAddFace filter (~5e-11 m^2) so it removes only absolute dregs, while staying
	// far below the smallest real feature in any meter-scale model.
	// Additionally, detect altitude-based slivers: triangles with very long edges but nearly collinear vertices
	// (minAltitude < toleranceVectorEquality). For these, snap the tip vertex onto the opposite edge
	// so that adjacent faces can later be split at the snap point (T-junction resolution).
uint32_t meshCleanup::RemoveDegeneratedTriangles(Geometry& workingMesh, std::string step,
	const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult) {
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
		return 0;
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
	const uint32_t removedFaces = n - tmp.numFaces;
	const uint32_t appliedChanges = removedFaces + static_cast<uint32_t>(snapEntries.size());

	if (TopologyStrictlyImprovedWithoutRegression(meshInfoInput, meshInfoRemovedSlivers)) {
		workingMesh = tmp;
		meshInfoResult = meshInfoRemovedSlivers;
		return appliedChanges;
	}
	else {
		meshInfoResult = meshInfoInput;
#ifdef DUMP_CSG_MESHES
		webifc::geometry::IfcGeometry  inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
		DumpDebugGeometry(inputWebIfc, "meshCleanup"+step+"-fail.obj");
#endif
		return 0;
	}
}

// ---------------------------------------------------------------
// Resolve T-junctions: find boundary vertices that lie on the surface
// of another face and re-triangulate that face to incorporate them.
// This closes open edges caused by boolean operations where one face's
// edge lies on another face's surface without topological connection.
// ---------------------------------------------------------------
uint32_t meshCleanup::ResolveTJunctions(Geometry& geom, std::string step,
	const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult) {
	const uint32_t nFaces = geom.numFaces;
	if (nFaces == 0) {
		meshInfoResult = meshInfoInput;
		return 0;
	}

	if (meshInfoInput.numOpenEdges == 0) {
		meshInfoResult = meshInfoInput;
		return 0;
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

	if (boundaryVids.empty()) {
		meshInfoResult = meshInfoInput;
		return 0;
	}

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

	if (faceInsertions.empty()) {
		meshInfoResult = meshInfoInput;
		return 0;
	}

	// Safety: if too many insertions, skip
	uint32_t totalInsertions = 0;
	for (auto& [fi, pts] : faceInsertions) totalInsertions += static_cast<uint32_t>(pts.size());
	if (totalInsertions > nFaces) {
		meshInfoResult = meshInfoInput;
		return 0;
	}

	// -- Re-triangulate affected faces and rebuild geometry --
	Geometry rebuilt;
	rebuilt.planes = geom.planes;
	rebuilt.hasPlanes = geom.hasPlanes;
	rebuilt.data = geom.data;
	uint32_t retriangulatedFaces = 0;

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
			else {
				++retriangulatedFaces;
			}
		}
		catch (...) {
			// CDT failed -- keep original face unmodified
			rebuilt.AddFace(faces[fi].verts[0], faces[fi].verts[1], faces[fi].verts[2], faces[fi].pId);
		}
	}

	MeshInfo infoRebuilt = meshCleanup::isMeshWatertight(rebuilt);
	if (MeshPenaltyImprovedNoRegression(meshInfoInput, infoRebuilt)) {
		geom = std::move(rebuilt);
		meshInfoResult = infoRebuilt;
		return retriangulatedFaces;
	}
	else {
		meshInfoResult = meshInfoInput;
#ifdef DUMP_CSG_MESHES
		webifc::geometry::IfcGeometry webifcGeom = webifc::geometry::booleanManager::convertToWebIfc(geom);
		webifc::geometry::IfcGeometry geomFail = webifc::geometry::booleanManager::convertToWebIfc(rebuilt);
		DumpDebugGeometry(webifcGeom, "meshCleanup" + step + "-input.obj");
		DumpDebugGeometry(geomFail, "meshCleanup" + step + "-fail.obj");
#endif
		return 0;
	}
}

uint32_t meshCleanup::RemoveTinyBoundaryBridgeFaces(Geometry& geom, std::string step,
	const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult) {
	(void)step;
	meshInfoResult = meshInfoInput;

	if (geom.numFaces == 0 || meshInfoInput.numOpenEdges == 0) {
		return 0;
	}

	const double cellSize = toleranceVectorEquality;
	const double cellSizeSq = cellSize * cellSize;
	const double areaThreshold = toleranceVectorEquality * toleranceVectorEquality;
	const double minAltitudeThreshold = toleranceVectorEquality * 0.25;

	using GridKey = std::tuple<int64_t, int64_t, int64_t>;
	struct GridKeyHash {
		size_t operator()(const GridKey& k) const {
			size_t h = std::hash<int64_t>()(std::get<0>(k));
			h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
			h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
			return h;
		}
	};

	auto getCell = [&](const Vec& p) -> GridKey {
		return { static_cast<int64_t>(std::floor(p.x / cellSize)),
				static_cast<int64_t>(std::floor(p.y / cellSize)),
				static_cast<int64_t>(std::floor(p.z / cellSize)) };
	};

	struct FaceInfo {
		std::array<uint32_t, 3> vid;
		Vec verts[3];
		uint32_t pId;
		double area;
		double maxEdge;
	};

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

	struct CandidateFace {
		uint32_t faceIndex;
		double area;
		double minAltitude;
	};

	Geometry current = geom;
	MeshInfo currentInfo = meshInfoInput;
	uint32_t removedFaces = 0;

	while (current.numFaces > 0 && currentInfo.numOpenEdges > 0) {
		const uint32_t nFaces = current.numFaces;
		std::unordered_map<GridKey, std::vector<std::pair<uint32_t, Vec>>, GridKeyHash> vtxGrid;
		uint32_t nextVid = 0;

		auto findOrAdd = [&](const Vec& p) -> uint32_t {
			auto center = getCell(p);
			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz) {
						GridKey nk = { std::get<0>(center) + dx,
										std::get<1>(center) + dy,
										std::get<2>(center) + dz };
						auto it = vtxGrid.find(nk);
						if (it == vtxGrid.end()) continue;
						for (auto& [id, pos] : it->second) {
							Vec d = p - pos;
							if (glm::dot(d, d) < cellSizeSq) {
								return id;
							}
						}
					}
			uint32_t id = nextVid++;
			vtxGrid[center].emplace_back(id, p);
			return id;
		};

		std::vector<FaceInfo> faces(nFaces);
		for (uint32_t i = 0; i < nFaces; i++) {
			Face f = current.GetFace(i);
			Vec a = current.GetPoint(f.i0);
			Vec b = current.GetPoint(f.i1);
			Vec c = current.GetPoint(f.i2);
			Vec crossP = glm::cross(b - a, c - a);
			double crossLen = glm::length(crossP);
			double e0 = glm::length(b - a);
			double e1 = glm::length(c - b);
			double e2 = glm::length(a - c);
			faces[i].vid = { findOrAdd(a), findOrAdd(b), findOrAdd(c) };
			faces[i].verts[0] = a;
			faces[i].verts[1] = b;
			faces[i].verts[2] = c;
			faces[i].pId = static_cast<uint32_t>(f.pId);
			faces[i].area = crossLen * 0.5;
			faces[i].maxEdge = std::max(e0, std::max(e1, e2));
		}

		std::unordered_map<EKey, std::vector<uint32_t>, EKeyHash> edgeFaces;
		for (uint32_t i = 0; i < nFaces; i++) {
			for (int e = 0; e < 3; e++) {
				edgeFaces[makeEKey(faces[i].vid[e], faces[i].vid[(e + 1) % 3])].push_back(i);
			}
		}

		std::vector<CandidateFace> candidates;
		for (uint32_t i = 0; i < nFaces; i++) {
			int boundaryEdges = 0;
			int crowdedEdges = 0;
			for (int e = 0; e < 3; e++) {
				const auto it = edgeFaces.find(makeEKey(faces[i].vid[e], faces[i].vid[(e + 1) % 3]));
				if (it == edgeFaces.end()) continue;
				if (it->second.size() == 1) boundaryEdges++;
				else if (it->second.size() > 2) crowdedEdges++;
			}

			const bool isBridgeFace = (boundaryEdges == 1 && crowdedEdges == 2);
			const bool isSeamFace = (boundaryEdges == 2 && crowdedEdges == 1);
			if (!isBridgeFace && !isSeamFace) continue;
			if (boundaryEdges + crowdedEdges != 3) continue;

			double minAltitude = faces[i].maxEdge > 1e-15 ? 2.0 * faces[i].area / faces[i].maxEdge : 0.0;
			if (faces[i].area > areaThreshold && minAltitude > minAltitudeThreshold) continue;

			candidates.push_back({ i, faces[i].area, minAltitude });
		}

		if (candidates.empty()) {
			break;
		}

		std::sort(candidates.begin(), candidates.end(), [](const CandidateFace& a, const CandidateFace& b) {
			if (a.area != b.area) return a.area < b.area;
			if (a.minAltitude != b.minAltitude) return a.minAltitude < b.minAltitude;
			return a.faceIndex < b.faceIndex;
			});

		bool acceptedAnyCandidate = false;
		for (const auto& candidate : candidates) {
			Geometry rebuilt;
			rebuilt.planes = current.planes;
			rebuilt.hasPlanes = current.hasPlanes;
			rebuilt.data = current.data;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (i == candidate.faceIndex) continue;
				rebuilt.AddFace(faces[i].verts[0], faces[i].verts[1], faces[i].verts[2], faces[i].pId);
			}

			if (rebuilt.numFaces == 0) continue;

			auto infoRebuilt = meshCleanup::isMeshWatertight(rebuilt);
			if (!TopologyStrictlyImprovedWithoutRegression(currentInfo, infoRebuilt)) {
				continue;
			}

			current = std::move(rebuilt);
			current.hasPlanes = false;
			currentInfo = infoRebuilt;
			removedFaces++;
			acceptedAnyCandidate = true;
			break;
		}

		if (!acceptedAnyCandidate) {
			break;
		}
	}

	if (removedFaces == 0) {
		return 0;
	}

	geom = std::move(current);
	geom.hasPlanes = false;
	meshInfoResult = currentInfo;
	return removedFaces;
}

// Patch coplanar holes: find boundary-edge loops and fill them with earcut triangulation 
// when all loop vertices are coplanar.
uint32_t meshCleanup::PatchCoplanarHoles(Geometry& geom, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult) {
	if (meshInfoInput.numOpenEdges == 0) {
		meshInfoResult = meshInfoInput;
		return 0;
	}
	const uint32_t nFaces = geom.numFaces;
	if (nFaces == 0) {
		meshInfoResult = meshInfoInput;
		return 0;
	}

	// Save backup for revert if patching doesn't improve the mesh
	Geometry backup = geom;

	// -- vertex deduplication (same spatial-hash as PostBooleanOperationMeshCleanup)
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
		double area;
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
		faces[i].area = len * 0.5;
		faces[i].pId = static_cast<uint32_t>(f.pId);
	}

	// Build edge -> face adjacency
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

	if (!hasBoundaryEdges) {
		return 0; // mesh is already closed
	}

	struct BoundaryEdge { uint32_t v0, v1; uint32_t faceIdx; };
	std::vector<BoundaryEdge> boundaryEdges;
	std::unordered_map<uint32_t, std::vector<uint32_t>> boundaryAdj;
	std::unordered_map<uint32_t, std::vector<size_t>> boundaryVertEdges;
	auto dirEdgeKey = [](uint32_t from, uint32_t to) -> uint64_t {
		return (static_cast<uint64_t>(from) << 32) | static_cast<uint64_t>(to);
	};
	auto undirectedEdgeKey = [&](uint32_t a, uint32_t b) -> uint64_t {
		EKey ek = makeEKey(a, b);
		return dirEdgeKey(ek.v0, ek.v1);
	};
	std::unordered_map<uint64_t, uint32_t> boundaryEdgeFace;
	for (auto& [ek, fl] : edgeFaces) {
		if (fl.size() == 1) {
			boundaryEdges.push_back({ ek.v0, ek.v1, fl[0] });
			size_t boundaryEdgeIdx = boundaryEdges.size() - 1;
			boundaryAdj[ek.v0].push_back(ek.v1);
			boundaryAdj[ek.v1].push_back(ek.v0);
			boundaryVertEdges[ek.v0].push_back(boundaryEdgeIdx);
			boundaryVertEdges[ek.v1].push_back(boundaryEdgeIdx);
			boundaryEdgeFace[undirectedEdgeKey(ek.v0, ek.v1)] = fl[0];
		}
	}

	const double PLANE_EPS = 5e-5;
	const double EDGE_LOOP_PLANE_EPS = std::max(PLANE_EPS * 4.0, toleranceVectorEquality * 3.0);
	const double EDGE_LOOP_PLANE_DOT = 0.985;

	std::vector<std::vector<uint32_t>> loops;
	std::vector<uint32_t> loopAdjacentFace;

	auto computePlaneFromVertices = [&](const std::vector<uint32_t>& componentVerts,
										Vec& planeNormal,
										Vec& planePoint) -> bool {
		if (componentVerts.size() < 3) return false;

		planePoint = canonPos[componentVerts[0]];
		planeNormal = Vec(0);
		for (size_t i = 1; i < componentVerts.size(); i++) {
			Vec e1 = canonPos[componentVerts[i]] - planePoint;
			for (size_t j = i + 1; j < componentVerts.size(); j++) {
				Vec e2 = canonPos[componentVerts[j]] - planePoint;
				Vec cr = glm::cross(e1, e2);
				double len = glm::length(cr);
				if (len > 1e-12) {
					planeNormal = cr / len;
					for (uint32_t vid : componentVerts) {
						double dist = std::abs(glm::dot(canonPos[vid] - planePoint, planeNormal));
						if (dist > PLANE_EPS) {
							return false;
						}
					}
					return true;
				}
			}
		}
		return false;
	};

	auto loopSignedArea2D = [](const std::vector<std::array<double, 2>>& polygon) -> double {
		if (polygon.size() < 3) return 0.0;
		double area = 0.0;
		for (size_t i = 0; i < polygon.size(); i++) {
			const auto& a = polygon[i];
			const auto& b = polygon[(i + 1) % polygon.size()];
			area += a[0] * b[1] - b[0] * a[1];
		}
		return area * 0.5;
	};

	auto canonicalizeCycle = [](const std::vector<uint32_t>& cycle) -> std::vector<uint32_t> {
		if (cycle.empty()) return {};

		std::vector<uint32_t> best;
		bool hasBest = false;
		auto consider = [&](const std::vector<uint32_t>& seq) {
			size_t n = seq.size();
			for (size_t offset = 0; offset < n; offset++) {
				std::vector<uint32_t> candidate;
				candidate.reserve(n);
				for (size_t i = 0; i < n; i++) {
					candidate.push_back(seq[(offset + i) % n]);
				}
				if (!hasBest || candidate < best) {
					best = std::move(candidate);
					hasBest = true;
				}
			}
		};

		consider(cycle);
		std::vector<uint32_t> reversed(cycle.rbegin(), cycle.rend());
		consider(reversed);
		return best;
	};

	auto facesArePlaneCompatible = [&](uint32_t faceAIdx, uint32_t faceBIdx) -> bool {
		const FaceInfo& faceA = faces[faceAIdx];
		const FaceInfo& faceB = faces[faceBIdx];
		if (std::abs(glm::dot(faceA.normal, faceB.normal)) < EDGE_LOOP_PLANE_DOT) {
			return false;
		}
		double planeDist = std::abs(glm::dot(
			canonPos[faceB.vid[0]] - canonPos[faceA.vid[0]],
			faceA.normal));
		return planeDist <= EDGE_LOOP_PLANE_EPS;
	};

	auto extractPlaneLocalLoops = [&](const std::unordered_set<uint32_t>& subSet) {
		std::unordered_set<size_t> eligibleEdges;
		for (size_t edgeIdx = 0; edgeIdx < boundaryEdges.size(); ++edgeIdx) {
			const BoundaryEdge& be = boundaryEdges[edgeIdx];
			if (!subSet.count(be.v0) || !subSet.count(be.v1)) continue;
			eligibleEdges.insert(edgeIdx);
		}
		if (eligibleEdges.empty()) {
			return;
		}

		std::unordered_set<size_t> visitedPlaneEdges;
		for (size_t seedEdgeIdx : eligibleEdges) {
			if (visitedPlaneEdges.count(seedEdgeIdx)) continue;

			std::vector<size_t> groupedEdges;
			std::queue<size_t> edgeQueue;
			edgeQueue.push(seedEdgeIdx);
			visitedPlaneEdges.insert(seedEdgeIdx);

			while (!edgeQueue.empty()) {
				size_t edgeIdx = edgeQueue.front();
				edgeQueue.pop();
				groupedEdges.push_back(edgeIdx);

				const BoundaryEdge& edge = boundaryEdges[edgeIdx];
				for (uint32_t vid : { edge.v0, edge.v1 }) {
					auto itTouching = boundaryVertEdges.find(vid);
					if (itTouching == boundaryVertEdges.end()) continue;
					for (size_t nbEdgeIdx : itTouching->second) {
						if (!eligibleEdges.count(nbEdgeIdx) || visitedPlaneEdges.count(nbEdgeIdx)) continue;
						if (!facesArePlaneCompatible(edge.faceIdx, boundaryEdges[nbEdgeIdx].faceIdx)) continue;
						visitedPlaneEdges.insert(nbEdgeIdx);
						edgeQueue.push(nbEdgeIdx);
					}
				}
			}

			if (groupedEdges.size() < 3) continue;

			std::unordered_map<uint32_t, std::vector<uint32_t>> localAdj;
			std::unordered_set<uint32_t> localVerts;
			for (size_t edgeIdx : groupedEdges) {
				const BoundaryEdge& edge = boundaryEdges[edgeIdx];
				localAdj[edge.v0].push_back(edge.v1);
				localAdj[edge.v1].push_back(edge.v0);
				localVerts.insert(edge.v0);
				localVerts.insert(edge.v1);
			}

			std::unordered_map<uint32_t, int> localDegree;
			std::queue<uint32_t> pruneQueue;
			std::unordered_set<uint32_t> prunedVerts;
			for (uint32_t vid : localVerts) {
				int deg = static_cast<int>(localAdj[vid].size());
				localDegree[vid] = deg;
				if (deg < 2) {
					pruneQueue.push(vid);
				}
			}

			while (!pruneQueue.empty()) {
				uint32_t vid = pruneQueue.front();
				pruneQueue.pop();
				if (!prunedVerts.insert(vid).second) continue;
				for (uint32_t nb : localAdj[vid]) {
					if (prunedVerts.count(nb)) continue;
					int& deg = localDegree[nb];
					if (--deg < 2) {
						pruneQueue.push(nb);
					}
				}
			}

			std::unordered_set<uint32_t> coreVisited;
			for (uint32_t startVid : localVerts) {
				if (prunedVerts.count(startVid) || coreVisited.count(startVid)) continue;

				std::vector<uint32_t> coreVerts;
				std::queue<uint32_t> coreQueue;
				coreQueue.push(startVid);
				coreVisited.insert(startVid);
				bool allDegreeTwo = true;

				while (!coreQueue.empty()) {
					uint32_t vid = coreQueue.front();
					coreQueue.pop();
					coreVerts.push_back(vid);

					int liveDegree = 0;
					for (uint32_t nb : localAdj[vid]) {
						if (prunedVerts.count(nb)) continue;
						++liveDegree;
						if (coreVisited.insert(nb).second) {
							coreQueue.push(nb);
						}
					}
					if (liveDegree != 2) {
						allDegreeTwo = false;
					}
				}

				if (!allDegreeTwo || coreVerts.size() < 3) continue;

				uint32_t loopStart = *std::min_element(coreVerts.begin(), coreVerts.end());
				std::vector<uint32_t> neighbors;
				for (uint32_t nb : localAdj[loopStart]) {
					if (prunedVerts.count(nb)) continue;
					neighbors.push_back(nb);
				}
				if (neighbors.size() != 2) continue;

				std::vector<uint32_t> orderedLoop;
				bool foundLoop = false;
				for (uint32_t firstNext : neighbors) {
					orderedLoop.clear();
					orderedLoop.push_back(loopStart);
					uint32_t prev = loopStart;
					uint32_t cur = firstNext;

					for (size_t safety = 0; safety <= coreVerts.size(); ++safety) {
						orderedLoop.push_back(cur);
						if (cur == loopStart) break;

						std::vector<uint32_t> nextCandidates;
						for (uint32_t nb : localAdj[cur]) {
							if (prunedVerts.count(nb) || nb == prev) continue;
							nextCandidates.push_back(nb);
						}
						if (nextCandidates.size() != 1) break;

						prev = cur;
						cur = nextCandidates.front();
					}

					if (orderedLoop.size() == coreVerts.size() + 1 && orderedLoop.back() == loopStart) {
						orderedLoop.pop_back();
						foundLoop = true;
						break;
					}
				}

				if (!foundLoop) continue;

				loops.push_back(std::move(orderedLoop));
				loopAdjacentFace.push_back(boundaryEdges[groupedEdges.front()].faceIdx);
			}
		}
	};

	std::unordered_set<uint32_t> visitedBoundaryVerts;
	for (const auto& [seed, _] : boundaryAdj) {
		if (visitedBoundaryVerts.count(seed)) continue;

		std::vector<uint32_t> componentVerts;
		std::queue<uint32_t> componentQueue;
		componentQueue.push(seed);
		visitedBoundaryVerts.insert(seed);

		while (!componentQueue.empty()) {
			uint32_t v = componentQueue.front();
			componentQueue.pop();
			componentVerts.push_back(v);

			for (uint32_t nb : boundaryAdj[v]) {
				if (visitedBoundaryVerts.insert(nb).second) {
					componentQueue.push(nb);
				}
			}
		}

		if (componentVerts.size() < 3) continue;

		std::unordered_set<uint32_t> componentSet(componentVerts.begin(), componentVerts.end());
		std::unordered_map<uint32_t, int> coreDegree;
		std::queue<uint32_t> pruneQueue;
		std::unordered_set<uint32_t> prunedVerts;

		for (uint32_t v : componentVerts) {
			int deg = 0;
			for (uint32_t nb : boundaryAdj[v]) {
				if (componentSet.count(nb)) deg++;
			}
			coreDegree[v] = deg;
			if (deg < 2) {
				pruneQueue.push(v);
			}
		}

		while (!pruneQueue.empty()) {
			uint32_t v = pruneQueue.front();
			pruneQueue.pop();
			if (!prunedVerts.insert(v).second) continue;

			for (uint32_t nb : boundaryAdj[v]) {
				if (!componentSet.count(nb) || prunedVerts.count(nb)) continue;
				int& deg = coreDegree[nb];
				if (--deg < 2) {
					pruneQueue.push(nb);
				}
			}
		}

		std::unordered_set<uint32_t> coreVisited;
		for (uint32_t coreSeed : componentVerts) {
			if (prunedVerts.count(coreSeed) || coreVisited.count(coreSeed)) continue;

			std::vector<uint32_t> subVerts;
			std::queue<uint32_t> subQueue;
			subQueue.push(coreSeed);
			coreVisited.insert(coreSeed);

			while (!subQueue.empty()) {
				uint32_t v = subQueue.front();
				subQueue.pop();
				subVerts.push_back(v);

				for (uint32_t nb : boundaryAdj[v]) {
					if (prunedVerts.count(nb) || !componentSet.count(nb)) continue;
					if (coreVisited.insert(nb).second) {
						subQueue.push(nb);
					}
				}
			}

			if (subVerts.size() < 3) continue;

			Vec componentPlaneNormal(0);
			Vec componentPlanePoint(0);
			bool componentIsCoplanar = computePlaneFromVertices(subVerts, componentPlaneNormal, componentPlanePoint);
			(void)componentPlanePoint;

			std::unordered_set<uint32_t> subSet(subVerts.begin(), subVerts.end());

			if (componentIsCoplanar) {
				int dropAxis = SelectProjectionAxis(componentPlaneNormal);
				std::unordered_map<uint32_t, std::array<double, 2>> projectedPos;
				projectedPos.reserve(subVerts.size());
				for (uint32_t vid : subVerts) {
					projectedPos[vid] = ProjectTo2D(canonPos[vid], dropAxis);
				}

				std::unordered_map<uint32_t, std::vector<uint32_t>> orderedNeighbors;
				orderedNeighbors.reserve(subVerts.size());
				for (uint32_t vid : subVerts) {
					std::vector<std::tuple<double, double, uint32_t>> around;
					around.reserve(boundaryAdj[vid].size());
					const auto& base = projectedPos[vid];
					for (uint32_t nb : boundaryAdj[vid]) {
						if (!subSet.count(nb)) continue;
						const auto& p = projectedPos[nb];
						double dx = p[0] - base[0];
						double dy = p[1] - base[1];
						double angle = std::atan2(dy, dx);
						double dist2 = dx * dx + dy * dy;
						around.emplace_back(angle, dist2, nb);
					}

					if (around.size() < 2) continue;

					std::sort(around.begin(), around.end(),
								[](const auto& a, const auto& b) {
									if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
									if (std::get<1>(a) != std::get<1>(b)) return std::get<1>(a) < std::get<1>(b);
									return std::get<2>(a) < std::get<2>(b);
								});

					auto& ordered = orderedNeighbors[vid];
					ordered.reserve(around.size());
					for (const auto& item : around) {
						ordered.push_back(std::get<2>(item));
					}
				}

				struct CandidateLoop {
					std::vector<uint32_t> vid;
					uint32_t adjFace;
				};
				std::vector<CandidateLoop> componentLoops;
				std::unordered_set<uint64_t> visitedDirected;

				for (const auto& be : boundaryEdges) {
					if (!subSet.count(be.v0) || !subSet.count(be.v1)) continue;

					for (int dir = 0; dir < 2; dir++) {
						uint32_t startFrom = dir == 0 ? be.v0 : be.v1;
						uint32_t startTo = dir == 0 ? be.v1 : be.v0;
						uint64_t startKey = dirEdgeKey(startFrom, startTo);
						if (visitedDirected.count(startKey)) continue;

						std::vector<uint32_t> loop;
						std::unordered_set<uint32_t> seenVerts;
						loop.reserve(subVerts.size() + 1);
						loop.push_back(startFrom);
						seenVerts.insert(startFrom);

						uint32_t prev = startFrom;
						uint32_t cur = startTo;
						bool valid = true;

						for (size_t safety = 0; safety < boundaryEdges.size() + 1; safety++) {
							if (!visitedDirected.insert(dirEdgeKey(prev, cur)).second) {
								valid = false;
								break;
							}

							if (cur == startFrom) {
								loop.push_back(cur);
							}
							else {
								if (!seenVerts.insert(cur).second) {
									valid = false;
									break;
								}
								loop.push_back(cur);
							}

							auto adjIt = orderedNeighbors.find(cur);
							if (adjIt == orderedNeighbors.end() || adjIt->second.size() < 2) {
								valid = false;
								break;
							}

							const auto& around = adjIt->second;
							auto revIt = std::find(around.begin(), around.end(), prev);
							if (revIt == around.end()) {
								valid = false;
								break;
							}

							size_t idx = static_cast<size_t>(std::distance(around.begin(), revIt));
							uint32_t next = around[(idx + around.size() - 1) % around.size()];

							prev = cur;
							cur = next;

							if (prev == startFrom && cur == startTo) {
								break;
							}
						}

						if (!valid || loop.size() < 4 || loop.back() != loop.front()) {
							continue;
						}

						loop.pop_back();
						if (loop.size() < 3) continue;

						std::vector<std::array<double, 2>> projectedLoop;
						projectedLoop.reserve(loop.size());
						for (uint32_t vid : loop) {
							projectedLoop.push_back(projectedPos[vid]);
						}

						double signedArea = loopSignedArea2D(projectedLoop);
						if (signedArea <= 1e-12) continue;

						componentLoops.push_back({ std::move(loop), be.faceIdx });
					}
				}

				for (auto& candidate : componentLoops) {
					loops.push_back(std::move(candidate.vid));
					loopAdjacentFace.push_back(candidate.adjFace);
				}
			}
			else {
				const size_t MAX_ENUM_VERTS = 48;
				const size_t MAX_ENUM_EDGES = 64;
				const size_t MAX_ENUM_CYCLES = 256;

				size_t subEdgeCount = 0;
				std::unordered_map<uint32_t, std::vector<uint32_t>> subAdj;
				subAdj.reserve(subVerts.size());
				for (uint32_t vid : subVerts) {
					auto& out = subAdj[vid];
					for (uint32_t nb : boundaryAdj[vid]) {
						if (!subSet.count(nb)) continue;
						out.push_back(nb);
					}
					subEdgeCount += out.size();
				}
				subEdgeCount /= 2;

				if (subVerts.size() > MAX_ENUM_VERTS || subEdgeCount > MAX_ENUM_EDGES) {
					extractPlaneLocalLoops(subSet);
					continue;
				}

				std::vector<uint32_t> sortedStarts = subVerts;
				std::sort(sortedStarts.begin(), sortedStarts.end());
				std::set<std::vector<uint32_t>> uniqueCycles;
				std::vector<uint32_t> path;
				std::unordered_set<uint32_t> pathVisited;

				std::function<void(uint32_t, uint32_t)> enumerateFrom = [&](uint32_t start, uint32_t current) {
					if (uniqueCycles.size() >= MAX_ENUM_CYCLES) return;

					for (uint32_t nb : subAdj[current]) {
						if (nb == start) {
							if (path.size() >= 3) {
								std::vector<uint32_t> canonical = canonicalizeCycle(path);
								if (uniqueCycles.insert(canonical).second) {
									Vec loopPlaneNormal(0);
									Vec loopPlanePoint(0);
									if (!computePlaneFromVertices(canonical, loopPlaneNormal, loopPlanePoint)) {
										continue;
									}

									int loopAxis = SelectProjectionAxis(loopPlaneNormal);
									std::vector<std::array<double, 2>> projectedLoop;
									projectedLoop.reserve(canonical.size());
									for (uint32_t vid : canonical) {
										projectedLoop.push_back(ProjectTo2D(canonPos[vid], loopAxis));
									}
									double area = std::abs(loopSignedArea2D(projectedLoop));
									if (area <= 1e-12) {
										continue;
									}

									uint32_t adjFace = 0;
									if (canonical.size() >= 2) {
										for (size_t edgeIdx = 0; edgeIdx < canonical.size(); edgeIdx++) {
											uint32_t a = canonical[edgeIdx];
											uint32_t b = canonical[(edgeIdx + 1) % canonical.size()];
											auto faceIt = boundaryEdgeFace.find(undirectedEdgeKey(a, b));
											if (faceIt != boundaryEdgeFace.end()) {
												adjFace = faceIt->second;
												break;
											}
										}
									}

									loops.push_back(std::move(canonical));
									loopAdjacentFace.push_back(adjFace);
								}
							}
							continue;
						}

						if (pathVisited.count(nb)) continue;
						if (nb < start) continue;
						if (path.size() >= subVerts.size()) continue;

						path.push_back(nb);
						pathVisited.insert(nb);
						enumerateFrom(start, nb);
						pathVisited.erase(nb);
						path.pop_back();
					}
				};

				for (uint32_t start : sortedStarts) {
					path.clear();
					pathVisited.clear();
					path.push_back(start);
					pathVisited.insert(start);
					enumerateFrom(start, start);
				}
			}
		}
	}

	if (loops.empty()) {
		meshInfoResult = meshInfoInput;
		return 0;
	}
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
		return 0;
	}

	// -- Global parent-child containment detection across ALL loops --
	// For each candidate parent, only compare loops that lie on the same
	// geometric plane. Overlapping projections on different planes can create
	// bogus polygons and membrane-like caps.
	struct ProjectedLoop {
		size_t loopIndex;
		int dropAxis;  // per-loop projection axis
		std::vector<std::array<double, 2>> projected;
		double area = 0.0;
		std::array<double, 2> bboxMin{ 0.0, 0.0 };
		std::array<double, 2> bboxMax{ 0.0, 0.0 };
		int planeGroup = -1;
		int parent = -1;
		int depth = 0;
	};

	const double LOOP_PLANE_GROUP_EPS = std::max(PLANE_EPS * 4.0, toleranceVectorEquality * 3.0);
	const double LOOP_PLANE_GROUP_DOT = 0.985;

	auto loopsArePlaneCompatible = [&](const LoopInfo& a, const LoopInfo& b) -> bool {
		if (SelectProjectionAxis(a.planeNormal) != SelectProjectionAxis(b.planeNormal)) return false;

		double normalDot = std::abs(glm::dot(a.planeNormal, b.planeNormal));
		if (normalDot < LOOP_PLANE_GROUP_DOT) return false;

		Vec delta = b.planePoint - a.planePoint;
		double separationA = std::abs(glm::dot(delta, a.planeNormal));
		double separationB = std::abs(glm::dot(delta, b.planeNormal));
		return std::max(separationA, separationB) <= LOOP_PLANE_GROUP_EPS;
	};

	auto loopLiesOnPlane = [&](const LoopInfo& loop, const Vec& planeNormal, const Vec& planePoint,
								double eps, bool useRawBoundaryPositions) -> bool {
		for (uint32_t vid : loop.vid) {
			Vec p = useRawBoundaryPositions ? getRawBoundaryPos(vid) : canonPos[vid];
			if (std::abs(glm::dot(p - planePoint, planeNormal)) > eps) {
				return false;
			}
		}
		return true;
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
			pl.bboxMin = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
			pl.bboxMax = { -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };
			for (const auto& p : pl.projected) {
				pl.bboxMin[0] = std::min(pl.bboxMin[0], p[0]);
				pl.bboxMin[1] = std::min(pl.bboxMin[1], p[1]);
				pl.bboxMax[0] = std::max(pl.bboxMax[0], p[0]);
				pl.bboxMax[1] = std::max(pl.bboxMax[1], p[1]);
			}
			projectedLoops.push_back(std::move(pl));
		}
	}

	int nextPlaneGroup = 0;
	for (size_t i = 0; i < projectedLoops.size(); i++) {
		int planeGroup = -1;
		const LoopInfo& loop = validLoops[projectedLoops[i].loopIndex];
		for (size_t j = 0; j < i; j++) {
			const LoopInfo& otherLoop = validLoops[projectedLoops[j].loopIndex];
			if (!loopsArePlaneCompatible(loop, otherLoop)) continue;
			planeGroup = projectedLoops[j].planeGroup;
			break;
		}
		if (planeGroup < 0) {
			planeGroup = nextPlaneGroup++;
		}
		projectedLoops[i].planeGroup = planeGroup;
	}

	// Find parent (smallest enclosing loop) for each loop.
	// For each candidate parent, re-project the child's sample point using
	// the PARENT's axis and test containment against the parent's polygon.
	// Restrict the search to loops in the same plane group.
	for (size_t loopIdx = 0; loopIdx < projectedLoops.size(); loopIdx++) {
		double bestParentArea = std::numeric_limits<double>::max();
		const LoopInfo& loop = validLoops[projectedLoops[loopIdx].loopIndex];

		for (size_t candidateIdx = 0; candidateIdx < projectedLoops.size(); candidateIdx++) {
			if (candidateIdx == loopIdx) continue;
			if (projectedLoops[candidateIdx].planeGroup != projectedLoops[loopIdx].planeGroup) continue;
			if (projectedLoops[candidateIdx].area <= projectedLoops[loopIdx].area + 1e-12) continue;
			const LoopInfo& candidateLoop = validLoops[projectedLoops[candidateIdx].loopIndex];
			if (!loopsArePlaneCompatible(candidateLoop, loop)) continue;

			// Project the child's first vertex using the CANDIDATE PARENT's axis
			int parentAxis = projectedLoops[candidateIdx].dropAxis;
			uint32_t sampleVid = loop.vid[0];
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

	// Skip paired opening mouths: two similar coplanar loops on opposite
	// host planes indicate an intentional through-opening corridor, not a
	// missing cap that should be patched.
	std::array<double, 3> dominantNormalArea = { 0.0, 0.0, 0.0 };
	for (const auto& face : faces) {
		dominantNormalArea[0] += face.area * std::abs(face.normal.x);
		dominantNormalArea[1] += face.area * std::abs(face.normal.y);
		dominantNormalArea[2] += face.area * std::abs(face.normal.z);
	}
	int dominantAxis = 0;
	if (dominantNormalArea[1] > dominantNormalArea[dominantAxis]) dominantAxis = 1;
	if (dominantNormalArea[2] > dominantNormalArea[dominantAxis]) dominantAxis = 2;

	auto axisComponent = [&](const Vec& v) -> double {
		return dominantAxis == 0 ? v.x : (dominantAxis == 1 ? v.y : v.z);
	};

	std::vector<bool> skipProjectedLoop(projectedLoops.size(), false);
	constexpr bool kSkipPairedOpeningMouths = false;
	if (kSkipPairedOpeningMouths) {
		for (size_t i = 0; i < projectedLoops.size(); i++) {
			const LoopInfo& loopA = validLoops[projectedLoops[i].loopIndex];
			if (std::abs(axisComponent(loopA.planeNormal)) < 0.85) continue;
			if (projectedLoops[i].dropAxis != dominantAxis) continue;

			for (size_t j = i + 1; j < projectedLoops.size(); j++) {
				const LoopInfo& loopB = validLoops[projectedLoops[j].loopIndex];
				if (std::abs(axisComponent(loopB.planeNormal)) < 0.85) continue;
				if (projectedLoops[j].dropAxis != dominantAxis) continue;
				if (std::abs(glm::dot(loopA.planeNormal, loopB.planeNormal)) < 0.95) continue;

				double minX = std::max(projectedLoops[i].bboxMin[0], projectedLoops[j].bboxMin[0]);
				double minY = std::max(projectedLoops[i].bboxMin[1], projectedLoops[j].bboxMin[1]);
				double maxX = std::min(projectedLoops[i].bboxMax[0], projectedLoops[j].bboxMax[0]);
				double maxY = std::min(projectedLoops[i].bboxMax[1], projectedLoops[j].bboxMax[1]);
				double overlapWidth = std::max(0.0, maxX - minX);
				double overlapHeight = std::max(0.0, maxY - minY);
				double overlapArea = overlapWidth * overlapHeight;
				double minLoopArea = std::min(projectedLoops[i].area, projectedLoops[j].area);
				if (minLoopArea <= 1e-12 || overlapArea < minLoopArea * 0.8) continue;

				double separation = std::abs(glm::dot(loopB.planePoint - loopA.planePoint, loopA.planeNormal));
				double transverseSpan = std::max({ projectedLoops[i].bboxMax[0] - projectedLoops[i].bboxMin[0],
													projectedLoops[i].bboxMax[1] - projectedLoops[i].bboxMin[1],
													projectedLoops[j].bboxMax[0] - projectedLoops[j].bboxMin[0],
													projectedLoops[j].bboxMax[1] - projectedLoops[j].bboxMin[1] });
				if (separation <= toleranceVectorEquality * 2.0) continue;
				if (transverseSpan <= 1e-12) continue;
				if (separation >= transverseSpan * 0.25) continue;

				skipProjectedLoop[i] = true;
				skipProjectedLoop[j] = true;
			}
		}
	}

	// Triangulate: even-depth loops are outer boundaries, odd-depth are holes
	for (size_t outerIdx = 0; outerIdx < projectedLoops.size(); outerIdx++) {
		if ((projectedLoops[outerIdx].depth & 1) != 0) continue;
		if (skipProjectedLoop[outerIdx]) continue;

		const LoopInfo& outerLoop = validLoops[projectedLoops[outerIdx].loopIndex];
		Vec outerPlanePoint = getRawBoundaryPos(outerLoop.vid[0]);
		if (!loopLiesOnPlane(outerLoop, outerLoop.planeNormal, outerPlanePoint, LOOP_PLANE_GROUP_EPS, true)) {
			continue;
		}

		std::vector<std::vector<std::array<double, 2>>> polygon;
		std::vector<uint32_t> polygonVertexIds;
		polygon.reserve(projectedLoops.size());
		polygonVertexIds.reserve(outerLoop.vid.size());

		polygon.push_back(projectedLoops[outerIdx].projected);
		polygonVertexIds.insert(polygonVertexIds.end(), outerLoop.vid.begin(), outerLoop.vid.end());

		for (size_t holeIdx = 0; holeIdx < projectedLoops.size(); holeIdx++) {
			if (projectedLoops[holeIdx].parent != static_cast<int>(outerIdx)) continue;
			if (skipProjectedLoop[holeIdx]) continue;

			const auto& holeLoop = validLoops[projectedLoops[holeIdx].loopIndex];
			if (projectedLoops[holeIdx].planeGroup != projectedLoops[outerIdx].planeGroup) continue;
			if (!loopsArePlaneCompatible(outerLoop, holeLoop)) continue;
			if (!loopLiesOnPlane(holeLoop, outerLoop.planeNormal, outerPlanePoint, LOOP_PLANE_GROUP_EPS, true)) {
				continue;
			}

			// Re-project hole using the OUTER loop's axis for consistent earcut input
			int outerAxis = projectedLoops[outerIdx].dropAxis;
			std::vector<std::array<double, 2>> holeProjected;
			holeProjected.reserve(holeLoop.vid.size());
			for (uint32_t vid : holeLoop.vid) {
				holeProjected.push_back(ProjectTo2D(canonPos[vid], outerAxis));
			}
			polygon.push_back(std::move(holeProjected));
			polygonVertexIds.insert(polygonVertexIds.end(), holeLoop.vid.begin(), holeLoop.vid.end());
		}

		bool polygonIsCoplanar = true;
		for (uint32_t vid : polygonVertexIds) {
			Vec p = getRawBoundaryPos(vid);
			if (std::abs(glm::dot(p - outerPlanePoint, outerLoop.planeNormal)) > LOOP_PLANE_GROUP_EPS) {
				polygonIsCoplanar = false;
				break;
			}
		}
		if (!polygonIsCoplanar) continue;

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
		const uint32_t addedPatchFaces = geom.numFaces - nFaces;

		if (MeshPenaltyImprovedNoRegression(meshInfoInput, infoPatched)) {
			meshInfoResult = infoPatched;
			return addedPatchFaces;
		}
		else {
			// patching didn't improve -- revert
			geom = std::move(backup);
			meshInfoResult = meshInfoInput;
			return 0;
		}
	}
	else {
		meshInfoResult = meshInfoInput;
		return 0;
	}
}

uint32_t meshCleanup::RemoveThinMembranes(Geometry& workingMesh, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult, double thinThresholdOverride) {
	meshInfoResult = meshInfoInput;

	const Geometry baseMesh = workingMesh;
	const uint32_t nFaces = baseMesh.numFaces;
	if (nFaces < 4) return 0;

	struct FV {
		Vec a, b, c;
		uint32_t pId;
		Vec center, normal;
		double area, maxEdge;
	};

	std::vector<FV> fv(nFaces);
	double totalAreaBefore = 0.0;
	for (uint32_t i = 0; i < nFaces; i++) {
		Face f = baseMesh.GetFace(i);
		Vec a = baseMesh.GetPoint(f.i0);
		Vec b = baseMesh.GetPoint(f.i1);
		Vec c = baseMesh.GetPoint(f.i2);
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
		totalAreaBefore += fv[i].area;
	}

	const double cellSize = toleranceVectorEquality * 1.5;
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
					if (it == vtxGrid.end()) continue;
					for (auto& [id, pos] : it->second) {
						Vec d = p - pos;
						if (glm::dot(d, d) < cellSizeSq) return id;
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
		for (int e = 0; e < 3; e++) {
			edgeFaces[makeEKey(fvid[i][e], fvid[i][(e + 1) % 3])].push_back(i);
		}
	}

	std::vector<std::vector<uint32_t>> faceAdj(nFaces);
	for (auto& [ek, fl] : edgeFaces) {
		for (size_t a = 0; a < fl.size(); a++) {
			for (size_t b = a + 1; b < fl.size(); b++) {
				faceAdj[fl[a]].push_back(fl[b]);
				faceAdj[fl[b]].push_back(fl[a]);
			}
		}
	}

	const double THIN_THRESHOLD = thinThresholdOverride > 0 ? thinThresholdOverride : 5e-3;
	constexpr double VOLUME_THRESHOLD = 1e-6;
	std::vector<bool> thinMarked(nFaces, false);

	BVH resultBVH = MakeBVH(baseMesh);

	for (uint32_t i = 0; i < nFaces; i++) {
		if (thinMarked[i]) continue;
		if (glm::length(fv[i].normal) < 0.5) continue;

		for (int sign = -1; sign <= 1; sign += 2) {
			if (thinMarked[i]) break;

			Vec rayDir = fv[i].normal * static_cast<double>(sign);
			Vec rayEnd = fv[i].center + rayDir * THIN_THRESHOLD;

			resultBVH.IntersectRay(fv[i].center, rayDir, [&](uint32_t faceIdx) -> bool {
				if (faceIdx == i) return false;
				if (glm::dot(fv[i].normal, fv[faceIdx].normal) > -0.7) return false;

				Vec hitPos;
				double t, d_plane;
				if (!intersect_ray_triangle(fv[i].center, rayEnd, fv[faceIdx].a, fv[faceIdx].b,
					fv[faceIdx].c, hitPos, t, d_plane, false)) {
					return false;
				}

				double dist = glm::length(hitPos - fv[i].center);
				if (dist > 1e-6 && dist < THIN_THRESHOLD) {
					thinMarked[i] = true;
					thinMarked[faceIdx] = true;
					return true;
				}
				return false;
			});
		}
	}

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
			uint32_t cur = q.front();
			q.pop();
			for (uint32_t nb : manifoldAdj[cur]) {
				if (mCompId[nb] < 0) {
					mCompId[nb] = cid;
					q.push(nb);
				}
			}
		}
	}

	if (mNumComp > 1) {
		struct MCInfo {
			uint32_t count = 0;
			Vec centroid{ 0 };
		};
		std::vector<MCInfo> mci(mNumComp);
		for (uint32_t i = 0; i < nFaces; i++) {
			mci[mCompId[i]].count++;
			mci[mCompId[i]].centroid += fv[i].center;
		}
		for (auto& info : mci) {
			if (info.count > 0) {
				info.centroid /= static_cast<double>(info.count);
			}
		}

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
				thinMarked[i] = true;
			}
		}
	}

	for (int iter = 0; iter < 20; iter++) {
		bool changed = false;
		for (uint32_t i = 0; i < nFaces; i++) {
			if (thinMarked[i]) continue;
			double minH = fv[i].maxEdge > 0 ? 2.0 * fv[i].area / fv[i].maxEdge : 0.0;
			if (minH >= THIN_THRESHOLD) continue;

			int thinNb = 0;
			int totalNb = 0;
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

	uint32_t thinCount = 0;
	double removedArea = 0.0;
	for (uint32_t i = 0; i < nFaces; i++) {
		if (!thinMarked[i]) continue;
		thinCount++;
		removedArea += fv[i].area;
	}

	// Compute total signed volume of the mesh. A near-zero volume indicates
	// the mesh is a degenerate membrane (e.g. boolean subtract fully consumed
	// the first operand). In that case, relax the safety guards below so the
	// membrane can be removed entirely.
	double totalVolume = 0.0;
	for (uint32_t i = 0; i < nFaces; i++)
		totalVolume += glm::dot(fv[i].a, glm::cross(fv[i].b, fv[i].c)) / 6.0;
	const bool isZeroVolume = std::abs(totalVolume) < VOLUME_THRESHOLD;

	if (thinCount == 0) {
		return 0;
	}
	if (thinCount >= nFaces * 3 / 4 && !isZeroVolume) {
		return 0;
	}

	const double areaLossRatio = totalAreaBefore > 1e-9 ? removedArea / totalAreaBefore : 0.0;
	if (areaLossRatio > 0.20 && !isZeroVolume) {
		return 0;
	}

	Geometry cleaned;
	cleaned.planes = baseMesh.planes;
	cleaned.hasPlanes = baseMesh.hasPlanes;
	cleaned.data = baseMesh.data;
	for (uint32_t i = 0; i < nFaces; i++) {
		if (thinMarked[i]) continue;
		cleaned.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
	}

	if (cleaned.numFaces >= nFaces) {
		return 0;
	}
	// An empty result is valid for zero-volume membranes (fully consumed operand).
	if (cleaned.numFaces == 0 && !isZeroVolume) {
		return 0;
	}

	// Zero-volume membrane fully removed: skip topology checks (empty mesh is correct).
	if (cleaned.numFaces == 0 && isZeroVolume) {
		workingMesh = std::move(cleaned);
		workingMesh.hasPlanes = false;
		meshInfoResult = meshCleanup::isMeshWatertight(workingMesh);
		return thinCount;
	}

	auto cleanedInfo = meshCleanup::isMeshWatertight(cleaned);
	if (cleanedInfo.numOpenEdges > meshInfoInput.numOpenEdges + 8) {
		return 0;
	}

	const bool improvesTopology = TopologyStrictlyImprovedWithoutRegression(meshInfoInput, cleanedInfo);
	const bool improvesPenalty = MeshPenaltyImprovedNoRegression(meshInfoInput, cleanedInfo);
	const bool reducesNonManifold =
		cleanedInfo.numNonManifoldEdges < meshInfoInput.numNonManifoldEdges;
	if (!(improvesTopology || improvesPenalty || reducesNonManifold)) {
		return 0;
	}

	workingMesh = std::move(cleaned);
	workingMesh.hasPlanes = false;
	meshInfoResult = cleanedInfo;

#ifdef DUMP_CSG_MESHES
	if (thinCount > 7 && cleanedInfo.numFaces > 87) {
		webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
		DumpDebugGeometry(inputWebIfc, "meshCleanup" + step + "-removedMembranes.obj");
	}
#endif
	return thinCount;
}

uint32_t meshCleanup::RemoveOpposedEdgeMembranes(Geometry& workingMesh, std::string step,
	const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult) {
	meshInfoResult = meshInfoInput;
	if (workingMesh.numFaces == 0 || meshInfoInput.numNonManifoldEdges == 0) {
		return 0;
	}

	const double cellSize = toleranceVectorEquality;
	const double cellSizeSq = cellSize * cellSize;
	const double maxMembraneThickness = 1e-2;
	const double maxMembraneThicknessSq = maxMembraneThickness * maxMembraneThickness;
	const double maxAltitude = 0.1;
	const double oppositeNormalDot = -0.8;
	const double coplanarNormalDot = 0.999;
	const double coplanarPlaneTol = cellSize * 4.0;

	struct FaceInfo {
		std::array<uint32_t, 3> vid;
		Vec verts[3];
		Vec normal{ 0.0 };
		uint32_t pId = 0;
		double area = 0.0;
		double maxEdge = 0.0;
		double minAltitude = 0.0;
	};

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

	using GridKey = std::tuple<int64_t, int64_t, int64_t>;
	struct GridKeyHash {
		size_t operator()(const GridKey& k) const {
			size_t h = std::hash<int64_t>()(std::get<0>(k));
			h ^= std::hash<int64_t>()(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
			h ^= std::hash<int64_t>()(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
			return h;
		}
	};

	struct CandidatePair {
		uint32_t faceA = UINT32_MAX;
		uint32_t faceB = UINT32_MAX;
		double thicknessSq = 0.0;
		double pairArea = 0.0;
	};
	struct CandidateFaceRemoval {
		uint32_t face = UINT32_MAX;
		double area = 0.0;
		double minAltitude = 0.0;
		uint32_t coplanarCrowdedHits = 0;
	};

	Geometry current = workingMesh;
	MeshInfo currentInfo = meshInfoInput;
	uint32_t removedFacesTotal = 0;

	auto getCell = [&](const Vec& p) -> GridKey {
		return { static_cast<int64_t>(std::floor(p.x / cellSize)),
				static_cast<int64_t>(std::floor(p.y / cellSize)),
				static_cast<int64_t>(std::floor(p.z / cellSize)) };
	};

	while (current.numFaces > 0 && currentInfo.numNonManifoldEdges > 0) {
		const uint32_t nFaces = current.numFaces;
		std::unordered_map<GridKey, std::vector<std::pair<uint32_t, Vec>>, GridKeyHash> vtxGrid;
		std::vector<Vec> canonicalPos;
		uint32_t nextVid = 0;

		auto findOrAdd = [&](const Vec& p) -> uint32_t {
			auto center = getCell(p);
			for (int dx = -1; dx <= 1; ++dx)
				for (int dy = -1; dy <= 1; ++dy)
					for (int dz = -1; dz <= 1; ++dz) {
						GridKey nk = { std::get<0>(center) + dx,
										std::get<1>(center) + dy,
										std::get<2>(center) + dz };
						auto it = vtxGrid.find(nk);
						if (it == vtxGrid.end()) continue;
						for (auto& [id, pos] : it->second) {
							Vec d = p - pos;
							if (glm::dot(d, d) < cellSizeSq) {
								return id;
							}
						}
					}
			uint32_t id = nextVid++;
			vtxGrid[center].emplace_back(id, p);
			canonicalPos.push_back(p);
			return id;
		};

		std::vector<FaceInfo> faces(nFaces);
		for (uint32_t i = 0; i < nFaces; ++i) {
			Face f = current.GetFace(i);
			Vec a = current.GetPoint(f.i0);
			Vec b = current.GetPoint(f.i1);
			Vec c = current.GetPoint(f.i2);
			Vec crossP = glm::cross(b - a, c - a);
			double crossLen = glm::length(crossP);
			double e0 = glm::length(b - a);
			double e1 = glm::length(c - b);
			double e2 = glm::length(a - c);
			double maxEdge = std::max(e0, std::max(e1, e2));
			faces[i].vid = { findOrAdd(a), findOrAdd(b), findOrAdd(c) };
			faces[i].verts[0] = a;
			faces[i].verts[1] = b;
			faces[i].verts[2] = c;
			faces[i].normal = crossLen > 1e-15 ? crossP / crossLen : Vec(0.0);
			faces[i].pId = static_cast<uint32_t>(f.pId);
			faces[i].area = crossLen * 0.5;
			faces[i].maxEdge = maxEdge;
			faces[i].minAltitude = maxEdge > 1e-15 ? (crossLen / maxEdge) : 0.0;
		}

		std::unordered_map<EKey, std::vector<uint32_t>, EKeyHash> edgeFaces;
		for (uint32_t i = 0; i < nFaces; ++i) {
			for (int e = 0; e < 3; ++e) {
				edgeFaces[makeEKey(faces[i].vid[e], faces[i].vid[(e + 1) % 3])].push_back(i);
			}
		}

		std::vector<uint32_t> faceCoplanarCrowdedHits(nFaces, 0);
		for (const auto& [edge, owners] : edgeFaces) {
			if (owners.size() < 3 || owners.size() > 4) continue;
			for (uint32_t faceIdx : owners) {
				const FaceInfo& face = faces[faceIdx];
				bool hasCoplanarPartner = false;
				for (uint32_t otherIdx : owners) {
					if (otherIdx == faceIdx) continue;
					const FaceInfo& other = faces[otherIdx];
					if (std::abs(glm::dot(face.normal, other.normal)) < coplanarNormalDot) continue;
					double planeDist = std::abs(glm::dot(other.verts[0] - face.verts[0], face.normal));
					if (planeDist > coplanarPlaneTol) continue;
					hasCoplanarPartner = true;
					break;
				}
				if (hasCoplanarPartner) {
					faceCoplanarCrowdedHits[faceIdx]++;
				}
			}
		}

		auto tryClosedNeighborhoodRepair = [&]() -> bool {
			if (currentInfo.numOpenEdges != 0 || currentInfo.numNonManifoldEdges == 0) {
				return false;
			}

			struct PlaneGroup {
				Vec normal{ 0.0 };
				Vec point{ 0.0 };
				uint32_t pId = 0;
				double representativeArea = 0.0;
				std::vector<uint32_t> faceIndices;
				std::unordered_set<uint32_t> vertexSet;
			};
			struct TriKey {
				uint32_t a = 0, b = 0, c = 0;
				bool operator==(const TriKey& o) const { return a == o.a && b == o.b && c == o.c; }
			};
			struct TriKeyHash {
				size_t operator()(const TriKey& k) const {
					size_t h = std::hash<uint32_t>()(k.a);
					h ^= std::hash<uint32_t>()(k.b) + 0x9e3779b9 + (h << 6) + (h >> 2);
					h ^= std::hash<uint32_t>()(k.c) + 0x9e3779b9 + (h << 6) + (h >> 2);
					return h;
				}
			};
			struct CandidateTri {
				std::array<uint32_t, 3> vid{ 0, 0, 0 };
				std::array<Vec, 3> verts{ Vec(0.0), Vec(0.0), Vec(0.0) };
				std::array<uint64_t, 3> edges{ 0, 0, 0 };
				Vec normal{ 0.0 };
				uint32_t pId = 0;
			};

			auto edgeKey64 = [](uint32_t a, uint32_t b) -> uint64_t {
				if (a > b) std::swap(a, b);
				return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
			};

			for (const auto& [crowdedEdge, owners] : edgeFaces) {
				if (owners.size() != 4) continue;

				std::unordered_set<uint32_t> neighborhoodFaces(owners.begin(), owners.end());
				std::unordered_set<uint32_t> neighborhoodVerts;
				for (uint32_t ownerIdx : owners) {
					for (uint32_t vid : faces[ownerIdx].vid) {
						neighborhoodVerts.insert(vid);
					}
				}

				for (int pass = 0; pass < 2; ++pass) {
					bool addedAny = false;
					for (uint32_t fi = 0; fi < nFaces; ++fi) {
						if (neighborhoodFaces.count(fi)) continue;
						int sharedVerts = 0;
						for (uint32_t vid : faces[fi].vid) {
							if (neighborhoodVerts.count(vid)) {
								sharedVerts++;
							}
						}
						if (sharedVerts < 2) continue;
						neighborhoodFaces.insert(fi);
						for (uint32_t vid : faces[fi].vid) {
							neighborhoodVerts.insert(vid);
						}
						addedAny = true;
						if (neighborhoodFaces.size() > 24) {
							break;
						}
					}
					if (!addedAny || neighborhoodFaces.size() > 24) {
						break;
					}
				}

				if (neighborhoodFaces.size() < 8 || neighborhoodFaces.size() > 24) {
					continue;
				}

				std::vector<PlaneGroup> planeGroups;
				for (uint32_t fi : neighborhoodFaces) {
					const FaceInfo& face = faces[fi];
					if (face.area < 1e-12 || glm::length(face.normal) < 1e-12) continue;

					int matchedGroup = -1;
					for (size_t gi = 0; gi < planeGroups.size(); ++gi) {
						double ndot = std::abs(glm::dot(face.normal, planeGroups[gi].normal));
						if (ndot < coplanarNormalDot) continue;
						double planeDist = std::abs(glm::dot(face.verts[0] - planeGroups[gi].point, planeGroups[gi].normal));
						if (planeDist > coplanarPlaneTol) continue;
						matchedGroup = static_cast<int>(gi);
						break;
					}

					if (matchedGroup < 0) {
						PlaneGroup group;
						group.normal = face.normal;
						group.point = face.verts[0];
						group.pId = face.pId;
						group.representativeArea = face.area;
						group.faceIndices.push_back(fi);
						for (uint32_t vid : face.vid) {
							group.vertexSet.insert(vid);
						}
						planeGroups.push_back(std::move(group));
					}
					else {
						PlaneGroup& group = planeGroups[matchedGroup];
						group.faceIndices.push_back(fi);
						for (uint32_t vid : face.vid) {
							group.vertexSet.insert(vid);
						}
						if (face.area > group.representativeArea) {
							group.point = face.verts[0];
							group.normal = face.normal;
							group.pId = face.pId;
							group.representativeArea = face.area;
						}
					}
				}

				if (planeGroups.empty()) {
					continue;
				}

				std::unordered_map<TriKey, CandidateTri, TriKeyHash> candidateMap;
				auto addCandidate = [&](std::array<uint32_t, 3> vid, const Vec& planeNormal, uint32_t pId) {
					std::array<uint32_t, 3> sortedVid = vid;
					std::sort(sortedVid.begin(), sortedVid.end());
					if (sortedVid[0] == sortedVid[1] || sortedVid[1] == sortedVid[2]) {
						return;
					}

					const Vec& a = canonicalPos[sortedVid[0]];
					const Vec& b = canonicalPos[sortedVid[1]];
					const Vec& c = canonicalPos[sortedVid[2]];
					double triArea = fuzzybools::areaOfTriangle(a, b, c);
					if (triArea <= 1e-10) {
						return;
					}

					TriKey key{ sortedVid[0], sortedVid[1], sortedVid[2] };
					if (candidateMap.count(key)) {
						return;
					}

					CandidateTri candidate;
					candidate.vid = sortedVid;
					candidate.verts = { a, b, c };
					candidate.edges = {
						edgeKey64(sortedVid[0], sortedVid[1]),
						edgeKey64(sortedVid[1], sortedVid[2]),
						edgeKey64(sortedVid[0], sortedVid[2])
					};
					candidate.normal = planeNormal;
					candidate.pId = pId;
					candidateMap.emplace(key, std::move(candidate));
				};

				for (uint32_t fi : neighborhoodFaces) {
					addCandidate(faces[fi].vid, faces[fi].normal, faces[fi].pId);
				}

				for (const PlaneGroup& group : planeGroups) {
					std::vector<uint32_t> groupVerts(group.vertexSet.begin(), group.vertexSet.end());
					if (groupVerts.size() > 12) {
						continue;
					}

					for (size_t ia = 0; ia < groupVerts.size(); ++ia) {
						for (size_t ib = ia + 1; ib < groupVerts.size(); ++ib) {
							for (size_t ic = ib + 1; ic < groupVerts.size(); ++ic) {
								uint32_t va = groupVerts[ia];
								uint32_t vb = groupVerts[ib];
								uint32_t vc = groupVerts[ic];
								if (std::abs(glm::dot(canonicalPos[va] - group.point, group.normal)) > coplanarPlaneTol) continue;
								if (std::abs(glm::dot(canonicalPos[vb] - group.point, group.normal)) > coplanarPlaneTol) continue;
								if (std::abs(glm::dot(canonicalPos[vc] - group.point, group.normal)) > coplanarPlaneTol) continue;
								addCandidate({ va, vb, vc }, group.normal, group.pId);
							}
						}
					}
				}

				if (candidateMap.empty() || candidateMap.size() > 220) {
					continue;
				}

				std::vector<CandidateTri> candidates;
				candidates.reserve(candidateMap.size());
				for (auto& [key, candidate] : candidateMap) {
					candidates.push_back(candidate);
				}

				std::unordered_map<uint64_t, std::vector<uint32_t>> edgeToCandidates;
				std::unordered_set<uint64_t> relevantEdgeSet;
				for (const CandidateTri& candidate : candidates) {
					for (uint64_t edgeKey : candidate.edges) {
						relevantEdgeSet.insert(edgeKey);
					}
				}
				for (uint32_t fi : neighborhoodFaces) {
					const auto& vid = faces[fi].vid;
					relevantEdgeSet.insert(edgeKey64(vid[0], vid[1]));
					relevantEdgeSet.insert(edgeKey64(vid[1], vid[2]));
					relevantEdgeSet.insert(edgeKey64(vid[2], vid[0]));
				}

				std::unordered_map<uint64_t, int> outsideCount;
				for (uint32_t fi = 0; fi < nFaces; ++fi) {
					if (neighborhoodFaces.count(fi)) continue;
					const auto& vid = faces[fi].vid;
					uint64_t e0 = edgeKey64(vid[0], vid[1]);
					uint64_t e1 = edgeKey64(vid[1], vid[2]);
					uint64_t e2 = edgeKey64(vid[2], vid[0]);
					if (relevantEdgeSet.count(e0)) outsideCount[e0]++;
					if (relevantEdgeSet.count(e1)) outsideCount[e1]++;
					if (relevantEdgeSet.count(e2)) outsideCount[e2]++;
				}

				bool outsideInvalid = false;
				for (uint64_t edgeKey : relevantEdgeSet) {
					if (outsideCount[edgeKey] > 2) {
						outsideInvalid = true;
						break;
					}
				}
				if (outsideInvalid) {
					continue;
				}

				for (uint32_t ci = 0; ci < candidates.size(); ++ci) {
					for (uint64_t edgeKey : candidates[ci].edges) {
						edgeToCandidates[edgeKey].push_back(ci);
					}
				}

				std::vector<uint64_t> relevantEdges(relevantEdgeSet.begin(), relevantEdgeSet.end());
				std::vector<uint8_t> selected(candidates.size(), 0);
				std::unordered_map<uint64_t, int> localCount;
				int selectedCount = 0;
				int bestCount = std::numeric_limits<int>::max();
				std::vector<uint32_t> bestSelection;
				int searchSteps = 0;
				const int kMaxSearchSteps = 20000;

				auto edgeStatus = [&](uint64_t edgeKey) -> int {
					int outside = outsideCount[edgeKey];
					int local = localCount[edgeKey];
					if (outside == 2) {
						return local == 0 ? 0 : -1;
					}
					if (outside == 1) {
						if (local == 1) return 0;
						if (local > 1) return -1;
						return 1;
					}
					if (local == 0 || local == 2) return 0;
					if (local > 2) return -1;
					return 2;
				};

				std::function<void()> dfs;
				dfs = [&]() {
					if (++searchSteps > kMaxSearchSteps) {
						return;
					}
					if (selectedCount >= bestCount) {
						return;
					}

					uint64_t chosenEdge = 0;
					std::vector<uint32_t> chosenOptions;
					bool haveChosenEdge = false;

					for (uint64_t edgeKey : relevantEdges) {
						int status = edgeStatus(edgeKey);
						if (status < 0) {
							return;
						}
						if (status == 0) {
							continue;
						}

						std::vector<uint32_t> options;
						const auto itCand = edgeToCandidates.find(edgeKey);
						if (itCand != edgeToCandidates.end()) {
							for (uint32_t candidateIdx : itCand->second) {
								if (selected[candidateIdx]) continue;
								bool canUse = true;
								for (uint64_t candidateEdge : candidates[candidateIdx].edges) {
									if (outsideCount[candidateEdge] + localCount[candidateEdge] + 1 > 2) {
										canUse = false;
										break;
									}
								}
								if (canUse) {
									options.push_back(candidateIdx);
								}
							}
						}

						if (!haveChosenEdge || options.size() < chosenOptions.size()) {
							haveChosenEdge = true;
							chosenEdge = edgeKey;
							chosenOptions = std::move(options);
							if (chosenOptions.size() <= 1) {
								break;
							}
						}
					}

					if (!haveChosenEdge) {
						bestCount = selectedCount;
						bestSelection.clear();
						for (uint32_t candidateIdx = 0; candidateIdx < candidates.size(); ++candidateIdx) {
							if (selected[candidateIdx]) {
								bestSelection.push_back(candidateIdx);
							}
						}
						return;
					}

					if (chosenOptions.empty()) {
						return;
					}

					for (uint32_t candidateIdx : chosenOptions) {
						selected[candidateIdx] = 1;
						selectedCount++;
						for (uint64_t candidateEdge : candidates[candidateIdx].edges) {
							localCount[candidateEdge]++;
						}
						dfs();
						for (uint64_t candidateEdge : candidates[candidateIdx].edges) {
							localCount[candidateEdge]--;
						}
						selectedCount--;
						selected[candidateIdx] = 0;
					}
				};

				dfs();

				if (bestSelection.empty() || bestCount > static_cast<int>(neighborhoodFaces.size())) {
					continue;
				}

				Geometry rebuilt;
				rebuilt.planes = current.planes;
				rebuilt.hasPlanes = current.hasPlanes;
				rebuilt.data = current.data;
				for (uint32_t fi = 0; fi < nFaces; ++fi) {
					if (neighborhoodFaces.count(fi)) continue;
					rebuilt.AddFace(faces[fi].verts[0], faces[fi].verts[1], faces[fi].verts[2], faces[fi].pId);
				}
				for (uint32_t candidateIdx : bestSelection) {
					const CandidateTri& candidate = candidates[candidateIdx];
					Vec a = candidate.verts[0];
					Vec b = candidate.verts[1];
					Vec c = candidate.verts[2];
					Vec triNormal = glm::cross(b - a, c - a);
					if (glm::length(triNormal) <= 1e-12) {
						continue;
					}
					if (glm::dot(triNormal, candidate.normal) < 0.0) {
						std::swap(b, c);
					}
					rebuilt.AddFace(a, b, c, candidate.pId);
				}

				if (rebuilt.numFaces == 0) {
					continue;
				}

				MeshInfo infoRebuilt = meshCleanup::isMeshWatertight(rebuilt);
				if (!TopologyStrictlyImprovedWithoutRegression(currentInfo, infoRebuilt) &&
					!MeshPenaltyImprovedNoRegression(currentInfo, infoRebuilt)) {
					continue;
				}

				current = std::move(rebuilt);
				current.hasPlanes = false;
				currentInfo = infoRebuilt;
				removedFacesTotal += std::max<int>(0, static_cast<int>(neighborhoodFaces.size()) - bestCount);
				return true;
			}

			return false;
		};

		if (tryClosedNeighborhoodRepair()) {
			continue;
		}

		std::vector<CandidateFaceRemoval> singleFaceCandidates;
		std::unordered_set<uint32_t> seenSingleFaces;
		std::vector<CandidatePair> pairCandidates;
		for (const auto& [edge, owners] : edgeFaces) {
			if (owners.size() < 3 || owners.size() > 4) continue;

			double edgeLen = glm::length(canonicalPos[edge.v1] - canonicalPos[edge.v0]);
			double allowedThicknessSq = std::max(maxMembraneThicknessSq, edgeLen * edgeLen * 1e-4);
			double maxOwnerArea = 0.0;
			for (uint32_t faceIdx : owners) {
				maxOwnerArea = std::max(maxOwnerArea, faces[faceIdx].area);
			}

			if (maxOwnerArea > 1e-12) {
				for (uint32_t faceIdx : owners) {
					const FaceInfo& face = faces[faceIdx];
					if (face.area < 1e-12) continue;
					const bool isCoplanarOverlapFace = faceCoplanarCrowdedHits[faceIdx] >= 2;
					const bool isTinySliverFace =
						face.minAltitude <= maxAltitude &&
						face.area <= maxOwnerArea * 0.25;
					if (!(isCoplanarOverlapFace || isTinySliverFace)) continue;
					if (!seenSingleFaces.insert(faceIdx).second) continue;
					singleFaceCandidates.push_back({ faceIdx, face.area, face.minAltitude, faceCoplanarCrowdedHits[faceIdx] });
				}
			}

			for (size_t ia = 0; ia < owners.size(); ++ia) {
				for (size_t ib = ia + 1; ib < owners.size(); ++ib) {
					uint32_t faceA = owners[ia];
					uint32_t faceB = owners[ib];
					const FaceInfo& a = faces[faceA];
					const FaceInfo& b = faces[faceB];
					if (a.area < 1e-12 || b.area < 1e-12) continue;
					if (glm::dot(a.normal, b.normal) > oppositeNormalDot) continue;
					if (a.minAltitude > maxAltitude || b.minAltitude > maxAltitude) continue;

					uint32_t offA = UINT32_MAX;
					uint32_t offB = UINT32_MAX;
					for (uint32_t vid : a.vid) {
						if (vid != edge.v0 && vid != edge.v1) {
							offA = vid;
							break;
						}
					}
					for (uint32_t vid : b.vid) {
						if (vid != edge.v0 && vid != edge.v1) {
							offB = vid;
							break;
						}
					}
					if (offA == UINT32_MAX || offB == UINT32_MAX) continue;

					Vec d = canonicalPos[offA] - canonicalPos[offB];
					double thicknessSq = glm::dot(d, d);
					if (thicknessSq > allowedThicknessSq) continue;

					pairCandidates.push_back({ faceA, faceB, thicknessSq, a.area + b.area });
				}
			}
		}

		if (singleFaceCandidates.empty() && pairCandidates.empty()) {
			break;
		}

		std::sort(singleFaceCandidates.begin(), singleFaceCandidates.end(), [](const CandidateFaceRemoval& a, const CandidateFaceRemoval& b) {
			if (a.coplanarCrowdedHits != b.coplanarCrowdedHits) return a.coplanarCrowdedHits > b.coplanarCrowdedHits;
			if (a.area != b.area) return a.area < b.area;
			if (a.minAltitude != b.minAltitude) return a.minAltitude < b.minAltitude;
			return a.face < b.face;
		});

		std::sort(pairCandidates.begin(), pairCandidates.end(), [](const CandidatePair& a, const CandidatePair& b) {
			if (a.thicknessSq != b.thicknessSq) return a.thicknessSq < b.thicknessSq;
			if (a.pairArea != b.pairArea) return a.pairArea < b.pairArea;
			if (a.faceA != b.faceA) return a.faceA < b.faceA;
			return a.faceB < b.faceB;
		});

		bool acceptedCandidate = false;
		auto tryRemoval = [&](const std::unordered_set<uint32_t>& removeFaces, uint32_t removedFacesCount) -> bool {
			Geometry rebuilt;
			rebuilt.planes = current.planes;
			rebuilt.hasPlanes = current.hasPlanes;
			rebuilt.data = current.data;
			for (uint32_t i = 0; i < nFaces; ++i) {
				if (removeFaces.count(i)) continue;
				rebuilt.AddFace(faces[i].verts[0], faces[i].verts[1], faces[i].verts[2], faces[i].pId);
			}
			if (rebuilt.numFaces == 0) return false;

			MeshInfo infoRebuilt = meshCleanup::isMeshWatertight(rebuilt);
			const bool improvesStrictly = MeshPenaltyImprovedNoRegression(currentInfo, infoRebuilt);
			const bool reducesNonManifoldWithBoundedCrack =
				currentInfo.numOpenEdges > 0 &&
				infoRebuilt.numNonManifoldEdges < currentInfo.numNonManifoldEdges &&
				infoRebuilt.numOpenEdges <= currentInfo.numOpenEdges + 4;
			if (!(improvesStrictly || reducesNonManifoldWithBoundedCrack)) {
				return false;
			}

			current = std::move(rebuilt);
			current.hasPlanes = false;
			currentInfo = infoRebuilt;
			removedFacesTotal += removedFacesCount;
			return true;
		};

		for (const auto& candidate : singleFaceCandidates) {
			if (tryRemoval({ candidate.face }, 1)) {
				acceptedCandidate = true;
				break;
			}
		}

		if (!acceptedCandidate) {
			for (const auto& candidate : pairCandidates) {
				if (tryRemoval({ candidate.faceA, candidate.faceB }, 2)) {
					acceptedCandidate = true;
					break;
				}
			}
		}

		if (!acceptedCandidate) {
			break;
		}
	}

	if (removedFacesTotal == 0) {
		return 0;
	}

	workingMesh = std::move(current);
	workingMesh.hasPlanes = false;
	meshInfoResult = currentInfo;
	return removedFacesTotal;
}

void meshCleanup::removeTempFiles() {
	if (g_debugDumpDirectory.empty()) {
		return;
	}

	const std::filesystem::path dir = g_debugDumpDirectory;
	std::error_code ec;
	if (!std::filesystem::exists(dir, ec) || ec) {
		return;
	}

	for (const auto& entry : std::filesystem::directory_iterator(dir, ec)) {
		if (ec) {
			break;
		}
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

// Post-boolean cleanup -- lightweight, subtractive only.
// Called after each Subtract/Union. Only removes bad faces, never adds.
void meshCleanup::PostBooleanOperationMeshCleanup(fuzzybools::Geometry& input) {
	if (input.numFaces == 0 || input.numFaces > 2000) return;
	auto infoEntry = isMeshWatertight(input);
	if (infoEntry.watertight) return;

	// Early zero-volume membrane detection: if the mesh is not watertight
	// and has near-zero signed volume, it is a degenerate membrane artifact
	// (e.g. from a boolean subtract that fully consumed the first operand).
	// Clear it to prevent cascading errors through subsequent boolean ops.
	// Use centroid-relative computation so the result is origin-independent.
	{
		Vec centroid(0.0);
		for (uint32_t i = 0; i < input.numFaces; i++) {
			Face f = input.GetFace(i);
			centroid += input.GetPoint(f.i0);
			centroid += input.GetPoint(f.i1);
			centroid += input.GetPoint(f.i2);
		}
		centroid /= static_cast<double>(input.numFaces * 3);

		double signedVol = 0.0;
		for (uint32_t i = 0; i < input.numFaces; i++) {
			Face f = input.GetFace(i);
			Vec a = input.GetPoint(f.i0) - centroid;
			Vec b = input.GetPoint(f.i1) - centroid;
			Vec c = input.GetPoint(f.i2) - centroid;
			signedVol += glm::dot(a, glm::cross(b, c)) / 6.0;
		}
		constexpr double MEMBRANE_VOLUME_THRESHOLD = 1e-6;
		if (std::abs(signedVol) < MEMBRANE_VOLUME_THRESHOLD) {
			input = fuzzybools::Geometry();
			return;
		}
	}

	fuzzybools::Geometry working = input;
	MeshInfo cur = infoEntry;
	const uint32_t maxAllowedFaces = input.numFaces * 2;

	for (int iter = 0; iter < 3 && !cur.watertight; ++iter) {
		MeshInfo next = cur;
		RemoveDegeneratedTriangles(working, "p1", cur, next); cur = next;
		RemoveDisconnectedFragments(working, "p2", cur, next); cur = next;
		RemoveThinMembranes(working, "p3", cur, next); cur = next;
		RemoveTinyBoundaryBridgeFaces(working, "p4", cur, next); cur = next;

		if (!cur.watertight && working.numFaces <= maxAllowedFaces) {
			Geometry beforeAdditive = working;
			MeshInfo beforeAdditiveInfo = cur;
			PatchCoplanarHoles(working, "p5", cur, next);
			if (working.numFaces > maxAllowedFaces) {
				working = std::move(beforeAdditive);
				cur = beforeAdditiveInfo;
			}
			else {
				cur = next;
			}
		}

		if (!cur.watertight && cur.numNonManifoldEdges > 0) {
			RemoveOpposedEdgeMembranes(working, "p6", cur, next); cur = next;
		}
	}

	auto infoExit = cur;
	if (MeshPenaltyImprovedNoRegression(infoEntry, infoExit)) {
		input = std::move(working);
		input.hasPlanes = false;
	}
}
