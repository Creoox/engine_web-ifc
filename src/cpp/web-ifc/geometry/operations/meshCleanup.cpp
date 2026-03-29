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
#ifdef DUMP_CSG_MESHES
	static void DumpDebugGeometry(const fuzzybools::Geometry& geom, const std::string& filename) {
		webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(geom);
		webifc::io::DumpIfcGeometry(inputWebIfc, filename);
	}
#endif

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

	static double Cross2D(const std::array<double, 2>& a, const std::array<double, 2>& b, const std::array<double, 2>& c) {
		return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
	}

	static bool PointInTriangle2D(const std::array<double, 2>& p,
		const std::array<double, 2>& a,
		const std::array<double, 2>& b,
		const std::array<double, 2>& c) {
		const double eps = 1e-9;
		double c1 = Cross2D(a, b, p);
		double c2 = Cross2D(b, c, p);
		double c3 = Cross2D(c, a, p);
		bool hasNeg = (c1 < -eps) || (c2 < -eps) || (c3 < -eps);
		bool hasPos = (c1 > eps) || (c2 > eps) || (c3 > eps);
		return !(hasNeg && hasPos);
	}

	static bool OnSegment2D(const std::array<double, 2>& a,
		const std::array<double, 2>& b,
		const std::array<double, 2>& p) {
		const double eps = 1e-9;
		return std::min(a[0], b[0]) - eps <= p[0] && p[0] <= std::max(a[0], b[0]) + eps &&
			std::min(a[1], b[1]) - eps <= p[1] && p[1] <= std::max(a[1], b[1]) + eps;
	}

	static bool SegmentsIntersect2D(const std::array<double, 2>& a0,
		const std::array<double, 2>& a1,
		const std::array<double, 2>& b0,
		const std::array<double, 2>& b1) {
		const double eps = 1e-9;
		double d1 = Cross2D(a0, a1, b0);
		double d2 = Cross2D(a0, a1, b1);
		double d3 = Cross2D(b0, b1, a0);
		double d4 = Cross2D(b0, b1, a1);

		bool intersectsProper =
			((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps)) &&
			((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps));
		if (intersectsProper) return true;

		if (std::abs(d1) <= eps && OnSegment2D(a0, a1, b0)) return true;
		if (std::abs(d2) <= eps && OnSegment2D(a0, a1, b1)) return true;
		if (std::abs(d3) <= eps && OnSegment2D(b0, b1, a0)) return true;
		if (std::abs(d4) <= eps && OnSegment2D(b0, b1, a1)) return true;
		return false;
	}

	static bool TriangleIntersectsPolygon2D(const std::array<std::array<double, 2>, 3>& tri,
		const std::vector<std::array<double, 2>>& polygon,
		const std::array<double, 2>& polyMin,
		const std::array<double, 2>& polyMax) {
		if (polygon.size() < 3) return false;

		std::array<double, 2> triMin{ std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
		std::array<double, 2> triMax{ -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };
		for (const auto& p : tri) {
			triMin[0] = std::min(triMin[0], p[0]);
			triMin[1] = std::min(triMin[1], p[1]);
			triMax[0] = std::max(triMax[0], p[0]);
			triMax[1] = std::max(triMax[1], p[1]);
		}
		if (triMax[0] < polyMin[0] || triMin[0] > polyMax[0] ||
			triMax[1] < polyMin[1] || triMin[1] > polyMax[1]) {
			return false;
		}

		for (const auto& p : tri) {
			if (PointInPolygon2D(p, polygon)) {
				return true;
			}
		}

		for (const auto& p : polygon) {
			if (PointInTriangle2D(p, tri[0], tri[1], tri[2])) {
				return true;
			}
		}

		for (int i = 0; i < 3; i++) {
			const auto& ta = tri[i];
			const auto& tb = tri[(i + 1) % 3];
			for (size_t j = 0, k = polygon.size() - 1; j < polygon.size(); k = j++) {
				if (SegmentsIntersect2D(ta, tb, polygon[k], polygon[j])) {
					return true;
				}
			}
		}

		return false;
	}

	static bool TraceSimpleLoop(const std::vector<std::pair<uint32_t, uint32_t>>& edges,
		std::vector<uint32_t>& loopOut) {
		loopOut.clear();
		if (edges.size() < 3) return false;

		std::unordered_map<uint32_t, std::vector<uint32_t>> adjacency;
		for (const auto& edge : edges) {
			adjacency[edge.first].push_back(edge.second);
			adjacency[edge.second].push_back(edge.first);
		}

		for (const auto& [vertex, neighbours] : adjacency) {
			if (neighbours.size() != 2) {
				return false;
			}
		}

		uint32_t start = edges.front().first;
		uint32_t prev = UINT32_MAX;
		uint32_t cur = start;

		for (size_t safety = 0; safety <= edges.size() + 1; safety++) {
			loopOut.push_back(cur);

			const auto& neighbours = adjacency[cur];
			uint32_t next = neighbours[0];
			if (next == prev) next = neighbours[1];

			prev = cur;
			cur = next;
			if (cur == start) {
				return loopOut.size() >= 3;
			}
		}

		loopOut.clear();
		return false;
	}

	static uint32_t RetriangulateHostPlaneOpeningPatches(Geometry& geom) {
		const uint32_t nFaces = geom.numFaces;
		if (nFaces == 0) return 0;

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

		auto getCell = [&](const Vec& p) -> GridKey {
			return { static_cast<int64_t>(std::floor(p.x / cellSize)),
					static_cast<int64_t>(std::floor(p.y / cellSize)),
					static_cast<int64_t>(std::floor(p.z / cellSize)) };
		};

		struct CanonFace {
			std::array<uint32_t, 3> vid;
			Vec normal;
			double area;
			uint32_t pId;
		};

		std::unordered_map<GridKey, std::vector<std::pair<uint32_t, Vec>>, GridKeyHash> vtxGrid;
		std::vector<Vec> canonPos;
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
						if (it != vtxGrid.end()) {
							for (const auto& [id, pos] : it->second) {
								Vec d = p - pos;
								if (glm::dot(d, d) < cellSizeSq) {
									return id;
								}
							}
						}
					}
			uint32_t id = nextVid++;
			vtxGrid[center].emplace_back(id, p);
			canonPos.push_back(p);
			return id;
		};

		std::vector<CanonFace> faces(nFaces);
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

		std::unordered_map<EKey, uint32_t, EKeyHash> allEdgeCount;
		for (const auto& face : faces) {
			for (int e = 0; e < 3; e++) {
				allEdgeCount[makeEKey(face.vid[e], face.vid[(e + 1) % 3])]++;
			}
		}

		std::array<double, 3> dominantNormalArea = { 0.0, 0.0, 0.0 };
		double meshAxisMin = std::numeric_limits<double>::max();
		double meshAxisMax = -std::numeric_limits<double>::max();
		for (const auto& face : faces) {
			dominantNormalArea[0] += face.area * std::abs(face.normal.x);
			dominantNormalArea[1] += face.area * std::abs(face.normal.y);
			dominantNormalArea[2] += face.area * std::abs(face.normal.z);
			for (uint32_t vid : face.vid) {
				meshAxisMin = std::min(meshAxisMin, canonPos[vid].x);
				meshAxisMax = std::max(meshAxisMax, canonPos[vid].x);
			}
		}

		int dominantAxis = 0;
		if (dominantNormalArea[1] > dominantNormalArea[dominantAxis]) dominantAxis = 1;
		if (dominantNormalArea[2] > dominantNormalArea[dominantAxis]) dominantAxis = 2;
		meshAxisMin = std::numeric_limits<double>::max();
		meshAxisMax = -std::numeric_limits<double>::max();
		for (const auto& point : canonPos) {
			double axisValue = dominantAxis == 0 ? point.x : (dominantAxis == 1 ? point.y : point.z);
			meshAxisMin = std::min(meshAxisMin, axisValue);
			meshAxisMax = std::max(meshAxisMax, axisValue);
		}

		auto axisComponent = [&](const Vec& v) -> double {
			return dominantAxis == 0 ? v.x : (dominantAxis == 1 ? v.y : v.z);
		};

		const double sideTol = std::max(1e-3 * 4.0, toleranceVectorEquality * 2.0);
		auto vertexOnSide = [&](uint32_t vid, bool minSide) -> bool {
			double axisValue = axisComponent(canonPos[vid]);
			return minSide ? std::abs(axisValue - meshAxisMin) <= sideTol : std::abs(axisValue - meshAxisMax) <= sideTol;
		};

		auto faceOnSide = [&](const CanonFace& face, bool minSide) -> bool {
			double faceMin = std::min({ axisComponent(canonPos[face.vid[0]]), axisComponent(canonPos[face.vid[1]]), axisComponent(canonPos[face.vid[2]]) });
			double faceMax = std::max({ axisComponent(canonPos[face.vid[0]]), axisComponent(canonPos[face.vid[1]]), axisComponent(canonPos[face.vid[2]]) });
			return minSide ? faceMax <= meshAxisMin + sideTol : faceMin >= meshAxisMax - sideTol;
		};

		auto isHostAligned = [&](const CanonFace& face) -> bool {
			return glm::length(face.normal) > 0.5 && std::abs(axisComponent(face.normal)) > 0.85;
		};

		auto isOrthogonalToHost = [&](const CanonFace& face) -> bool {
			return glm::length(face.normal) > 0.5 && std::abs(axisComponent(face.normal)) < 0.35;
		};

		struct PatchCandidate {
			std::vector<uint32_t> removedFaces;
			std::vector<uint32_t> polygonVertexIds;
			std::vector<uint32_t> triangulation;
			Vec planeNormal;
			Vec normalHint;
			uint32_t planeId = 0;
		};

		std::vector<PatchCandidate> patches;
		std::vector<bool> faceReserved(nFaces, false);

		for (int side = 0; side < 2; side++) {
			bool minSide = side == 0;
			bool hasBoundaryEdgesOnSide = false;
			for (const auto& [edge, count] : allEdgeCount) {
				if (count != 1) continue;
				if (vertexOnSide(edge.v0, minSide) && vertexOnSide(edge.v1, minSide)) {
					hasBoundaryEdgesOnSide = true;
					break;
				}
			}
			if (!hasBoundaryEdgesOnSide) continue;
			std::array<double, 2> sideHostMin{ std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
			std::array<double, 2> sideHostMax{ -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };
			bool hasHostSideBounds = false;
			for (const auto& face : faces) {
				if (!isHostAligned(face)) continue;
				if (!faceOnSide(face, minSide)) continue;
				for (uint32_t vid : face.vid) {
					auto projected = ProjectTo2D(canonPos[vid], dominantAxis);
					sideHostMin[0] = std::min(sideHostMin[0], projected[0]);
					sideHostMin[1] = std::min(sideHostMin[1], projected[1]);
					sideHostMax[0] = std::max(sideHostMax[0], projected[0]);
					sideHostMax[1] = std::max(sideHostMax[1], projected[1]);
				}
				hasHostSideBounds = true;
			}
			if (!hasHostSideBounds) continue;

			std::unordered_set<EKey, EKeyHash> holeEdgeSet;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (!isOrthogonalToHost(faces[i])) continue;
				for (int e = 0; e < 3; e++) {
					uint32_t v0 = faces[i].vid[e];
					uint32_t v1 = faces[i].vid[(e + 1) % 3];
					if (vertexOnSide(v0, minSide) && vertexOnSide(v1, minSide)) {
						holeEdgeSet.insert(makeEKey(v0, v1));
					}
				}
			}

			if (holeEdgeSet.empty()) continue;

			std::unordered_map<uint32_t, std::vector<uint32_t>> holeAdj;
			for (const auto& edge : holeEdgeSet) {
				holeAdj[edge.v0].push_back(edge.v1);
				holeAdj[edge.v1].push_back(edge.v0);
			}

			std::unordered_set<uint32_t> visitedHoleVertices;
			for (const auto& entry : holeAdj) {
				uint32_t startVertex = entry.first;
				if (visitedHoleVertices.count(startVertex)) continue;

				std::vector<uint32_t> componentVertices;
				std::queue<uint32_t> queue;
				queue.push(startVertex);
				visitedHoleVertices.insert(startVertex);

				while (!queue.empty()) {
					uint32_t cur = queue.front();
					queue.pop();
					componentVertices.push_back(cur);
					for (uint32_t nb : holeAdj[cur]) {
						if (visitedHoleVertices.insert(nb).second) {
							queue.push(nb);
						}
					}
				}

				std::unordered_set<uint32_t> componentSet(componentVertices.begin(), componentVertices.end());
				std::vector<std::pair<uint32_t, uint32_t>> componentEdges;
				for (const auto& edge : holeEdgeSet) {
					if (componentSet.count(edge.v0) && componentSet.count(edge.v1)) {
						componentEdges.push_back({ edge.v0, edge.v1 });
					}
				}

				std::vector<uint32_t> holeLoop;
				if (!TraceSimpleLoop(componentEdges, holeLoop)) continue;
				if (holeLoop.size() < 3) continue;

				std::vector<std::array<double, 2>> holeProjected;
				holeProjected.reserve(holeLoop.size());
				std::array<double, 2> holeMin{ std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
				std::array<double, 2> holeMax{ -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max() };
				for (uint32_t vid : holeLoop) {
					auto projected = ProjectTo2D(canonPos[vid], dominantAxis);
					holeProjected.push_back(projected);
					holeMin[0] = std::min(holeMin[0], projected[0]);
					holeMin[1] = std::min(holeMin[1], projected[1]);
					holeMax[0] = std::max(holeMax[0], projected[0]);
					holeMax[1] = std::max(holeMax[1], projected[1]);
				}
				if (std::abs(SignedArea2D(holeProjected)) <= 1e-10) continue;
				const double edgeMargin = toleranceVectorEquality * 4.0;
				if (holeMin[0] <= sideHostMin[0] + edgeMargin ||
					holeMin[1] <= sideHostMin[1] + edgeMargin ||
					holeMax[0] >= sideHostMax[0] - edgeMargin ||
					holeMax[1] >= sideHostMax[1] - edgeMargin) {
					continue;
				}

				std::vector<uint32_t> selectedFaces;
				Vec normalHint(0);
				std::unordered_map<uint32_t, double> planeArea;

				for (uint32_t i = 0; i < nFaces; i++) {
					if (!isHostAligned(faces[i])) continue;
					if (!faceOnSide(faces[i], minSide)) continue;

					std::array<std::array<double, 2>, 3> triProjected = {
						ProjectTo2D(canonPos[faces[i].vid[0]], dominantAxis),
						ProjectTo2D(canonPos[faces[i].vid[1]], dominantAxis),
						ProjectTo2D(canonPos[faces[i].vid[2]], dominantAxis)
					};

					if (!TriangleIntersectsPolygon2D(triProjected, holeProjected, holeMin, holeMax)) {
						continue;
					}

					selectedFaces.push_back(i);
					normalHint += faces[i].normal * faces[i].area;
					planeArea[faces[i].pId] += faces[i].area;
				}

				if (selectedFaces.empty()) continue;

				bool overlapsExistingPatch = false;
				for (uint32_t fi : selectedFaces) {
					if (faceReserved[fi]) {
						overlapsExistingPatch = true;
						break;
					}
				}
				if (overlapsExistingPatch) continue;

				std::unordered_map<EKey, uint32_t, EKeyHash> patchEdgeCount;
				for (uint32_t fi : selectedFaces) {
					for (int e = 0; e < 3; e++) {
						patchEdgeCount[makeEKey(faces[fi].vid[e], faces[fi].vid[(e + 1) % 3])]++;
					}
				}

				std::unordered_map<uint32_t, std::vector<uint32_t>> patchBoundaryAdj;
				for (const auto& [edge, count] : patchEdgeCount) {
					if (count != 1) continue;
					patchBoundaryAdj[edge.v0].push_back(edge.v1);
					patchBoundaryAdj[edge.v1].push_back(edge.v0);
				}
				if (patchBoundaryAdj.empty()) continue;
				bool patchBoundaryIsSimple = true;
				for (const auto& [vertex, neighbours] : patchBoundaryAdj) {
					if (neighbours.size() != 2) {
						patchBoundaryIsSimple = false;
						break;
					}
				}
				if (!patchBoundaryIsSimple) continue;

				std::unordered_set<uint32_t> visitedPatchVertices;
				std::vector<std::vector<uint32_t>> boundaryLoops;
				boundaryLoops.reserve(2);

				for (const auto& entry : patchBoundaryAdj) {
					uint32_t vertex = entry.first;
					if (visitedPatchVertices.count(vertex)) continue;

					std::vector<uint32_t> componentVerticesPatch;
					std::queue<uint32_t> patchQueue;
					patchQueue.push(vertex);
					visitedPatchVertices.insert(vertex);
					while (!patchQueue.empty()) {
						uint32_t cur = patchQueue.front();
						patchQueue.pop();
						componentVerticesPatch.push_back(cur);
						for (uint32_t nb : patchBoundaryAdj[cur]) {
							if (visitedPatchVertices.insert(nb).second) {
								patchQueue.push(nb);
							}
						}
					}

					std::unordered_set<uint32_t> patchComponentSet(componentVerticesPatch.begin(), componentVerticesPatch.end());
					std::vector<std::pair<uint32_t, uint32_t>> patchEdges;
					for (const auto& [edge, count] : patchEdgeCount) {
						if (count != 1) continue;
						if (patchComponentSet.count(edge.v0) && patchComponentSet.count(edge.v1)) {
							patchEdges.push_back({ edge.v0, edge.v1 });
						}
					}

					std::vector<uint32_t> boundaryLoop;
					if (!TraceSimpleLoop(patchEdges, boundaryLoop)) continue;
					boundaryLoops.push_back(std::move(boundaryLoop));
				}

				if (boundaryLoops.size() != 2) continue;

				std::unordered_set<uint32_t> holeLoopSet(holeLoop.begin(), holeLoop.end());
				auto loopMatchesHole = [&](const std::vector<uint32_t>& candidate) -> bool {
					if (candidate.size() != holeLoop.size()) return false;
					for (uint32_t vid : candidate) {
						if (!holeLoopSet.count(vid)) return false;
					}
					return true;
				};

				int holeBoundaryIndex = -1;
				for (int idx = 0; idx < static_cast<int>(boundaryLoops.size()); idx++) {
					if (loopMatchesHole(boundaryLoops[idx])) {
						holeBoundaryIndex = idx;
						break;
					}
				}
				if (holeBoundaryIndex < 0) continue;

				int outerBoundaryIndex = holeBoundaryIndex == 0 ? 1 : 0;
				const auto& outerLoop = boundaryLoops[outerBoundaryIndex];
				const auto& patchHoleLoop = boundaryLoops[holeBoundaryIndex];

				std::vector<std::array<double, 2>> outerProjected;
				outerProjected.reserve(outerLoop.size());
				for (uint32_t vid : outerLoop) {
					outerProjected.push_back(ProjectTo2D(canonPos[vid], dominantAxis));
				}
				double outerArea = std::abs(SignedArea2D(outerProjected));
				if (outerLoop.size() < 3 || outerArea <= 1e-10) continue;

				std::vector<std::array<double, 2>> patchHoleProjected;
				patchHoleProjected.reserve(patchHoleLoop.size());
				for (uint32_t vid : patchHoleLoop) {
					patchHoleProjected.push_back(ProjectTo2D(canonPos[vid], dominantAxis));
				}
				double patchHoleArea = std::abs(SignedArea2D(patchHoleProjected));
				if (patchHoleArea <= 1e-10 || outerArea <= patchHoleArea) continue;

				std::vector<std::vector<std::array<double, 2>>> polygon;
				std::vector<uint32_t> polygonVertexIds;
				polygon.reserve(2);
				polygonVertexIds.reserve(outerLoop.size() + patchHoleLoop.size());

				polygon.emplace_back();
				for (uint32_t vid : outerLoop) {
					polygon.back().push_back(ProjectTo2D(canonPos[vid], dominantAxis));
					polygonVertexIds.push_back(vid);
				}

				polygon.emplace_back();
				for (uint32_t vid : patchHoleLoop) {
					polygon.back().push_back(ProjectTo2D(canonPos[vid], dominantAxis));
					polygonVertexIds.push_back(vid);
				}

				std::vector<uint32_t> triangulation = mapbox::earcut<uint32_t>(polygon);
				if (triangulation.size() < 3) continue;

				uint32_t planeId = selectedFaces.front() < faces.size() ? faces[selectedFaces.front()].pId : 0;
				double bestPlaneArea = -1.0;
				for (const auto& [pId, area] : planeArea) {
					if (area > bestPlaneArea) {
						bestPlaneArea = area;
						planeId = pId;
					}
				}

				Vec planeNormal(0);
				if (dominantAxis == 0) planeNormal = Vec(minSide ? -1.0 : 1.0, 0.0, 0.0);
				if (dominantAxis == 1) planeNormal = Vec(0.0, minSide ? -1.0 : 1.0, 0.0);
				if (dominantAxis == 2) planeNormal = Vec(0.0, 0.0, minSide ? -1.0 : 1.0);
				if (glm::length(normalHint) < 1e-9) {
					normalHint = planeNormal;
				}

				for (uint32_t fi : selectedFaces) {
					faceReserved[fi] = true;
				}

				patches.push_back({
					std::move(selectedFaces),
					std::move(polygonVertexIds),
					std::move(triangulation),
					planeNormal,
					normalHint,
					planeId
					});
			}
		}

		if (patches.empty()) {
			return 0;
		}

		std::vector<bool> removeFace(nFaces, false);
		for (const auto& patch : patches) {
			for (uint32_t fi : patch.removedFaces) {
				removeFace[fi] = true;
			}
		}

		Geometry rebuilt;
		rebuilt.planes = geom.planes;
		rebuilt.hasPlanes = geom.hasPlanes;
		rebuilt.data = geom.data;

		for (uint32_t i = 0; i < nFaces; i++) {
			if (removeFace[i]) continue;
			Face f = geom.GetFace(i);
			rebuilt.AddFace(geom.GetPoint(f.i0), geom.GetPoint(f.i1), geom.GetPoint(f.i2), f.pId);
		}

		uint32_t rebuiltFaces = 0;
		for (const auto& patch : patches) {
			Vec fillNormal(0);
			for (size_t tri = 0; tri < patch.triangulation.size(); tri += 3) {
				Vec a = canonPos[patch.polygonVertexIds[patch.triangulation[tri]]];
				Vec b = canonPos[patch.polygonVertexIds[patch.triangulation[tri + 1]]];
				Vec c = canonPos[patch.polygonVertexIds[patch.triangulation[tri + 2]]];
				fillNormal += glm::cross(b - a, c - a);
			}

			bool flipWinding = glm::dot(fillNormal, patch.normalHint) < 0;
			for (size_t tri = 0; tri < patch.triangulation.size(); tri += 3) {
				Vec a = canonPos[patch.polygonVertexIds[patch.triangulation[tri]]];
				Vec b = canonPos[patch.polygonVertexIds[patch.triangulation[tri + 1]]];
				Vec c = canonPos[patch.polygonVertexIds[patch.triangulation[tri + 2]]];

				if (flipWinding) {
					rebuilt.AddFace(a, c, b, patch.planeId);
				}
				else {
					rebuilt.AddFace(a, b, c, patch.planeId);
				}
				rebuiltFaces++;
			}
		}

		if (rebuiltFaces == 0) {
			return 0;
		}

		geom = std::move(rebuilt);
		geom.hasPlanes = false;
		return rebuiltFaces;
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
			else if (count > 2) {
				++info.numNonManifoldEdges;
			}
		}

		info.watertight = (info.numOpenEdges == 0 && info.numNonManifoldEdges == 0);
		return info;
	}

	static bool TopologyStrictlyImproved(const MeshWatertightInfo& before, const MeshWatertightInfo& after) {
		return after.numNonManifoldEdges < before.numNonManifoldEdges ||
			(after.numNonManifoldEdges == before.numNonManifoldEdges &&
			 after.numOpenEdges < before.numOpenEdges);
	}

	static bool TopologyNotWorse(const MeshWatertightInfo& before, const MeshWatertightInfo& after) {
		return after.numNonManifoldEdges <= before.numNonManifoldEdges &&
			after.numOpenEdges <= before.numOpenEdges;
	}

	static bool TopologyStrictlyImprovedWithoutRegression(const MeshWatertightInfo& before,
		const MeshWatertightInfo& after) {
		return TopologyNotWorse(before, after) &&
			(after.numNonManifoldEdges < before.numNonManifoldEdges ||
			 after.numOpenEdges < before.numOpenEdges);
	}

	// ---------------------------------------------------------------
	// Remove disconnected non-manifold zero-volume components.
	// Safe to run at any point -- purely topological, no heuristics.
	// ---------------------------------------------------------------
	static uint32_t RemoveDisconnectedFragments(Geometry& workingMesh, std::string step,
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
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
			return 0;
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
	static uint32_t removeDegeneratedTriangles(Geometry& workingMesh, std::string step,
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
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
			webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup"+step+"-fail.obj");
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
	static uint32_t ResolveTJunctions(Geometry& geom, std::string step, 
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
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

		MeshWatertightInfo infoRebuilt = meshCleanup::isMeshWatertight(rebuilt);
		if (infoRebuilt.numOpenEdges < meshInfoInput.numOpenEdges) {
			geom = std::move(rebuilt);
			meshInfoResult = infoRebuilt;
			return retriangulatedFaces;
		}
		else {
			meshInfoResult = meshInfoInput;
#ifdef _DEBUG
			webifc::geometry::IfcGeometry webifcGeom = webifc::geometry::booleanManager::convertToWebIfc(geom);
			webifc::geometry::IfcGeometry geomFail = webifc::geometry::booleanManager::convertToWebIfc(rebuilt);
			webifc::io::DumpIfcGeometry(webifcGeom, "meshCleanup" + step + "-input.obj");
			webifc::io::DumpIfcGeometry(geomFail, "meshCleanup" + step + "-fail.obj");
#endif
			return 0;
		}
	}

	static uint32_t RemoveTinyBoundaryBridgeFaces(Geometry& geom, std::string step,
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
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
		MeshWatertightInfo currentInfo = meshInfoInput;
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
	static uint32_t PatchCoplanarHoles(Geometry& geom, std::string step, 
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
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
				boundaryAdj[ek.v0].push_back(ek.v1);
				boundaryAdj[ek.v1].push_back(ek.v0);
				boundaryEdgeFace[undirectedEdgeKey(ek.v0, ek.v1)] = fl[0];
			}
		}

		const double PLANE_EPS = 1e-5;

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

		const double LOOP_PLANE_GROUP_EPS = std::max(PLANE_EPS * 4.0, toleranceVectorEquality * 2.0);
		const double LOOP_PLANE_GROUP_DOT = 0.995;

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

			if (infoPatched.numOpenEdges < meshInfoInput.numOpenEdges) {
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

	static uint32_t RemoveThinMembranes(Geometry& workingMesh, std::string step, 
		const MeshWatertightInfo& meshInfoInput, MeshWatertightInfo& meshInfoResult) {
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

		constexpr double THIN_THRESHOLD = 1e-3;
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

		if (thinCount == 0 || thinCount >= nFaces * 3 / 4) {
			return 0;
		}

		const uint32_t maxRemovedFaces = std::max<uint32_t>(24u, nFaces / 12);
		if (thinCount > maxRemovedFaces) {
			return 0;
		}

		const double areaLossRatio = totalAreaBefore > 1e-9 ? removedArea / totalAreaBefore : 0.0;
		if (areaLossRatio > 0.03) {
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

		if (cleaned.numFaces == 0 || cleaned.numFaces >= nFaces) {
			return 0;
		}

		auto cleanedInfo = meshCleanup::isMeshWatertight(cleaned);
		if (cleanedInfo.numOpenEdges > meshInfoInput.numOpenEdges + 8) {
			return 0;
		}

		const bool improvesTopology = TopologyStrictlyImprovedWithoutRegression(meshInfoInput, cleanedInfo);
		const bool reducesNonManifold =
			cleanedInfo.numNonManifoldEdges < meshInfoInput.numNonManifoldEdges;
		if (!(improvesTopology || reducesNonManifold)) {
			return 0;
		}

		workingMesh = std::move(cleaned);
		workingMesh.hasPlanes = false;
		meshInfoResult = cleanedInfo;

#ifdef _DEBUG
		if (thinCount > 7 && cleanedInfo.numFaces > 87) {
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
			webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup" + step + "-removedMembranes.obj");
		}
#endif
		return thinCount;
#if 0
		// -- Shared setup for Phases B & C 
		const Geometry baseMesh = workingMesh;
		const uint32_t nFaces = baseMesh.numFaces;

		// -- Step 1: cache face vertices + per-face geometric info
		struct FV {
			Vec a, b, c;
			uint32_t pId;
			Vec center, normal;
			double area, maxEdge;
		};
		std::vector<FV> fv(nFaces);
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
		}

		// -- Step 2: spatial-hash vertex deduplication
		// Mirror SharedPosition::AddPoint: grid cell = floor(coord/cellSize), 27-cell neighbourhood search, exact Euclidean distance check.
		const double cellSize = toleranceVectorEquality * 1.5;          // align with legacy membrane cleanup
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
		BVH resultBVH = MakeBVH(baseMesh);

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

		// Keep B.2 marks separate for pass accounting, but also feed them into
		// the thin-mask so B.3 erosion can grow bridge faces the same way the
		// legacy cleanup did.
		for (uint32_t i = 0; i < nFaces; i++) {
			if (b2Marked[i]) {
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

		// B.7: host-plane membrane detection for incomplete opening cuts.
		// A failed subtract can leave a cap-like sheet on one host wall plane.
		// Detect the opening corridor from abnormal side-wall faces, then grow
		// only host-aligned faces that lie on an outer host plane and stay
		// inside that corridor. This targets the membrane itself instead of the
		// tunnel walls of a valid opening.
		{
			std::array<double, 3> dominantNormalArea = { 0.0, 0.0, 0.0 };
			for (uint32_t i = 0; i < nFaces; i++) {
				dominantNormalArea[0] += fv[i].area * std::abs(fv[i].normal.x);
				dominantNormalArea[1] += fv[i].area * std::abs(fv[i].normal.y);
				dominantNormalArea[2] += fv[i].area * std::abs(fv[i].normal.z);
			}

			int dominantAxis = 0;
			if (dominantNormalArea[1] > dominantNormalArea[dominantAxis]) dominantAxis = 1;
			if (dominantNormalArea[2] > dominantNormalArea[dominantAxis]) dominantAxis = 2;

			double dominantArea = dominantNormalArea[dominantAxis];
			double secondaryArea = 0.0;
			for (int axis = 0; axis < 3; axis++) {
				if (axis == dominantAxis) continue;
				secondaryArea = std::max(secondaryArea, dominantNormalArea[axis]);
			}

			if (dominantArea > secondaryArea * 2.0) {
				auto faceBoundaryEdges = [&](uint32_t fi) -> int {
					int count = 0;
					for (int e = 0; e < 3; e++) {
						uint32_t va = fvid[fi][e], vb = fvid[fi][(e + 1) % 3];
						EKey ek = { std::min(va, vb), std::max(va, vb) };
						auto it = edgeFaces.find(ek);
						if (it != edgeFaces.end() && it->second.size() == 1) {
							count++;
						}
					}
					return count;
				};

				auto faceNonManifoldEdges = [&](uint32_t fi) -> int {
					int count = 0;
					for (int e = 0; e < 3; e++) {
						uint32_t va = fvid[fi][e], vb = fvid[fi][(e + 1) % 3];
						EKey ek = { std::min(va, vb), std::max(va, vb) };
						auto it = edgeFaces.find(ek);
						if (it != edgeFaces.end() && it->second.size() > 2) {
							count++;
						}
					}
					return count;
				};

				auto axisComponent = [&](const Vec& v) -> double {
					return dominantAxis == 0 ? v.x : (dominantAxis == 1 ? v.y : v.z);
				};

				auto isOrthogonalToHost = [&](uint32_t fi) -> bool {
					if (glm::length(fv[fi].normal) < 0.5) return false;
					return std::abs(axisComponent(fv[fi].normal)) < 0.35;
				};

				auto isHostAligned = [&](uint32_t fi) -> bool {
					if (glm::length(fv[fi].normal) < 0.5) return false;
					return std::abs(axisComponent(fv[fi].normal)) > 0.85;
				};

				Vec corridorMin(std::numeric_limits<double>::max());
				Vec corridorMax(-std::numeric_limits<double>::max());
				std::vector<uint32_t> corridorFaces;
				auto expandBounds = [&](const Vec& p, Vec& minP, Vec& maxP) {
					minP.x = std::min(minP.x, p.x);
					minP.y = std::min(minP.y, p.y);
					minP.z = std::min(minP.z, p.z);
					maxP.x = std::max(maxP.x, p.x);
					maxP.y = std::max(maxP.y, p.y);
					maxP.z = std::max(maxP.z, p.z);
				};

				for (uint32_t i = 0; i < nFaces; i++) {
					if (!isOrthogonalToHost(i)) continue;
					if (faceBoundaryEdges(i) == 0 && faceNonManifoldEdges(i) == 0) continue;
					corridorFaces.push_back(i);
					expandBounds(fv[i].a, corridorMin, corridorMax);
					expandBounds(fv[i].b, corridorMin, corridorMax);
					expandBounds(fv[i].c, corridorMin, corridorMax);
				}

				if (!corridorFaces.empty()) {
					double meshAxisMin = std::numeric_limits<double>::max();
					double meshAxisMax = -std::numeric_limits<double>::max();
					for (uint32_t i = 0; i < nFaces; i++) {
						meshAxisMin = std::min(meshAxisMin, std::min({ axisComponent(fv[i].a), axisComponent(fv[i].b), axisComponent(fv[i].c) }));
						meshAxisMax = std::max(meshAxisMax, std::max({ axisComponent(fv[i].a), axisComponent(fv[i].b), axisComponent(fv[i].c) }));
					}

					const double corridorTol = toleranceVectorEquality * 2.0;
					const double sideTol = std::max(THIN_THRESHOLD * 4.0, toleranceVectorEquality * 2.0);

					auto withinCorridor = [&](const Vec& p) -> bool {
						if (dominantAxis != 0 && (p.x < corridorMin.x - corridorTol || p.x > corridorMax.x + corridorTol)) return false;
						if (dominantAxis != 1 && (p.y < corridorMin.y - corridorTol || p.y > corridorMax.y + corridorTol)) return false;
						if (dominantAxis != 2 && (p.z < corridorMin.z - corridorTol || p.z > corridorMax.z + corridorTol)) return false;
						return true;
					};

					auto isOnOuterHostSide = [&](uint32_t fi, bool minSide) -> bool {
						double faceMinAxis = std::min({ axisComponent(fv[fi].a), axisComponent(fv[fi].b), axisComponent(fv[fi].c) });
						double faceMaxAxis = std::max({ axisComponent(fv[fi].a), axisComponent(fv[fi].b), axisComponent(fv[fi].c) });
						if (minSide) return faceMaxAxis <= meshAxisMin + sideTol;
						return faceMinAxis >= meshAxisMax - sideTol;
					};

					auto isHostMembraneCandidate = [&](uint32_t fi, bool minSide) -> bool {
						if (!isHostAligned(fi)) return false;
						if (!isOnOuterHostSide(fi, minSide)) return false;
						return withinCorridor(fv[fi].a) && withinCorridor(fv[fi].b) && withinCorridor(fv[fi].c);
					};

					for (int side = 0; side < 2; side++) {
						std::vector<bool> visitedCandidate(nFaces, false);
						for (uint32_t corridorFace : corridorFaces) {
							for (uint32_t seed : faceAdj[corridorFace]) {
								if (visitedCandidate[seed]) continue;
								if (thinMarked[seed]) continue;
								if (!isHostMembraneCandidate(seed, side == 0)) continue;

								std::vector<uint32_t> region;
								std::queue<uint32_t> q;
								q.push(seed);
								visitedCandidate[seed] = true;
								uint32_t boundaryLinkedFaces = 0;
								uint32_t nonManifoldLinkedFaces = 0;

								while (!q.empty()) {
									uint32_t cur = q.front();
									q.pop();
									region.push_back(cur);

									bool touchesCorridor = false;
									for (uint32_t nb : faceAdj[cur]) {
										if (isOrthogonalToHost(nb)) {
											touchesCorridor = true;
											if (faceBoundaryEdges(nb) > 0) boundaryLinkedFaces++;
											if (faceNonManifoldEdges(nb) > 0) nonManifoldLinkedFaces++;
										}
										if (visitedCandidate[nb]) continue;
										if (thinMarked[nb]) continue;
										if (!isHostMembraneCandidate(nb, side == 0)) continue;
										visitedCandidate[nb] = true;
										q.push(nb);
									}

									if (!touchesCorridor && faceBoundaryEdges(cur) > 0) {
										boundaryLinkedFaces++;
									}
								}

								bool membraneRegion =
									!region.empty() &&
									region.size() < nFaces / 4 &&
									(boundaryLinkedFaces > 0 || nonManifoldLinkedFaces > 0);

								if (membraneRegion) {
									for (uint32_t fi : region) {
										thinMarked[fi] = true;
									}
								}
							}
						}
					}
				}
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
		uint32_t removedMembraneFaces = 0;
		uint32_t pass2RemovedMembraneFaces = 0;
		uint32_t pass3RemovedMembraneFaces = 0;
		for (uint32_t i = 0; i < nFaces; i++) {
			if (badComp[compId[i]]) continue;
			if (b2Marked[i]) {
				++pass2RemovedMembraneFaces;
			}
			if (b2Marked[i] || thinMarked[i]) {
				++pass3RemovedMembraneFaces;
			}
		}

		// Pass 1: remove bad components only (Phase C)
		if (removedByComp > 0 && removedByComp < nFaces) {
			Geometry pass1;
			pass1.planes = baseMesh.planes;
			pass1.hasPlanes = baseMesh.hasPlanes;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (badComp[compId[i]]) continue;
				pass1.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
			}
			pass1.data = baseMesh.data;

			auto infoPass1 = meshCleanup::isMeshWatertight(pass1);
			if (TopologyStrictlyImprovedWithoutRegression(currentInfo, infoPass1)) {
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
			pass2.planes = baseMesh.planes;
			pass2.hasPlanes = baseMesh.hasPlanes;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (badComp[compId[i]]) continue;
				if (b2Marked[i]) continue;
				pass2.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
			}
			pass2.data = baseMesh.data;
			if (pass2RemovedMembraneFaces > 0) {
				RetriangulateHostPlaneOpeningPatches(pass2);
			}

			MeshWatertightInfo infoPass2Before = meshCleanup::isMeshWatertight(pass2);
			MeshWatertightInfo infoPass2 = infoPass2Before;
			if (infoPass2Before.numOpenEdges > 0) {
				PatchCoplanarHoles(pass2, step + "b", infoPass2Before, infoPass2);
			}
			if (TopologyStrictlyImprovedWithoutRegression(currentInfo, infoPass2) ||
				(pass2RemovedMembraneFaces > 0 && TopologyNotWorse(currentInfo, infoPass2))) {
				workingMesh = std::move(pass2);
				currentInfo = infoPass2;
				removedMembraneFaces = pass2RemovedMembraneFaces;
#ifdef DUMP_CSG_MESHES
				DumpDebugGeometry(workingMesh, "meshCleanup" + step + "b-patched.obj");
#endif
			}
		}

		// Pass 3: also remove B.1+B.3+B.5 thin-marked faces
		if (thinEnabled && thinCount > 0) {
			Geometry pass3;
			pass3.planes = baseMesh.planes;
			pass3.hasPlanes = baseMesh.hasPlanes;
			for (uint32_t i = 0; i < nFaces; i++) {
				if (badComp[compId[i]]) continue;
				if (b2Marked[i]) continue;
				if (thinMarked[i]) continue;
				pass3.AddFace(fv[i].a, fv[i].b, fv[i].c, fv[i].pId);
			}
			pass3.data = baseMesh.data;
			if (pass3RemovedMembraneFaces > 0) {
				RetriangulateHostPlaneOpeningPatches(pass3);
			}

			MeshWatertightInfo infoPass3Before = meshCleanup::isMeshWatertight(pass3);
			MeshWatertightInfo infoPass3 = infoPass3Before;
			if (infoPass3Before.numOpenEdges > 0) {
				PatchCoplanarHoles(pass3, step + "c", infoPass3Before, infoPass3);
			}

			if (TopologyStrictlyImprovedWithoutRegression(currentInfo, infoPass3) ||
				(pass3RemovedMembraneFaces > 0 && TopologyNotWorse(currentInfo, infoPass3))) {
				workingMesh = std::move(pass3);
				currentInfo = infoPass3;
				removedMembraneFaces = pass3RemovedMembraneFaces;
#ifdef DUMP_CSG_MESHES
				DumpDebugGeometry(workingMesh, "meshCleanup" + step + "c-patched.obj");
#endif
			}
		}

		meshInfoResult = currentInfo;

#ifdef _DEBUG
		if (currentInfo.numOpenEdges > 0 && currentInfo.numOpenEdges >= meshInfoInput.numOpenEdges) {
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
			webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup" + step + "-fail.obj");
		}
		if (removedMembraneFaces > 7) {
			if (currentInfo.numFaces > 87) {
				webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(workingMesh);
				webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup" + step + "-removedMembranes.obj");
			}
		}
#endif
		return removedMembraneFaces;
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

	size_t cleanUpCallCount = 0;
	// Post-boolean cleanup
	void PostBooleanOperationMeshCleanup(fuzzybools::Geometry& input) {
		struct CleanupMeshCounts {
			uint32_t removedDegeneratedTriangleChanges = 0;
			uint32_t removedDisconnectedFragmentFaces = 0;
			uint32_t patchedHoleFacesPass1 = 0;
			uint32_t removedMembraneFaces = 0;
			uint32_t resolvedTJunctionFaces = 0;
			uint32_t patchedHoleFacesPass2 = 0;
			uint32_t removedDegeneratedTriangleChangesFinal = 0;
			uint32_t removedTinyBoundaryBridgeFaces = 0;
		} changeCounts;

		if (input.numFaces > 8000) {
#ifdef _DEBUG
			std::cout << "PostBooleanOperationMeshCleanup: skipping mesh with " << input.numFaces << " faces" << std::endl;
#endif
			return;
		}
		fuzzybools::Geometry workingMesh = input;
		auto meshInfoOnEntry = meshCleanup::isMeshWatertight(input);
		if (meshInfoOnEntry.watertight) {
			// TODO: test return here, and do a more aggressive membrane removal before rendering/export
			//return;
		}
#ifdef DUMP_CSG_MESHES
		// remove all E:\work\creoox\cxconverter\meshCleanup*.obj
		removeTempFiles();

		++cleanUpCallCount;
		if (meshInfoOnEntry.numFaces > 190 && cleanUpCallCount > 13) {
			int wait = 0;
		}
		
		webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(input);
		webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup1.obj");
#endif

		// 1: remove degenerated triangles
		MeshWatertightInfo meshInfoBeforeRemoveDegeneratedTriangles = meshInfoOnEntry;
		MeshWatertightInfo meshInfoAfterRemoveDegeneratedTriangles = meshInfoOnEntry;
		changeCounts.removedDegeneratedTriangleChanges =
			removeDegeneratedTriangles(workingMesh, "1", meshInfoBeforeRemoveDegeneratedTriangles, meshInfoAfterRemoveDegeneratedTriangles);
		
		// 2: remove disconnected fragments (Phase C only -- safe before hole patching)
		MeshWatertightInfo meshInfoAfterRemoveFragments = meshInfoAfterRemoveDegeneratedTriangles;
		changeCounts.removedDisconnectedFragmentFaces =
			RemoveDisconnectedFragments(workingMesh, "2", meshInfoAfterRemoveDegeneratedTriangles, meshInfoAfterRemoveFragments);

		// 3: remove thin membranes before generic hole patching so the membrane
		// detector still sees the raw boolean scars around openings.
		MeshWatertightInfo meshInfoBeforeRemoveThinMembranes = meshInfoAfterRemoveFragments;
		MeshWatertightInfo meshInfoAfterRemoveThinMembranes = meshInfoAfterRemoveFragments;
		changeCounts.removedMembraneFaces =
			RemoveThinMembranes(workingMesh, "3", meshInfoBeforeRemoveThinMembranes, meshInfoAfterRemoveThinMembranes);

		// 4: patch coplanar holes after membrane cleanup
		MeshWatertightInfo meshInfoBeforePatchCoplanarHoles = meshInfoAfterRemoveThinMembranes;
		MeshWatertightInfo meshInfoAfterPatchCoplanarHoles = meshInfoAfterRemoveThinMembranes;
		changeCounts.patchedHoleFacesPass1 =
			PatchCoplanarHoles(workingMesh, "4", meshInfoBeforePatchCoplanarHoles, meshInfoAfterPatchCoplanarHoles);
#ifdef DUMP_CSG_MESHES
		if (changeCounts.patchedHoleFacesPass1 > 0) {
			DumpDebugGeometry(workingMesh, "meshCleanup4-patched.obj");
		}
#endif

		// 5: resolve T-junctions
		MeshWatertightInfo meshInfoBeforeResolveTJunctions = meshInfoAfterPatchCoplanarHoles;
		MeshWatertightInfo meshInfoAfterResolveTJunctions = meshInfoAfterPatchCoplanarHoles;
		changeCounts.resolvedTJunctionFaces =
			ResolveTJunctions(workingMesh, "5", meshInfoBeforeResolveTJunctions, meshInfoAfterResolveTJunctions);

		// 6: PatchCoplanarHoles re-run
		MeshWatertightInfo meshInfoBeforePatchCoplanarHoles2 = meshInfoAfterResolveTJunctions;
		MeshWatertightInfo meshInfoAfterPatchCoplanarHoles2 = meshInfoAfterResolveTJunctions;
		changeCounts.patchedHoleFacesPass2 =
			PatchCoplanarHoles(workingMesh, "6", meshInfoBeforePatchCoplanarHoles2, meshInfoAfterPatchCoplanarHoles2);
#ifdef DUMP_CSG_MESHES
		if (changeCounts.patchedHoleFacesPass2 > 0) {
			DumpDebugGeometry(workingMesh, "meshCleanup6-patched.obj");
		}
#endif

		// 7: final sliver cleanup for tiny crack-producing triangles that can
		// survive or be created by the previous repair passes.
		MeshWatertightInfo meshInfoBeforeRemoveDegeneratedTrianglesFinal = meshInfoAfterPatchCoplanarHoles2;
		MeshWatertightInfo meshInfoAfterRemoveDegeneratedTrianglesFinal = meshInfoAfterPatchCoplanarHoles2;
		changeCounts.removedDegeneratedTriangleChangesFinal =
			removeDegeneratedTriangles(workingMesh, "7", meshInfoBeforeRemoveDegeneratedTrianglesFinal, meshInfoAfterRemoveDegeneratedTrianglesFinal);
#ifdef DUMP_CSG_MESHES
		if (changeCounts.removedDegeneratedTriangleChangesFinal > 0) {
			DumpDebugGeometry(workingMesh, "meshCleanup7-slivers.obj");
		}
#endif

		// 8: remove tiny boundary bridge faces that leave a single crack edge
		MeshWatertightInfo meshInfoBeforeRemoveTinyBoundaryBridgeFaces = meshInfoAfterRemoveDegeneratedTrianglesFinal;
		MeshWatertightInfo meshInfoAfterRemoveTinyBoundaryBridgeFaces = meshInfoAfterRemoveDegeneratedTrianglesFinal;
		changeCounts.removedTinyBoundaryBridgeFaces =
			RemoveTinyBoundaryBridgeFaces(workingMesh, "8", meshInfoBeforeRemoveTinyBoundaryBridgeFaces, meshInfoAfterRemoveTinyBoundaryBridgeFaces);
#ifdef DUMP_CSG_MESHES
		if (changeCounts.removedTinyBoundaryBridgeFaces > 0) {
			DumpDebugGeometry(workingMesh, "meshCleanup8-bridges.obj");
		}
#endif

		auto meshInfoOnExit = meshInfoAfterRemoveTinyBoundaryBridgeFaces;
		const bool improvedTopology = TopologyStrictlyImproved(meshInfoOnEntry, meshInfoOnExit);
		const bool removedAnyMembranes = changeCounts.removedMembraneFaces > 0;
		const bool keptMembraneImprovement =
			removedAnyMembranes &&
			TopologyNotWorse(meshInfoOnEntry, meshInfoOnExit);
		if (improvedTopology || keptMembraneImprovement) {
			input = std::move(workingMesh);
			input.hasPlanes = false;

#ifdef DUMP_CSG_MESHES
			webifc::geometry::IfcGeometry inputWebIfc = webifc::geometry::booleanManager::convertToWebIfc(input);
			webifc::io::DumpIfcGeometry(inputWebIfc, "meshCleanup9-improved.obj");
#endif
		}
	}
}
