#pragma once

#include <cstdint>
#include <filesystem>
#include "web-ifc/geometry/operations/boolean-utils/geometry.h"
#include "web-ifc/geometry/representation/IfcGeometry.h"

namespace meshCleanup
{
	struct MeshInfo {
		bool watertight = false;        // true if there are no boundary or non-manifold edges
		uint32_t numOpenEdges = 0;      // edges shared by exactly 1 face (actual holes)
		uint32_t numNonManifoldEdges = 0; // edges shared by more than 2 faces
		uint32_t numTotalEdges = 0;
		uint32_t numUniqueVertices = 0;
		uint32_t numFaces = 0;
	};

	static inline bool TopologyStrictlyImproved(const MeshInfo& before, const MeshInfo& after) {  // keep this here!
		return after.numNonManifoldEdges < before.numNonManifoldEdges ||
			(after.numNonManifoldEdges == before.numNonManifoldEdges &&
				after.numOpenEdges < before.numOpenEdges);
	}
	static inline bool TopologyNotWorse(const MeshInfo& before, const MeshInfo& after) {  // keep this here!
		return after.numNonManifoldEdges <= before.numNonManifoldEdges &&
			after.numOpenEdges <= before.numOpenEdges;
	}
	static inline bool TopologyStrictlyImprovedWithoutRegression(const MeshInfo& before, const MeshInfo& after) {  // keep this here!
		return TopologyNotWorse(before, after) &&
			(after.numNonManifoldEdges < before.numNonManifoldEdges ||
				after.numOpenEdges < before.numOpenEdges);
	}
	static inline uint64_t MeshPenaltyScore(const MeshInfo& info) {
		return static_cast<uint64_t>(info.numNonManifoldEdges) * 8ull +
			static_cast<uint64_t>(info.numOpenEdges) * 2ull;
	}
	static inline bool MeshPenaltyImprovedNoRegression(const MeshInfo& before, const MeshInfo& after) {
		return after.numNonManifoldEdges <= before.numNonManifoldEdges &&
			after.numOpenEdges <= before.numOpenEdges &&
			MeshPenaltyScore(after) < MeshPenaltyScore(before);
	}

	MeshInfo isMeshWatertight(const fuzzybools::Geometry& geom);
	void SetDebugDumpDirectory(const std::filesystem::path& dir);

	uint32_t RemoveDegeneratedTriangles(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t RemoveDisconnectedFragments(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t ResolveTJunctions(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t PatchCoplanarHoles(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t RemoveTinyBoundaryBridgeFaces(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t RemoveThinMembranes(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult, double thinThresholdOverride = 0);
	uint32_t RemoveOpposedEdgeMembranes(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);

	// Stage-A dedup: removes exact same-winding duplicate faces (keeps one copy) and faces whose
	// vertices collapse at watertight accounting precision (1e-7). The kept surface is identical,
	// so this is safe to run on operands (BoolProcess ingress) and on chain-exit results. Do NOT
	// call it on results that feed further booleans mid-chain: fuzzy-bools outputs are chaotically
	// sensitive to input changes (validated in v5 runs engine1-engine3).
	uint32_t RemoveExactDuplicateFaces(fuzzybools::Geometry& input);

	// Orientation heal for boolean operands: a closed shell whose total signed
	// volume is negative is wound inside-out (inward normals). The fuzzy-bools
	// classifier derives inside/outside from face winding, so an inside-out
	// operand inverts every verdict and the boolean output is garbage (Allplan
	// IFCTRIANGULATEDFACESET walls arrive like this). Flips every triangle and
	// returns the number flipped; no-op on open or correctly wound meshes.
	uint32_t FlipIfInsideOut(fuzzybools::Geometry& input);

	void PostBooleanOperationMeshCleanup(fuzzybools::Geometry& input);

	void DumpDebugGeometry(const webifc::geometry::IfcGeometry& webifcGeom, const std::string& filename);
	void DumpDebugGeometry(const fuzzybools::Geometry& geom, const std::string& filename);
	void removeTempFiles();
}
