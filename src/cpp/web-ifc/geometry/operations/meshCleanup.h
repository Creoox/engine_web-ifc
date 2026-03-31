#pragma once

#include <cstdint>
#include <filesystem>
#include <web-ifc/geometry/operations/bim-geometry/geometry.h>
#include <web-ifc/geometry/operations/boolean-utils/geometry.h>

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

	MeshInfo isMeshWatertight(const fuzzybools::Geometry& geom);
	void SetDebugDumpDirectory(const std::filesystem::path& dir);

	uint32_t RemoveDegeneratedTriangles(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t RemoveDisconnectedFragments(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t ResolveTJunctions(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t PatchCoplanarHoles(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);
	uint32_t RemoveTinyBoundaryBridgeFaces(fuzzybools::Geometry& input, std::string step, const MeshInfo& meshInfoInput, MeshInfo& meshInfoResult);

	void PostBooleanOperationMeshCleanup(fuzzybools::Geometry& input);
	void FullMeshCleanup(fuzzybools::Geometry& input);
	void DumpDebugGeometry(const fuzzybools::Geometry& geom, const std::string& filename);
	void removeTempFiles();
}
