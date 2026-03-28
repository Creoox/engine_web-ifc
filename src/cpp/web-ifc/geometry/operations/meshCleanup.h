#pragma once

#include <cstdint>
#include <web-ifc/geometry/operations/bim-geometry/geometry.h>

namespace meshCleanup
{
	struct MeshWatertightInfo {
		bool watertight = false;        // true if no open edges (no holes)
		uint32_t numOpenEdges = 0;      // edges shared by exactly 1 face (actual holes)
		uint32_t numTotalEdges = 0;
		uint32_t numUniqueVertices = 0;
		uint32_t numFaces = 0;
	};

	MeshWatertightInfo isMeshWatertight(const fuzzybools::Geometry& geom);

	void PostBooleanOperationMeshCleanup(fuzzybools::Geometry& input);
}
