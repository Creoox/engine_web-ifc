#pragma once

#include <web-ifc/geometry/operations/bim-geometry/geometry.h>

namespace fuzzybools
{
	bool isMeshWatertight(const Geometry& geom);

	void CleanNonManifoldShells(fuzzybools::Geometry& input);
}
