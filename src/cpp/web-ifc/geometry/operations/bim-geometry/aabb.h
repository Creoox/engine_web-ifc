#pragma once

#include "../boolean-utils/aabb.h"
#include "buffers.h"

namespace bimGeometry {

    // bimGeometry::AABB is fuzzybools::AABB extended with the two wasm
    // bindings (SetValues, GetBuffers). The dedicated intersects /
    // contains / merge / Intersect overloads that used to live here
    // duplicated the fuzzybools implementations with a slightly stricter
    // tolerance; none of the bim-geometry C++ code actually called them,
    // so they are inherited from fuzzybools::AABB (with its toleranceAABB
    // / _TOLERANCE_BOUNDING_BOX constants) and the bim-geometry variants
    // are dropped.
    struct AABB : fuzzybools::AABB
    {
        void SetValues(double minX, double minY, double minZ,
                       double maxX, double maxY, double maxZ);
        Buffers GetBuffers();
    };
}
