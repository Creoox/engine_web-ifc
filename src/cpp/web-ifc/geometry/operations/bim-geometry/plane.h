#pragma once

#include "../boolean-utils/geometry.h" // fuzzybools::SimplePlane

namespace bimGeometry
{
    // bimGeometry::Plane had { double distance; Vec normal; size_t id; }
    // plus IsEqualTo. After the id field was added to fuzzybools::SimplePlane,
    // the two are layout- and method- equivalent, so we alias.
    using Plane = fuzzybools::SimplePlane;
}
