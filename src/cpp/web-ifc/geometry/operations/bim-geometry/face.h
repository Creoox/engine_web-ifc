#pragma once

#include "../boolean-utils/geometry.h" // fuzzybools::Face

namespace bimGeometry {
    // bimGeometry::Face was a struct { int i0, i1, i2, pId; } which is
    // field-for-field equivalent to fuzzybools::Face. They are the same
    // concept, so we alias rather than duplicating.
    using Face = fuzzybools::Face;
}
