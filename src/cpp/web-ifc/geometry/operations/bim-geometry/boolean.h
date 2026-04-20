#include <vector>
#include <algorithm>
#include <glm/glm.hpp>
#include "web-ifc/geometry/operations/boolean-utils/geometry.h"
#include "buffers.h"
#include <string>

using Vec = glm::dvec3;

#pragma once

namespace bimGeometry {
    
    struct Boolean
    { 
        int type;
        std::string op;
        fuzzybools::Geometry geometry;
        std::vector<fuzzybools::Geometry> seconds;
        std::vector<std::vector<glm::dvec3>> triangles;

        Buffers GetBuffers();
        void clear();
        void SetValues(std::vector<double> triangles_, std::string op_);
        void SetSecond(std::vector<double> triangles_);
    };
}