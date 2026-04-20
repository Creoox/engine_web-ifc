#include <vector>
#include <algorithm>
#include <glm/glm.hpp>
#include "curve.h"
#include "buffers.h"

#pragma once

namespace bimGeometry {
	struct Profile
	{
        uint16_t pType = 0;
        double width = 0.0;
        double depth = 0.0;
        double thickness = 0.0;
        double flangeThickness = 0.0;
        bool hasFillet = false;
        double filletRadius = 0.0;
        double radius = 0.0;
        double slope = 0.0;
        uint16_t numSegments = 0;
        std::vector<double> placement;
        Curve profile;

        void SetValues(uint16_t _pType, double _width, double _depth, double _webThickness, double _flangeThickness, bool _hasFillet, double _filletRadius, double _radius, double _slope, uint16_t _numSegments, std::vector<double> _placement);
        Buffers GetBuffers();
    };
}