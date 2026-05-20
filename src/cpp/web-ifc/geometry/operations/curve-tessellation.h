/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace webifc::geometry
{
    constexpr double CURVE_TESSELLATION_PI = 3.141592653589793238462643383279502884;
    constexpr double CURVE_TESSELLATION_EPSILON = 1.0E-12;

    struct CurveTessellationSettings
    {
        CurveTessellationSettings() = default;
        explicit CurveTessellationSettings(uint16_t segments)
            : circleSegments(segments)
        {
        }

        uint16_t circleSegments = 12;
        uint16_t minArcSegments = 3;
        double maxSegmentLength = 0.0;
        double maxSagittaError = 0.0;
        uint32_t maxAdaptiveSegments = 4096;
        double lengthUnitScale = 1.0;
    };

    inline double NormalizeCurveTessellationScale(double lengthUnitScale)
    {
        if (!std::isfinite(lengthUnitScale) || lengthUnitScale <= CURVE_TESSELLATION_EPSILON)
        {
            return 1.0;
        }

        return lengthUnitScale;
    }

    inline int CapAdaptiveCurveSegments(int segments, int legacySegments, const CurveTessellationSettings& settings)
    {
        if (settings.maxAdaptiveSegments == 0 || segments <= legacySegments)
        {
            return segments;
        }

        int cap = static_cast<int>(std::min<uint32_t>(
            settings.maxAdaptiveSegments,
            static_cast<uint32_t>(std::numeric_limits<int>::max())));

        return std::min(segments, std::max(legacySegments, cap));
    }

    inline int ComputeCurveSegmentsByLength(double length, const CurveTessellationSettings& settings, int legacySegments, int minimumSegments = 1)
    {
        minimumSegments = std::max(1, minimumSegments);
        legacySegments = std::max(legacySegments, minimumSegments);
        int segments = legacySegments;

        double scaledLength = std::abs(length) * NormalizeCurveTessellationScale(settings.lengthUnitScale);
        if (std::isfinite(scaledLength)
            && settings.maxSegmentLength > CURVE_TESSELLATION_EPSILON
            && scaledLength > CURVE_TESSELLATION_EPSILON)
        {
            segments = std::max(segments, static_cast<int>(std::ceil(scaledLength / settings.maxSegmentLength)));
        }

        return std::max(minimumSegments, CapAdaptiveCurveSegments(segments, legacySegments, settings));
    }

    inline int ComputeArcSegments(double radius, double sweepRad, const CurveTessellationSettings& settings, int minimumSegments = 1)
    {
        double absSweep = std::abs(sweepRad);
        minimumSegments = std::max(1, minimumSegments);

        if (!std::isfinite(absSweep) || absSweep <= CURVE_TESSELLATION_EPSILON)
        {
            return minimumSegments;
        }

        int legacySegments = std::max(
            minimumSegments,
            static_cast<int>(static_cast<double>(settings.circleSegments) * absSweep / (2.0 * CURVE_TESSELLATION_PI)));

        int segments = legacySegments;
        double scaledRadius = std::abs(radius) * NormalizeCurveTessellationScale(settings.lengthUnitScale);

        if (std::isfinite(scaledRadius) && scaledRadius > CURVE_TESSELLATION_EPSILON)
        {
            if (settings.maxSegmentLength > CURVE_TESSELLATION_EPSILON)
            {
                double arcLength = scaledRadius * absSweep;
                segments = std::max(segments, static_cast<int>(std::ceil(arcLength / settings.maxSegmentLength)));
            }

            if (settings.maxSagittaError > CURVE_TESSELLATION_EPSILON && settings.maxSagittaError < scaledRadius)
            {
                double cosine = 1.0 - settings.maxSagittaError / scaledRadius;
                cosine = std::clamp(cosine, -1.0, 1.0);

                double maxAngle = 2.0 * std::acos(cosine);
                if (std::isfinite(maxAngle) && maxAngle > CURVE_TESSELLATION_EPSILON)
                {
                    segments = std::max(segments, static_cast<int>(std::ceil(absSweep / maxAngle)));
                }
            }
        }

        return std::max(minimumSegments, CapAdaptiveCurveSegments(segments, legacySegments, settings));
    }

    inline double PositiveAngleDifference(double startRad, double endRad)
    {
        double diff = std::fmod(endRad - startRad, 2.0 * CURVE_TESSELLATION_PI);
        if (diff < 0.0)
        {
            diff += 2.0 * CURVE_TESSELLATION_PI;
        }

        return diff;
    }

    inline double SweepAngleThroughMidpoint(double startRad, double midRad, double endRad)
    {
        double startToMid = PositiveAngleDifference(startRad, midRad);
        double startToEnd = PositiveAngleDifference(startRad, endRad);

        if (startToMid <= startToEnd)
        {
            return startToEnd;
        }

        return -(2.0 * CURVE_TESSELLATION_PI - startToEnd);
    }
}
