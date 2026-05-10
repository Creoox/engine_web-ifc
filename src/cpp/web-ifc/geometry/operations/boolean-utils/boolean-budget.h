#pragma once

#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace fuzzybools
{
    inline constexpr double CSG_DEFAULT_BUDGET_FLOOR_S = 10.0;
    inline constexpr double CSG_DEFAULT_BUDGET_PER_FACE_S = 0.001;
    inline constexpr double CSG_DEFAULT_BUDGET_CEILING_S = 600.0;
    inline constexpr uint32_t CSG_DEFAULT_INPROGRESS_FACE_FACTOR = 50;

    struct BooleanAbortedException : std::runtime_error
    {
        explicit BooleanAbortedException(const std::string& reason)
            : std::runtime_error(reason)
        {
        }
    };

    struct BooleanBudget
    {
        using Clock = std::chrono::steady_clock;

        Clock::time_point deadline = Clock::time_point::max();
        uint64_t inputFaces = 0;
        uint64_t maxInProgressFaces = 0;

        static BooleanBudget Unlimited(uint64_t inputFaces = 0)
        {
            BooleanBudget budget;
            budget.inputFaces = inputFaces;
            return budget;
        }

        static BooleanBudget FromSeconds(uint64_t inputFaces, double seconds, uint32_t inProgressFaceFactor)
        {
            BooleanBudget budget;
            budget.inputFaces = inputFaces;
            budget.maxInProgressFaces = inputFaces * static_cast<uint64_t>(inProgressFaceFactor);
            budget.deadline = Clock::now() + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(seconds));
            return budget;
        }

        void CheckDeadline(const char* phase) const
        {
            if (deadline == Clock::time_point::max())
            {
                return;
            }

            if (Clock::now() > deadline)
            {
                throw BooleanAbortedException(std::string("timeout in ") + PhaseName(phase));
            }
        }

        void CheckFaceCount(uint64_t faces, const char* phase) const
        {
            if (maxInProgressFaces == 0)
            {
                return;
            }

            if (faces > maxInProgressFaces)
            {
                throw BooleanAbortedException(
                    std::string("face explosion in ") + PhaseName(phase) +
                    ": faces=" + std::to_string(faces) +
                    " max=" + std::to_string(maxInProgressFaces) +
                    " inputFaces=" + std::to_string(inputFaces));
            }
        }

    private:
        static const char* PhaseName(const char* phase)
        {
            return phase != nullptr ? phase : "CSG";
        }
    };

    inline double DefaultBooleanBudgetSeconds(uint64_t inputFaces)
    {
        double seconds = CSG_DEFAULT_BUDGET_FLOOR_S + static_cast<double>(inputFaces) * CSG_DEFAULT_BUDGET_PER_FACE_S;
        if (seconds < CSG_DEFAULT_BUDGET_FLOOR_S)
        {
            return CSG_DEFAULT_BUDGET_FLOOR_S;
        }
        if (seconds > CSG_DEFAULT_BUDGET_CEILING_S)
        {
            return CSG_DEFAULT_BUDGET_CEILING_S;
        }
        return seconds;
    }

    inline BooleanBudget MakeDefaultBooleanBudget(uint64_t inputFaces)
    {
        return BooleanBudget::FromSeconds(inputFaces, DefaultBooleanBudgetSeconds(inputFaces), CSG_DEFAULT_INPROGRESS_FACE_FACTOR);
    }
}
