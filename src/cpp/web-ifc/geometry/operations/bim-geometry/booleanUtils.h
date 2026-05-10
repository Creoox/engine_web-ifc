#pragma once

#include "web-ifc/geometry/operations/boolean-utils/fuzzy-bools.h"
#include "web-ifc/geometry/operations/boolean-utils/geometry.h"

namespace bimGeometry
{
    // Sequentially apply a list of second operands to a running first operand. 
    // Kept here because a couple of places compose CSG chains outside of IfcGeometryProcessor (not the main boolean path, which
    // is in IfcGeometryProcessor::BoolProcess).
    static fuzzybools::Geometry BoolProcess(fuzzybools::Geometry firstOperand, std::vector<fuzzybools::Geometry>& secondGeoms, std::string op)
    {
        fuzzybools::Geometry finalResult;

        for (auto& secondGeom : secondGeoms)
        {
            bool doit = true;
            if (secondGeom.numFaces == 0)
            {
                // empty second operand -- skip
                doit = false;
            }
            if (firstOperand.numFaces == 0 && op != "UNION")
            {
                break;
            }
            if (doit)
            {
                firstOperand.buildPlanes();
                secondGeom.buildPlanes();
                if (op == "DIFFERENCE")
                {
                    auto budget = fuzzybools::MakeDefaultBooleanBudget(static_cast<uint64_t>(firstOperand.numFaces) + static_cast<uint64_t>(secondGeom.numFaces));
                    try
                    {
                        firstOperand = fuzzybools::Subtract(firstOperand, secondGeom, budget);
                    }
                    catch (const fuzzybools::BooleanAbortedException&)
                    {
                        continue;
                    }
                }
                else if (op == "UNION")
                {
                    auto budget = fuzzybools::MakeDefaultBooleanBudget(static_cast<uint64_t>(firstOperand.numFaces) + static_cast<uint64_t>(secondGeom.numFaces));
                    try
                    {
                        firstOperand = fuzzybools::Union(firstOperand, secondGeom, budget);
                    }
                    catch (const fuzzybools::BooleanAbortedException&)
                    {
                        continue;
                    }
                }
            }
        }
        finalResult.AddGeometry(firstOperand);
        return finalResult;
    }
}
