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
                    firstOperand = fuzzybools::Subtract(firstOperand, secondGeom);
                }
                else if (op == "UNION")
                {
                    firstOperand = fuzzybools::Union(firstOperand, secondGeom);
                }
            }
        }
        finalResult.AddGeometry(firstOperand);
        return finalResult;
    }
}
