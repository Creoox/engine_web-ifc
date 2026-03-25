#include "../boolean-utils/fuzzy-bools.h"

#pragma once
using namespace fuzzybools;

namespace bimGeometry 
{
    Geometry Union(Geometry firstOperator, Geometry secondOperator)
    {        
        return fuzzybools::Union(firstOperator, secondOperator);
    }

    Geometry Subtract(Geometry firstOperator, Geometry secondOperator)
    {
        return fuzzybools::Subtract(firstOperator, secondOperator);
    }

    Geometry BoolProcess(Geometry firstOperator, std::vector<Geometry> &secondGeoms, std::string op)
    {
        Geometry finalResult;

        for (auto &secondGeom : secondGeoms)
        {
            bool doit = true;
            if (secondGeom.numFaces == 0)
            {
                // bail out because we will get strange meshes
                // if this happens, probably there's an issue parsing the mesh that occurred earlier
                doit = false;
            }

            if (firstOperator.numFaces == 0 && op != "UNION")
            {
                // bail out because we will get strange meshes
                // if this happens, probably there's an issue parsing the mesh that occurred earlier
                break;
            }

            if (doit)
            {
                Geometry secondOperator;
                secondOperator = secondGeom;
                firstOperator.buildPlanes();
                secondOperator.buildPlanes();

                if (op == "DIFFERENCE")
                {
                    firstOperator = fuzzybools::Subtract(firstOperator, secondOperator);
                }
                else if (op == "UNION")
                {
                    firstOperator = fuzzybools::Union(firstOperator, secondOperator);
                }
            }
        }
        finalResult.AddGeometry(firstOperator);

        return finalResult;
    }
}