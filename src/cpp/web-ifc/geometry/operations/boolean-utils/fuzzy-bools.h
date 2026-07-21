#pragma once

#include <cstdlib>
#include <fstream>
#include <string>

#include "boolean-budget.h"
#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"
#include <web-ifc/geometry/operations/meshCleanup.h>

namespace fuzzybools
{
	inline void SetEpsilons(double TOLERANCE_PLANE_INTERSECTION, double TOLERANCE_PLANE_DEVIATION, double TOLERANCE_BACK_DEVIATION_DISTANCE, double TOLERANCE_INSIDE_OUTSIDE_PERIMETER, double TOLERANCE_BOUNDING_BOX, double BOOLSTATUS)
	{
		_TOLERANCE_PLANE_INTERSECTION = TOLERANCE_PLANE_INTERSECTION;
		_TOLERANCE_PLANE_DEVIATION = TOLERANCE_PLANE_DEVIATION;
		_TOLERANCE_BACK_DEVIATION_DISTANCE = TOLERANCE_BACK_DEVIATION_DISTANCE;
		_TOLERANCE_INSIDE_OUTSIDE_PERIMETER = TOLERANCE_INSIDE_OUTSIDE_PERIMETER;
		_TOLERANCE_BOUNDING_BOX = TOLERANCE_BOUNDING_BOX;
		_BOOLSTATUS = BOOLSTATUS;
	}

	inline Geometry Subtract(const Geometry &A, const Geometry &B, const BooleanBudget& budget)
	{
		budget.CheckDeadline("Subtract start");

		fuzzybools::SharedPosition sp;
		sp.Construct(A, B, false, budget);
		budget.CheckDeadline("Subtract Construct");

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);
		budget.CheckDeadline("Subtract BVH");

		auto geom = Normalize(A, B, sp, false, budget);
		budget.CheckFaceCount(geom.numFaces, "Subtract Normalize");

		if (const char* dumpDir = std::getenv("CXDIAG_CSG_NORM_DUMP"))
		{
			static int normDumpCounter = 0;
			std::string path = std::string(dumpDir) + "/norm" + std::to_string(normDumpCounter++) + ".obj";
			std::ofstream out(path);
			out.precision(17);
			for (uint32_t i = 0; i < geom.numPoints; i++)
			{
				const glm::dvec3 p = geom.GetPoint(i);
				out << "v " << p.x << " " << p.y << " " << p.z << "\n";
			}
			for (uint32_t i = 0; i < geom.numFaces; i++)
			{
				const Face f = geom.GetFace(i);
				out << "f " << (f.i0 + 1) << " " << (f.i1 + 1) << " " << (f.i2 + 1) << "\n";
			}
		}

#ifdef CSG_DEBUG_OUTPUT
//	DumpGeometry(geom, L"Post-normalize.obj");
#endif

		Geometry result = fuzzybools::clipSubtract(geom, bvh1, bvh2, budget, (uint32_t)A.planes.size());
		budget.CheckFaceCount(result.numFaces, "Subtract result");
		// Track that this mesh is a boolean result -- cleanup is allowed to
		// be more aggressive on it.
		result.mBoolOpCount = std::max(A.mBoolOpCount, B.mBoolOpCount) + 1;
		return result;
	}

	inline Geometry Subtract(const Geometry &A, const Geometry &B)
	{
		auto budget = MakeDefaultBooleanBudget(static_cast<uint64_t>(A.numFaces) + static_cast<uint64_t>(B.numFaces));
		try
		{
			return Subtract(A, B, budget);
		}
		catch (const BooleanAbortedException&)
		{
			return A;
		}
	}

	inline Geometry Union(const Geometry &A, const Geometry &B, const BooleanBudget& budget)
	{
		budget.CheckDeadline("Union start");

		fuzzybools::SharedPosition sp;
		sp.Construct(A, B, true, budget);
		budget.CheckDeadline("Union Construct");

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);
		budget.CheckDeadline("Union BVH");

		auto geom = Normalize(A, B, sp, true, budget);
		budget.CheckFaceCount(geom.numFaces, "Union Normalize");

		Geometry result = fuzzybools::clipJoin(geom, bvh1, bvh2, budget, (uint32_t)A.planes.size());
		budget.CheckFaceCount(result.numFaces, "Union result");
		result.mBoolOpCount = std::max(A.mBoolOpCount, B.mBoolOpCount) + 1;
		return result;
	}

	inline Geometry Union(const Geometry &A, const Geometry &B)
	{
		auto budget = MakeDefaultBooleanBudget(static_cast<uint64_t>(A.numFaces) + static_cast<uint64_t>(B.numFaces));
		try
		{
			return Union(A, B, budget);
		}
		catch (const BooleanAbortedException&)
		{
			return A;
		}
	}
}
