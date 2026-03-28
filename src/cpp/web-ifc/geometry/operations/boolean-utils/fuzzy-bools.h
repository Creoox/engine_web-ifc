#pragma once

#include "geometry.h"
#include "shared-position.h"
#include "clip-mesh.h"
#include <web-ifc/geometry/operations/meshCleanup.h>

#include <queue>
#include <array>
#include <tuple>
#include <functional>
#include <algorithm>

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

	inline Geometry Subtract(const Geometry &A, const Geometry &B)
	{
		auto bboxA = A.GetAABB();
		auto bboxB = B.GetAABB();
		if (!bboxA.intersects(bboxB)) {
			return A;
		}
#ifdef _DEBUG
		if (A.vertexData.size() == 720) {
			int wait = 0;
		}
#endif
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B, false);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(A, B, sp, false);

#ifdef CSG_DEBUG_OUTPUT
//	DumpGeometry(geom, L"Post-normalize.obj");
#endif

		auto result = fuzzybools::clipSubtract(geom, bvh1, bvh2);
		meshCleanup::PostBooleanOperationMeshCleanup(result);
		return result;
	}

	inline Geometry Union(const Geometry &A, const Geometry &B)
	{
		fuzzybools::SharedPosition sp;
		sp.Construct(A, B, true);

		auto bvh1 = fuzzybools::MakeBVH(A);
		auto bvh2 = fuzzybools::MakeBVH(B);

		auto geom = Normalize(A, B, sp, true);

		auto result = fuzzybools::clipJoin(geom, bvh1, bvh2);
		meshCleanup::PostBooleanOperationMeshCleanup(result);
		return result;
	}
}
