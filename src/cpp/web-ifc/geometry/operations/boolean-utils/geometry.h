#pragma once

#include <vector>
#include <algorithm>
#include <glm/glm.hpp>

#include "math.h"
#include "aabb.h"

namespace fuzzybools
{
	constexpr int VERTEX_FORMAT_SIZE_FLOATS = 6;

	struct SimplePlane
	{
		// id is written by bimGeometry::Geometry::buildPlanes/AddPlane and
		// read by bim-geometry callers; fuzzybools itself ignores it. We
		// keep it on the shared struct so bimGeometry::Plane can alias to
		// SimplePlane without layout changes.
		size_t id = 0;
		double distance = 0.0;
		Vec normal = Vec(0.0);

		bool IsEqualTo(const Vec &n, double d)
		{
			return (equals(normal, n, toleranceVectorEquality) && equals(distance, d, TOLERANCE_SCALAR_EQUALITY));
		}
	};

	struct Face
	{
		int pId;
		int i0;
		int i1;
		int i2;
	};

	enum PrimitiveType
	{
		TRIANGLES,
		POLYLINE
	};

	struct Geometry
	{
		std::vector<float> fvertexData;
		std::vector<double> vertexData;
		std::vector<uint32_t> indexData;
		std::vector<uint32_t> planeData;
		std::vector<SimplePlane> planes;

		bool hasPlanes = false;
		uint32_t numPoints = 0;
		uint32_t numFaces = 0;
		uint32_t data = 0;
		uint32_t entityID = UINT32_MAX;
		PrimitiveType primitiveType = PrimitiveType::TRIANGLES;

		// Number of boolean operations (Subtract / Union) that produced this mesh. Incremented by fuzzybools::Subtract / Union on the result.
		// A value of 0 means the mesh is as-imported: cleanup must NOT modify it (open surface models, trees, sheets, IFCFACEBASEDSURFACEMODEL
		// etc. are legitimately open / zero-volume and must be preserved exactly). A value > 0 means the mesh is the result of CSG and may
		// carry artifacts worth cleaning aggressively.
		uint32_t mBoolOpCount = 0;

		bool isPolygon() const { return primitiveType == PrimitiveType::POLYLINE; }

		void BuildFromVectors(std::vector<double> &d, std::vector<uint32_t> &i)
		{
			vertexData = d;
			indexData = i;

			numPoints = indexData.size();
			numFaces = indexData.size() / 3;
		}

		inline void AddPoint(glm::dvec4 &pt, glm::dvec3 &n)
		{
			glm::dvec3 p = pt;
			AddPoint(p, n);
		}

		// const overload for callers built against bimGeometry::Geometry.
		inline void AddPoint(const glm::dvec3& pt, const glm::dvec3& n)
		{
			glm::dvec3 ptCopy = pt;
			glm::dvec3 nCopy = n;
			AddPoint(ptCopy, nCopy);
		}

		AABB GetAABB() const
		{
			AABB aabb;

			for (uint32_t i = 0; i < numPoints; i++)
			{
				aabb.min = glm::min(aabb.min, GetPoint(i));
				aabb.max = glm::max(aabb.max, GetPoint(i));
			}

			return aabb;
		}

		inline void AddPoint(glm::dvec3 &pt, glm::dvec3 &n)
		{
			vertexData.push_back(pt.x);
			vertexData.push_back(pt.y);
			vertexData.push_back(pt.z);

			vertexData.push_back(n.x);
			vertexData.push_back(n.y);
			vertexData.push_back(n.z);

			if (std::isnan(pt.x) || std::isnan(pt.y) || std::isnan(pt.z))
			{
				if (messages)
				{
					printf("NaN in geom!\n");
				}
			}

			if (std::isnan(n.x) || std::isnan(n.y) || std::isnan(n.z))
			{
				if (messages)
				{
					printf("NaN in geom!\n");
				}
			}

			numPoints += 1;
		}

		inline void AddFace(glm::dvec3 a, glm::dvec3 b, glm::dvec3 c, uint32_t pId = UINT32_MAX)
		{
			glm::dvec3 normal;

			double area = areaOfTriangle(a, b, c);
			if (!computeSafeNormal(a, b, c, normal, toleranceAddFace))
			{
				// bail out, zero area triangle
				if (messages)
				{
					printf("zero triangle, AddFace(vec, vec, vec)\n");
				}
				return;
			}

			AddPoint(a, normal);
			AddPoint(b, normal);
			AddPoint(c, normal);

			AddFace(numPoints - 3, numPoints - 2, numPoints - 1, pId);
		}

		// Copy-faithful variant: keeps zero-area (degenerate) triangles instead of dropping
		// them. Copy operations (flatten transform, batch concatenation) must not change the
		// face set -- the direct triangulation path (indexed AddFace) keeps degenerates, and
		// silently dropping them in copies tears watertight shells into open T-junction seams
		// (test61 subtractor class). The kernel already skips degenerate faces internally.
		inline void AddFaceKeepDegenerate(const glm::dvec3& a, const glm::dvec3& b, const glm::dvec3& c, uint32_t pId = UINT32_MAX)
		{
			glm::dvec3 normal(0.0);
			computeSafeNormal(a, b, c, normal, toleranceAddFace);

			AddPoint(a, normal);
			AddPoint(b, normal);
			AddPoint(c, normal);

			AddFace(numPoints - 3, numPoints - 2, numPoints - 1, pId);
		}

		inline void AddFace(uint32_t a, uint32_t b, uint32_t c, uint32_t pId = UINT32_MAX)
		{
			//			indexData.reserve((numFaces + 1) * 3);  // TODO: check if this is faster
			//			indexData[numFaces * 3 + 0] = a;
			//			indexData[numFaces * 3 + 1] = b;
			//			indexData[numFaces * 3 + 2] = c;
			indexData.push_back(a);
			indexData.push_back(b);
			indexData.push_back(c);
			planeData.push_back(pId);

			double area = areaOfTriangle(GetPoint(a), GetPoint(b), GetPoint(c));

			glm::dvec3 normal;
			if (!computeSafeNormal(GetPoint(a), GetPoint(b), GetPoint(c), normal, toleranceAddFace))
			{
				// bail out, zero area triangle
				
				// TODO: we are not actually bailing out here. We should either remove the degenerate face or not add it in the first place. For now we just log it.
				if (messages)
				{
					printf("zero triangle, AddFace(int, int, int)\n");
				}
			}

			numFaces++;
		}

		inline Face GetFace(size_t index) const
		{
			Face f;
			f.i0 = indexData[index * 3 + 0];
			f.i1 = indexData[index * 3 + 1];
			f.i2 = indexData[index * 3 + 2];
			f.pId = planeData[index];
			return f;
		}

		inline fuzzybools::AABB GetFaceBox(uint32_t index) const
		{
			fuzzybools::AABB aabb;
			aabb.index = index;

			glm::dvec3 a = GetPoint(indexData[index * 3 + 0]);
			glm::dvec3 b = GetPoint(indexData[index * 3 + 1]);
			glm::dvec3 c = GetPoint(indexData[index * 3 + 2]);

			aabb.min = glm::min(a, aabb.min);
			aabb.min = glm::min(b, aabb.min);
			aabb.min = glm::min(c, aabb.min);

			aabb.max = glm::max(a, aabb.max);
			aabb.max = glm::max(b, aabb.max);
			aabb.max = glm::max(c, aabb.max);

			aabb.center = (aabb.max + aabb.min) / 2.0;

			return aabb;
		}

		inline glm::dvec3 GetPoint(size_t index) const
		{
			return glm::dvec3(
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0],
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1],
				vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2]);
		}

		// Migrated from bimGeometry::Geometry during TODO 1 consolidation.
		inline void SetPoint(double x, double y, double z, size_t index)
		{
			vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0] = x;
			vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1] = y;
			vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2] = z;
		}

		// Merge another geometry into this one. The face set is copied verbatim
		// (including zero-area faces): dropping them here tore watertight shells
		// into open seams when operands were concatenated for batched subtraction.
		// The destination preserves the highest mBoolOpCount so the CSG
		// provenance survives merges.
		inline void AddGeometry(const Geometry& other)
		{
			mBoolOpCount = std::max(mBoolOpCount, other.mBoolOpCount);
			for (uint32_t i = 0; i < other.numFaces; i++)
			{
				Face f = other.GetFace(i);
				AddFaceKeepDegenerate(other.GetPoint(f.i0), other.GetPoint(f.i1), other.GetPoint(f.i2), UINT32_MAX);
			}
			// Plane table: preserve the offset bimGeometry::AddGeometry used.
			uint32_t planeDataOffset = static_cast<uint32_t>(planes.size());
			for (uint32_t i = 0; i < other.planeData.size(); i++)
			{
				planeData.push_back(planeDataOffset + other.planeData[i]);
			}
			for (uint32_t i = 0; i < other.planes.size(); i++)
			{
				planes.push_back(other.planes[i]);
			}
		}

		// Deduplicate-on-insert plane table.
		inline size_t AddPlane(const glm::dvec3& normal, double d)
		{
			for (auto& plane : planes)
			{
				if (plane.IsEqualTo(normal, d))
				{
					return plane.id;
				}
			}
			SimplePlane p;
			p.id = planes.size();
			p.normal = glm::normalize(normal);
			p.distance = d;
			planes.push_back(p);
			return p.id;
		}

		// Classify every face to a (merged) plane and build the face->plane table. Vertices are
		// never moved (see the NOTE inside). Iterates _PLANE_REFIT_ITERATIONS times.
		inline void buildPlanes()
		{
			if (hasPlanes) return;

			for (uint32_t r = 0; r < _PLANE_REFIT_ITERATIONS; r++)
			{
				planes.clear();
				planeData.clear();
				for (size_t i = 0; i < numFaces; i++)
				{
					planeData.push_back(UINT32_MAX);
				}

				glm::dvec3 centroid(0.0, 0.0, 0.0);
				for (size_t i = 0; i < numFaces; i++)
				{
					Face f = GetFace(i);
					centroid += (GetPoint(f.i0) + GetPoint(f.i1) + GetPoint(f.i2)) / 3.0;
				}
				if (numFaces > 0) centroid /= static_cast<double>(numFaces);

				for (size_t i = 0; i < numFaces; i++)
				{
					Face f = GetFace(i);
					glm::dvec3 a = GetPoint(f.i0);
					glm::dvec3 b = GetPoint(f.i1);
					glm::dvec3 c = GetPoint(f.i2);
					glm::dvec3 norm;
					if (computeSafeNormal(a, b, c, norm, EPS_SMALL))
					{
						double da = glm::dot(norm, a - centroid);
						double db = glm::dot(norm, b - centroid);
						double dc = glm::dot(norm, c - centroid);
						size_t id = AddPlane(norm, (da + db + dc) / 3.0);
						planeData[i] = static_cast<uint32_t>(id);
						hasPlanes = true;
					}
				}

				// NOTE: the historical "planar refit" vertex-snapping that lived here was removed.
				// It pulled every face's PRIVATE copies of its corners onto that face's own merged
				// plane, so copies of one shared vertex drifted up to the plane-merge tolerance
				// apart -- tearing watertight operands into open T-junction seams before every
				// boolean (test61 class: closed chain intermediates re-measured with hundreds of
				// open edges after buildPlanes). Position-consistent variants (per-position plane
				// intersection projection, with pruning or a trust region) were measured as well:
				// every mutation of the operand vertices redraws fragment-classification luck on
				// the chaos-sensitive suite files, and none beat simply not moving vertices at
				// all. The kernel only needs the face->plane table; vertices stay untouched.
			}

			for (size_t i = 0; i < numFaces; i++)
			{
				Face f = GetFace(i);
				if (f.pId > -1)
				{
					double da = glm::dot(planes[f.pId].normal, GetPoint(f.i0));
					planes[f.pId].distance = da;
				}
			}
		}

		void GetCenterExtents(glm::dvec3 &center, glm::dvec3 &extents) const
		{
			glm::dvec3 min(DBL_MAX, DBL_MAX, DBL_MAX);
			glm::dvec3 max(DBL_MIN, DBL_MIN, DBL_MIN);

			for (size_t i = 0; i < numPoints; i++)
			{
				auto pt = GetPoint(i);
				min = glm::min(min, pt);
				max = glm::max(max, pt);
			}

			extents = (max - min);
			center = min + extents / 2.0;
		}

		Geometry Normalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			Geometry newGeom;

			double scale = std::max(extents.x, std::max(extents.y, extents.z)) / 10.0;

			for (size_t i = 0; i < numFaces; i++)
			{
				auto face = GetFace(i);
				auto pa = GetPoint(face.i0);
				auto pb = GetPoint(face.i1);
				auto pc = GetPoint(face.i2);

				auto a = (pa - center) / scale;
				auto b = (pb - center) / scale;
				auto c = (pc - center) / scale;

				// std::cout << areaOfTriangle(pa, pb, pc) << std::endl;

				newGeom.AddFace(pa, pb, pc, face.pId);
			}

			return newGeom;
		}

		Geometry DeNormalize(glm::dvec3 center, glm::dvec3 extents) const
		{
			Geometry newGeom;

			double scale = std::max(extents.x, std::max(extents.y, extents.z)) / 10.0;

			for (size_t i = 0; i < numFaces; i++)
			{
				auto face = GetFace(i);
				auto pa = GetPoint(face.i0);
				auto pb = GetPoint(face.i1);
				auto pc = GetPoint(face.i2);

				// std::cout << areaOfTriangle(pa, pb, pc) << std::endl;

				auto a = pa * scale + center;
				auto b = pb * scale + center;
				auto c = pc * scale + center;

				newGeom.AddFace(a, b, c, face.pId);
			}

			return newGeom;
		}

		bool IsEmpty()
		{
			return vertexData.empty();
		}

		double Volume(const glm::dmat4 &trans = glm::dmat4(1))
		{
			double totalVolume = 0;

			for (uint32_t i = 0; i < numFaces; i++)
			{
				Face f = GetFace(i);

				glm::dvec3 a = trans * glm::dvec4(GetPoint(f.i0), 1);
				glm::dvec3 b = trans * glm::dvec4(GetPoint(f.i1), 1);
				glm::dvec3 c = trans * glm::dvec4(GetPoint(f.i2), 1);

				glm::dvec3 norm;

				//				if (computeSafeNormal(a, b, c, norm))
				if (computeSafeNormal(a, b, c, norm, EPS_SMALL))
				{
					double area = areaOfTriangle(a, b, c);
					double height = glm::dot(norm, a);

					double tetraVolume = area * height / 3;

					totalVolume += tetraVolume;
				}
			}

			return totalVolume;
		}
	};
}
