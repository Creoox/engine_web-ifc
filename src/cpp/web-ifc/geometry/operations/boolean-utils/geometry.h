#pragma once

#include <vector>
#include <glm/glm.hpp>

#include "math.h"
#include "aabb.h"

namespace fuzzybools
{
	constexpr int VERTEX_FORMAT_SIZE_FLOATS = 6;
	constexpr double reconstructTolerance = 1.0E-01;

	struct SimplePlane
	{
		double distance;
		Vec normal;
		size_t id;

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

		void BuildFromVectors(std::vector<double> &d, std::vector<uint32_t> &i)
		{
			vertexData = d;
			indexData = i;

			numPoints = indexData.size();
			numFaces = indexData.size() / 3;
		}

		inline void buildPlanes()
		{
			if (!hasPlanes)
			{
				std::vector<double> storedVertexData = vertexData;
				auto GetStoredPoint = [&](size_t index) -> glm::dvec3
					{
						return glm::dvec3(
							storedVertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0],
							storedVertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1],
							storedVertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2]);
					};

				uint32_t PLANE_REFIT_ITERATIONS = 3;
				for (uint32_t r = 0; r < PLANE_REFIT_ITERATIONS; r++)
				{
					planes.clear();
					planeData.clear();

					for (size_t i = 0; i < numFaces; i++)
					{
						planeData.push_back(UINT32_MAX);
					}

					Vec centroid = Vec(0, 0, 0);

					for (size_t i = 0; i < numFaces; i++)
					{
						Face f = GetFace(i);

						auto a = GetPoint(f.i0);
						auto b = GetPoint(f.i1);
						auto c = GetPoint(f.i2);

						centroid = centroid + (a + b + c) / 3.0;
					}

					centroid /= numFaces;

					for (size_t i = 0; i < numFaces; i++) {
						Face f = GetFace(i);

						auto a = GetPoint(f.i0);
						auto b = GetPoint(f.i1);
						auto c = GetPoint(f.i2);

						glm::dvec3 norm;

						if (computeSafeNormal(a, b, c, norm, EPS_SMALL)) {
							double da = glm::dot(norm, a - centroid);
							double db = glm::dot(norm, b - centroid);
							double dc = glm::dot(norm, c - centroid);

							size_t id = AddPlane(norm, (da + db + dc) / 3.0);
							planeData[i] = id;
							hasPlanes = true;
						}
						else {
							// if other planes set hasPlanes to true, we can not leave faces with invalid plane IDs
							glm::dvec3 defaultNorm(0,0,1);
							size_t id = AddPlane(defaultNorm, 0);
							planeData[i] = id;
						}
					}

					for (size_t i = 0; i < numFaces; i++) {
						Face f = GetFace(i);

						auto a = GetPoint(f.i0);
						auto b = GetPoint(f.i1);
						auto c = GetPoint(f.i2);

						if (f.pId != UINT32_MAX)
						{
							SimplePlane p = planes[f.pId];

							double da = glm::dot(p.normal, a - centroid);
							double db = glm::dot(p.normal, b - centroid);
							double dc = glm::dot(p.normal, c - centroid);

							da = p.distance - da;
							db = p.distance - db;
							dc = p.distance - dc;

							glm::dvec3 va = a + p.normal * da;
							glm::dvec3 vb = b + p.normal * db;
							glm::dvec3 vc = c + p.normal * dc;

							glm::dvec3 dsa = GetStoredPoint(f.i0) - va;
							glm::dvec3 dsb = GetStoredPoint(f.i1) - vb;
							glm::dvec3 dsc = GetStoredPoint(f.i2) - vc;

							double fa = glm::length(dsa) / reconstructTolerance;
							double fb = glm::length(dsb) / reconstructTolerance;
							double fc = glm::length(dsc) / reconstructTolerance;

							if (fa > 1)
							{
								fa = glm::length(dsa) / fa;
								dsa = glm::normalize(dsa) * fa;
								va = va + dsa;
							}
							if (fb > 1)
							{
								fb = glm::length(dsb) / fb;
								dsb = glm::normalize(dsb) * fb;
								vb = vb + dsb;
							}
							if (fc > 1)
							{
								fc = glm::length(dsc) / fc;
								dsc = glm::normalize(dsc) * fc;
								vc = vc + dsc;
							}
							SetPoint(va.x, va.y, va.z, f.i0);
							SetPoint(vb.x, vb.y, vb.z, f.i1);
							SetPoint(vc.x, vc.y, vc.z, f.i2);
						}
					}
				}

				for (size_t i = 0; i < numFaces; i++)
				{
					Face f = GetFace(i);

					if (f.pId > -1)
					{
						auto a = GetPoint(f.i0);
						auto b = GetPoint(f.i1);
						auto c = GetPoint(f.i2);

						double da = glm::dot(planes[f.pId].normal, a);

						planes[f.pId].distance = da;
					}
				}
			}
			// TODO: Remove unused planes
		}

		inline void AddPoint(glm::dvec4 &pt, glm::dvec3 &n)
		{
			glm::dvec3 p = pt;
			AddPoint(p, n);
		}

		inline void SetPoint(double x, double y, double z, size_t index)
		{
			vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 0] = x;
			vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 1] = y;
			vertexData[index * VERTEX_FORMAT_SIZE_FLOATS + 2] = z;
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

		//void Geometry::AddPoint(const glm::dvec3& pt, const glm::dvec3& n)
		//{
		//	vertexData.push_back(pt.x);
		//	vertexData.push_back(pt.y);
		//	vertexData.push_back(pt.z);

		//	vertexData.push_back(n.x);
		//	vertexData.push_back(n.y);
		//	vertexData.push_back(n.z);

		//	numPoints += 1;
		//}

		inline void AddPoint(const glm::dvec3 &pt, const glm::dvec3 &n)
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

		inline void AddFace(uint32_t a, uint32_t b, uint32_t c, uint32_t pId = UINT32_MAX)
		{
			//			indexData.reserve((numFaces + 1) * 3);
			//			indexData[numFaces * 3 + 0] = a;
			//			indexData[numFaces * 3 + 1] = b;
			//			indexData[numFaces * 3 + 2] = c;
			indexData.push_back(a);
			indexData.push_back(b);
			indexData.push_back(c);
			planeData.push_back(pId);

			double area = areaOfTriangle(GetPoint(a), GetPoint(b), GetPoint(c));

			glm::dvec3 normal;
			//			if (!computeSafeNormal(GetPoint(a), GetPoint(b), GetPoint(c), normal, EPS_SMALL))
			if (!computeSafeNormal(GetPoint(a), GetPoint(b), GetPoint(c), normal, toleranceAddFace))
			{
				// bail out, zero area triangle
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

		inline AABB GetFaceBox(size_t index) const
		{
			AABB aabb;
			aabb.index = static_cast<uint32_t>(index);

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

		inline size_t AddPlane( const glm::dvec3& normal, double d)
		{
			for (SimplePlane& plane : planes)
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

		void AddGeometry(Geometry geom)
		{
			for (uint32_t i = 0; i < geom.numFaces; i++)
			{
				Face f = geom.GetFace(i);
				glm::dvec3 a = geom.GetPoint(f.i0);
				glm::dvec3 b = geom.GetPoint(f.i1);
				glm::dvec3 c = geom.GetPoint(f.i2);
				AddFace(a, b, c);
			}
			uint32_t planeDataOffset = geom.planes.size();
			for (uint32_t i = 0; i < geom.planeData.size(); i++)
			{
				planeData.push_back(planeDataOffset + geom.planeData[i]);
			}

			for (uint32_t i = 0; i < geom.planes.size(); i++)
			{
				planes.push_back(geom.planes[i]);
			}
		}

		void GetCenterExtents(glm::dvec3 &center, glm::dvec3 &extents) const
		{
			glm::dvec3 min(DBL_MAX, DBL_MAX, DBL_MAX);
			glm::dvec3 max(-DBL_MAX, -DBL_MAX, -DBL_MAX);

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
namespace fb = fuzzybools;