/* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

// Implementation for IfcGeometry

#include "IfcGeometry.h"
#include "web-ifc/geometry/operations/geometryutils.h"

namespace webifc::geometry {
	void IfcGeometry::ReverseFace(uint32_t index)
	{
		fuzzybools::Face f = GetFace(index);
		indexData[index * 3 + 0] = f.i2;
		indexData[index * 3 + 1] = f.i1;
		indexData[index * 3 + 2] = f.i0;
	}

	void IfcGeometry::ReverseFaces()
	{
		for (size_t i = 0; i < numFaces; i++)
		{
			ReverseFace(i);
		}
	}

	glm::dmat4 IfcGeometry::Normalize()
	{
		glm::dvec3 center = normalizationCenter;
		if (!normalized)
		{
			glm::dvec3 extents(0, 0, 0);
			GetCenterExtents(center, extents);
			for (size_t i = 0; i < vertexData.size(); i += 6)
			{
				vertexData[i + 0] = vertexData[i + 0] - center.x;
				vertexData[i + 1] = vertexData[i + 1] - center.y;
				vertexData[i + 2] = vertexData[i + 2] - center.z;
			}
			for (size_t i = 0; i < sweptDiskSolid.axis.size(); i++)
			{
				for (size_t j = 0; j < sweptDiskSolid.axis[i].points.size(); j++)
				{
					sweptDiskSolid.axis[i].points[j].x -= center.x;
					sweptDiskSolid.axis[i].points[j].y -= center.y;
					sweptDiskSolid.axis[i].points[j].z -= center.z;
				}
			}

			for (size_t i = 0; i < sweptDiskSolid.profiles.size(); i++)
			{
				for (size_t j = 0; j < sweptDiskSolid.profiles[i].curve.points.size(); j++)
				{
					sweptDiskSolid.profiles[i].curve.points[j].x -= center.x;
					sweptDiskSolid.profiles[i].curve.points[j].y -= center.y;
					sweptDiskSolid.profiles[i].curve.points[j].z -= center.z;
				}
				for (size_t j = 0; j < sweptDiskSolid.profiles[i].holes.size(); j++)
				{
					for (size_t k = 0; k < sweptDiskSolid.profiles[i].holes[j].points.size(); k++)
					{
						sweptDiskSolid.profiles[i].holes[j].points[k].x -= center.x;
						sweptDiskSolid.profiles[i].holes[j].points[k].y -= center.y;
						sweptDiskSolid.profiles[i].holes[j].points[k].z -= center.z;
					}
				}
			}
			normalizationCenter = center;
			normalized = true;
		}
		glm::dmat4 resultMat = glm::dmat4(1.0);
		resultMat[3][0] = center[0];
		resultMat[3][1] = center[1];
		resultMat[3][2] = center[2];

		return  resultMat;
	}

	uintptr_t IfcGeometry::GetVertexData()
	{
		// unfortunately webgl can't do doubles
		if (fvertexData.size() != vertexData.size())
		{
			fvertexData.resize(vertexData.size());
			for (size_t i = 0; i < vertexData.size(); i++)
			{
				// The vector was previously copied in batches of 6, but
				// copying single entry at a time is more resilient if the 
				// underlying geometry lib changes the treatment of normals
				fvertexData[i] = vertexData[i];
			}
		}
		if (fvertexData.empty())
		{
			return 0;
		}
		return (uintptr_t)(size_t)&fvertexData[0];
	}

	uint32_t IfcGeometry::GetVertexDataSize()
	{
		return (uint32_t)fvertexData.size();
	}

	uintptr_t IfcGeometry::GetIndexData()
	{
		return (uintptr_t)(size_t)&indexData[0];
	}

	uint32_t IfcGeometry::GetIndexDataSize()
	{
		return (uint32_t)indexData.size();
	}

	SweptDiskSolid IfcGeometry::GetSweptDiskSolid()
	{
		return sweptDiskSolid;
	}

	void IfcGeometry::AddPart(IfcGeometry geom)
	{
		part.push_back(geom);
	}

	void IfcGeometry::AddPart(Geometry geom)
	{
		IfcGeometry newGeom;
		newGeom.MergeGeometry(geom);
		part.push_back(newGeom);
	}

	void IfcGeometry::AddGeometry(Geometry geom, glm::dmat4 trans, double scx, double scy, double scz, glm::dvec3 origin)
	{
		mBoolOpCount = std::max(mBoolOpCount, geom.mBoolOpCount);
		for (uint32_t i = 0; i < geom.numFaces; i++)
		{
			fuzzybools::Face f = geom.GetFace(i);
			glm::dvec3 a = geom.GetPoint(f.i0);
			glm::dvec3 b = geom.GetPoint(f.i1);
			glm::dvec3 c = geom.GetPoint(f.i2);
			if (scx != 1 || scy != 1 || scz != 1)
			{
				double aax = glm::dot(trans[0], glm::dvec4(a - origin, 1)) * scx;
				double aay = glm::dot(trans[1], glm::dvec4(a - origin, 1)) * scy;
				double aaz = glm::dot(trans[2], glm::dvec4(a - origin, 1)) * scz;
				a = origin + glm::dvec3(aax * trans[0]) + glm::dvec3(aay * trans[1]) + glm::dvec3(aaz * trans[2]);
				double bbx = glm::dot(trans[0], glm::dvec4(b - origin, 1)) * scx;
				double bby = glm::dot(trans[1], glm::dvec4(b - origin, 1)) * scy;
				double bbz = glm::dot(trans[2], glm::dvec4(b - origin, 1)) * scz;
				b = origin + glm::dvec3(bbx * trans[0]) + glm::dvec3(bby * trans[1]) + glm::dvec3(bbz * trans[2]);
				double ccx = glm::dot(trans[0], glm::dvec4(c - origin, 1)) * scx;
				double ccy = glm::dot(trans[1], glm::dvec4(c - origin, 1)) * scy;
				double ccz = glm::dot(trans[2], glm::dvec4(c - origin, 1)) * scz;
				c = origin + glm::dvec3(ccx * trans[0]) + glm::dvec3(ccy * trans[1]) + glm::dvec3(ccz * trans[2]);
			}
			// Copy operation: keep degenerate faces, dropping them tears watertight
			// shells open (see fuzzybools::Geometry::AddFaceKeepDegenerate).
			AddFaceKeepDegenerate(a, b, c);
		}
		AddPart(geom);
	}

	void IfcGeometry::MergeGeometry(Geometry geom)
	{
		mBoolOpCount = std::max(mBoolOpCount, geom.mBoolOpCount);
		for (uint32_t i = 0; i < geom.numFaces; i++)
		{
			fuzzybools::Face f = geom.GetFace(i);
			glm::dvec3 a = geom.GetPoint(f.i0);
			glm::dvec3 b = geom.GetPoint(f.i1);
			glm::dvec3 c = geom.GetPoint(f.i2);
			AddFaceKeepDegenerate(a, b, c);
		}
	}
}