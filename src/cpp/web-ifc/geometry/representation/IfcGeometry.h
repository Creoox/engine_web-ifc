/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.  */

// Represents a single piece of IFC Geometry

#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <glm/glm.hpp>
#include "web-ifc/geometry/operations/boolean-utils/geometry.h"
#include "../operations/bim-geometry/utils.h"

#include "geometry.h"

namespace webifc::geometry {

	constexpr int VERTEX_FORMAT_SIZE_FLOATS = fuzzybools::VERTEX_FORMAT_SIZE_FLOATS;

	struct IfcGeometry : fuzzybools::Geometry
	{
		bool halfSpace = false;
		std::vector<IfcGeometry>  part;
		Vec halfSpaceX = Vec(1, 0, 0);
		Vec halfSpaceY = Vec(0, 1, 0);
		Vec halfSpaceZ = Vec(0, 0, 1);
		Vec halfSpaceOrigin = Vec(0, 0, 0);
		Vec normalizationCenter = Vec(0, 0, 0);

		void ReverseFaces();
		void AddPart(IfcGeometry geom);
		void AddPart(fuzzybools::Geometry geom);
		void AddGeometry(fuzzybools::Geometry geom, glm::dmat4 trans = glm::dmat4(1), double scx = 1, double scy = 1, double scz = 1, glm::dvec3 origin = glm::dvec3(0, 0, 0));
		void MergeGeometry(fuzzybools::Geometry geom);
		uintptr_t GetVertexData();
		uint32_t GetVertexDataSize();
		uintptr_t GetIndexData();
		uint32_t GetIndexDataSize();
		SweptDiskSolid GetSweptDiskSolid();
		glm::dmat4 Normalize();
		SweptDiskSolid sweptDiskSolid;
	private:
		void ReverseFace(uint32_t index);
		bool normalized = false;
	};
}