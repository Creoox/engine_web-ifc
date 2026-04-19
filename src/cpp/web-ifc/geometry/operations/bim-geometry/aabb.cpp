#include <glm/glm.hpp>
#include "aabb.h"

namespace bimGeometry {

    Buffers AABB::GetBuffers()
    {
        Buffers buffers;

        buffers.AddPoint(glm::dvec3(max.x,max.y,min.z));
        buffers.AddPoint(glm::dvec3(max.x,min.y,min.z));
        buffers.AddPoint(glm::dvec3(max.x,min.y,max.z));
        buffers.AddPoint(glm::dvec3(max.x,max.y,max.z));
        buffers.AddPoint(glm::dvec3(min.x,max.y,min.z));
        buffers.AddPoint(glm::dvec3(min.x,min.y,min.z));
        buffers.AddPoint(glm::dvec3(min.x,min.y,max.z));
        buffers.AddPoint(glm::dvec3(min.x,max.y,max.z));

        buffers.AddTri(0, 1, 3);
        buffers.AddTri(3, 1, 2);
        buffers.AddTri(5, 2, 1);
        buffers.AddTri(2, 5, 6);
        buffers.AddTri(7, 0, 4);
        buffers.AddTri(3, 0, 7);
        buffers.AddTri(7, 4, 5);
        buffers.AddTri(5, 6, 7);
        buffers.AddTri(6, 7, 3);
        buffers.AddTri(6, 2, 3);
        buffers.AddTri(5, 1, 4);
        buffers.AddTri(1, 0, 4);

        return buffers;
    }

    void AABB::SetValues(double minX, double minY, double minZ, double maxX, double maxY, double maxZ) {
        min = glm::dvec3(minX, minY, minZ);
        max = glm::dvec3(maxX, maxY, maxZ);
        center = glm::dvec3((minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2);
    }
}
