#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <ranges>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <glm/glm.hpp>

#pragma warning(push)
#pragma warning(disable : 4267)
#include <CDT.h>
#pragma warning(pop)

#include "geometry.h"
#include "aabb.h"
#include "bvh.h"
#include "boolean-budget.h"
#include "util.h"
#include "svg.h"
#include "math.h"
#include "loop-finder.h"
#include "obj-exporter.h"
#include "is-inside-mesh.h"
#include "is-inside-boundary.h"

#ifndef GLM_ENABLE_EXPERIMENTAL
#define GLM_ENABLE_EXPERIMENTAL
#endif
#include <glm/gtx/norm.hpp>

using Vec2 = glm::dvec2;
using Vec3 = glm::dvec3;

// Custom hash for std::tuple<int64_t, int64_t, int64_t>
// Combine the three 64-bit integers in a reasonably collision-resistant way
template<>
struct std::hash<std::tuple<int64_t, int64_t, int64_t>>
{
    using argument_type = std::tuple<int64_t, int64_t, int64_t>;
    using result_type = std::size_t;

    result_type operator()(const argument_type& t) const noexcept
    {
        auto [x, y, z] = t;  // C++17 structured binding (preferred)

        // Simple, fast mixing - inspired by boost::hash_combine + murmur-like finalizer
        std::size_t seed = 0x517cc1b727220a95ULL;  // random 64-bit constant

        seed ^= static_cast<std::size_t>(x) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(y) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        seed ^= static_cast<std::size_t>(z) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);

        return seed;
    }
};

namespace fuzzybools
{
    struct PlaneBasis
    {
        Vec3 origin;

        Vec3 up;
        Vec3 left;
        Vec3 right;

        Vec2 project(const Vec3 &pt)
        {
            auto relative = pt - origin;
            return Vec2(glm::dot(relative, left), glm::dot(relative, right));
        }
    };

    struct ReferencePlane
    {
        size_t planeID;
        size_t pointID;
        size_t lineID;
        Vec2 location;
    };

    struct ReferenceLine
    {
        size_t lineID;
        size_t pointID;
        double location;
    };

    /*
        The vector direction is not necessarily a unit vector!  It is the vector pointing
        from the start of the line (identified with origin) to the end.  Thus,
            end = origin + direction
    */
    struct Line
    {
        size_t id;
        size_t globalID;
        Vec3 origin;
        Vec3 direction;

        Line()
        {
            static size_t idcounter = 0;
            idcounter++;
            globalID = idcounter;
        }

        /*
                Does point a lie on the line defined by vectors origin and direction?

        */
        bool IsPointOnLine(const Vec3 &a) const
        {
            /*
                        Vector d is the unit vector pointing along the line from origin towards
                            origin + direction.
            */
            Vec3 d = glm::normalize(direction);
            Vec3 v = a - origin;
            /*
                        Drop a perpendicular from point a to the line.  The quantity t
                        is the distance from origin along the line to the foot of the perpendicular,
                        denoted by p.
            */
            double t = glm::dot(v, d);
            Vec3 p = origin + t * d;
            /*
                        If point a is sufficently close to point p, then decide that point a
                        lies on the line.
            */
            return glm::distance(p, a) < tolerancePointOnLine;
        }

        /*
                The original version of GetPosOnline seemed to assume that direction is a unit vector.
                This version normalises direction first.
        */

        double GetPosOnLine(const Vec3 &pos) const
        {
            Vec3 unitDirection = glm::normalize(direction);
            return glm::dot(pos - origin, unitDirection);
        }

        Vec3 GetPosOnLine(const double dist) const
        {
            Vec3 unitDirection = glm::normalize(direction);
            return origin + dist * unitDirection;
        }

        bool IsCollinear(const Line &other) const
        {
            Vec3 unitDirection = glm::normalize(direction);
            Vec3 unitOtherDirection = glm::normalize(other.direction);
            return (equals(unitOtherDirection, unitDirection, toleranceCollinear) || equals(unitOtherDirection, -unitDirection, toleranceCollinear));
            //          return (equals(other.direction,    direction,     toleranceCollinear) || equals(other.direction,    -direction,     toleranceCollinear));
        }

        /*
                The original version of IsEqualTo compared dir and direction.  This version compares
                normalised copies of dir and direction.
        */
        bool IsEqualTo(const Vec3 &pos, const Vec3 &dir) const
        {
            // check dir
            Vec3 unitDir = glm::normalize(dir);
            Vec3 unitDirection = glm::normalize(direction);
            if (!(equals(unitDir, unitDirection, EPS_SMALL) || equals(unitDir, -unitDirection, EPS_SMALL)))
            {
                return false;
            }

            // check pos
            if (!IsPointOnLine(pos))
            {
                return false;
            }

            return true;
        }

        /*
            Insert a point (identified by id) at the given distance along the line,
            maintaining sorted order by distance. Duplicate point ids are ignored.
        */
        void AddPointToLine(double dist, size_t id)
        {
            // check existing
            for (auto &p : points)
            {
                if (p.second == id)
                    return;
            }

            // binary-search insert to maintain sorted order (O(log N) search + O(N) shift)
            // instead of push_back + full re-sort (O(N log N))
            auto it = std::lower_bound(points.begin(), points.end(), dist,
                [](const std::pair<double, size_t> &p, double val) { return p.first < val; });
            points.insert(it, std::make_pair(dist, id));
        }

        /*
            Return a lazy view of consecutive point-pair segments along the line.
            Each segment is a pair (pointId_i, pointId_{i+1}) derived from the
            sorted points list.
        */
        auto GetSegments() const
        {
            const auto makeSegments = [&](int i)
            {
                return std::make_pair(points[i - 1].second, points[i].second);
            };

            return std::views::iota(size_t(1), points.size()) | std::views::transform(makeSegments);
        }

        std::vector<std::pair<double, size_t>> points;

        std::vector<ReferencePlane> planes;
    };

    struct Point
    {
        size_t id;
        size_t globalID;
        Vec3 location3D;

        Point()
        {
            static size_t idcounter = 0;
            idcounter++;
            globalID = idcounter;
        }

        bool operator==(const Vec3 &pt)
        {
            return equals(location3D, pt, toleranceVectorEquality);
        }

        std::vector<ReferenceLine> lines;
        std::vector<ReferencePlane> planes;
    };

    struct Plane
    {
        int refPlane = -1;
        size_t id;
        size_t globalID;
        double distance;
        Vec3 normal;

        std::vector<Line> lines;
        AABB aabb;

        void AddPoint(const Vec3 &pt)
        {
            aabb.merge(pt);
        }

        Plane()
        {
            static size_t idcounter = 0;
            idcounter++;
            globalID = idcounter;
        }

        double round(double input)
        {
            input = std::fabs(input) < EPS_BIG ? 0.0 : input;
            input = std::fabs(input) < (1.0 - EPS_BIG) ? input : input > 0.0 ? 1.0
                                                                             : -1.0;
            return input;
        }

        Vec3 round(Vec3 in)
        {
            in.x = round(in.x);
            in.y = round(in.y);
            in.z = round(in.z);

            return in;
        }

        Vec3 GetDirection(Vec3 a, Vec3 b)
        {
            auto dir = b - a;
            return glm::normalize(dir);
        }

        /*
            Add a line defined by two Points to this plane. If an equivalent line
            already exists (same direction and coincident), it is reused. Both
            points are registered on the line at their projected distances.
            Returns the line's index and whether a new line was created (true)
            or an existing one was reused (false).
        */
        std::pair<size_t, bool> AddLine(const Point &a, const Point &b)
        {
            Vec3 pos = a.location3D;
            Vec3 dir = GetDirection(pos, b.location3D);

            auto lineId = AddLine(pos, dir);

            if (!lines[lineId.first].IsPointOnLine(a.location3D))
            {
                if (messages)
                {
                    printf("bad point in AddLine\n");
                }
            }
            if (!lines[lineId.first].IsPointOnLine(b.location3D))
            {
                if (messages)
                {
                    printf("bad point in AddLine\n");
                }
            }

            if (!aabb.contains(a.location3D))
            {
                if (messages)
                {
                    printf("bad points in AddLine\n");
                }
            }
            if (!aabb.contains(b.location3D))
            {
                if (messages)
                {
                    printf("bad points in AddLine\n");
                }
            }

            lines[lineId.first].AddPointToLine(lines[lineId.first].GetPosOnLine(a.location3D), a.id);
            lines[lineId.first].AddPointToLine(lines[lineId.first].GetPosOnLine(b.location3D), b.id);

            return lineId;
        }

        /*
            Add a line defined by a position and direction vector to this plane.
            If an equivalent line already exists, returns its index with false.
            Otherwise creates a new line (normalizing the direction) and returns
            its index with true.
        */
        std::pair<size_t, bool> AddLine(const Vec3 &pos, const Vec3 &dir)
        {
            for (auto &line : lines)
            {
                if (line.IsEqualTo(pos, dir))
                {
                    return {line.id, false};
                }
            }

            Line l;
            l.id = lines.size();
            l.origin = pos;
            //          l.direction = dir;
            Vec3 temp = glm::normalize(dir);
            l.direction = temp;

            lines.push_back(l);

            return {l.id, true};
        }

        void RemoveLastLine()
        {
            lines.pop_back();
        }

        bool IsEqualTo(const Vec3 &n, double d)
        {
            return (equals(normal, n, toleranceVectorEquality) && equals(distance, d, TOLERANCE_SCALAR_EQUALITY));
        }

        glm::dvec2 GetPosOnPlane(const glm::dvec3 &pos)
        {
            return {};
        }

        bool HasOverlap(const std::pair<size_t, size_t> &A, const std::pair<size_t, size_t> &B)
        {
            return (A.first == B.first || A.first == B.second || A.second == B.first || A.second == B.second);
        }

        /*
            Associate point p with every line on this plane that it lies on.
            For each matching line, a ReferenceLine entry is added to the point,
            recording the line id and the point's projected distance along it.
        */
        void PutPointOnLines(Point &p)
        {
            for (auto &l : lines)
            {
                if (l.IsPointOnLine(p.location3D))
                {
                    ReferenceLine ref;
                    ref.pointID = p.id;
                    ref.lineID = l.id;
                    ref.location = l.GetPosOnLine(p.location3D);
                    p.lines.push_back(ref);
                }
            }
        }

        /*
                Normal is assumed to be normalised.
        */
        bool IsPointOnPlane(const glm::dvec3 &pos)
        {
            double d = glm::dot(normal, pos);
            double posLength = pos.length();
            //          return equals(distance, d, toleranceVectorEquality * posLength);
            return equals(distance, d, toleranceVectorEquality);
        }

        /*
            Construct an orthonormal 2D coordinate basis on the plane surface.
            The origin is the closest point on the plane to the world origin.
            The 'left' and 'right' axes span the plane and are perpendicular
            to the plane normal. Used for projecting 3D points onto the plane
            for 2D triangulation.
        */
        PlaneBasis MakeBasis()
        {
            glm::dvec3 origin = normal * distance;
            glm::dvec3 up = normal;

            glm::dvec3 worldUp = glm::dvec3(0, 1, 0);
            glm::dvec3 worldRight = glm::dvec3(1, 0, 0);

            bool normalIsUp = equals(up, worldUp, EPS_SMALL) || equals(-up, worldUp, EPS_SMALL);
            glm::dvec3 left = normalIsUp ? glm::cross(up, worldRight) : glm::cross(up, worldUp);
            glm::dvec3 right = glm::cross(left, up);

            PlaneBasis basis;

            basis.origin = origin;
            basis.up = glm::normalize(up);
            basis.left = glm::normalize(left);
            basis.right = glm::normalize(right);

            return basis;
        }
    };

    struct Triangle
    {
        size_t id;

        size_t a;
        size_t b;
        size_t c;

        void Flip()
        {
            auto temp = a;
            a = b;
            b = temp;
        }

        bool HasPoint(size_t p)
        {
            return a == p || b == p || c == p;
        }

        bool IsNeighbour(Triangle &t)
        {
            return HasPoint(t.a) || HasPoint(t.b) || HasPoint(t.c);
        }

        bool SamePoints(Triangle &other)
        {
            return other.HasPoint(a) && other.HasPoint(b) && other.HasPoint(c);
        }

        size_t GetNotShared(Triangle &other)
        {
            if (!other.HasPoint(a))
            {
                return a;
            }
            else if (!other.HasPoint(b))
            {
                return b;
            }
            else if (!other.HasPoint(c))
            {
                return c;
            }

            // throw std::exception("Missing common point in GetNotShared");
        }
    };

    struct SegmentSet
    {
        std::vector<std::pair<size_t, size_t>> segments;
        std::vector<Triangle> triangles;
        std::map<std::pair<size_t, size_t>, size_t> segmentCounts;
        std::vector<size_t> irrelevantFaces;
        std::vector<size_t> irrelevantFaces_toTest;

        std::map<size_t, std::vector<std::pair<size_t, size_t>>> planeSegments;
        std::map<size_t, std::map<std::pair<size_t, size_t>, size_t>> planeSegmentCounts;

        /*
            Register an edge segment (a, b) associated with the given plane.
            The segment is stored in canonical order (smaller id first) and its
            occurrence count is incremented -- used later to distinguish boundary
            edges (count == 1) from interior edges (count == 2).
        */
        void AddSegment(size_t planeId, size_t a, size_t b)
        {
            if (a == b)
            {
                if (messages)
                {
                    printf("a == b in AddSegment\n");
                }
                return;
            }

            auto seg = a < b ? std::make_pair(a, b) : std::make_pair(b, a);

            segments.emplace_back(seg);
            segmentCounts[seg]++;

            planeSegments[planeId].emplace_back(seg);
            planeSegmentCounts[planeId][seg]++;
        }

        /*
            Register a triangle face (a, b, c) by adding its three edges as
            segments and storing the triangle in the triangles list.
        */
        void AddFace(size_t planeId, size_t a, size_t b, size_t c)
        {
            AddSegment(planeId, a, b);
            AddSegment(planeId, b, c);
            AddSegment(planeId, c, a);

            Triangle t;
            t.id = triangles.size();
            t.a = a;
            t.b = b;
            t.c = c;

            triangles.push_back(t);
        }

        std::vector<size_t> GetTrianglesWithPoint(size_t p)
        {
            std::vector<size_t> returnTriangles;

            for (auto &t : triangles)
            {
                if (t.HasPoint(p))
                {
                    returnTriangles.push_back(t.id);
                }
            }

            return returnTriangles;
        }

        std::vector<size_t> GetTrianglesWithEdge(size_t a, size_t b)
        {
            std::vector<size_t> returnTriangles;

            for (auto &t : triangles)
            {
                if (t.HasPoint(a) && t.HasPoint(b))
                {
                    returnTriangles.push_back(t.id);
                }
            }

            return returnTriangles;
        }

        std::vector<std::pair<std::pair<size_t, size_t>, std::vector<size_t>>> GetNeighbourTriangles(Triangle &triangle)
        {
            std::vector<std::pair<std::pair<size_t, size_t>, std::vector<size_t>>> returnTriangles;

            {
                auto tris = GetTrianglesWithEdge(triangle.a, triangle.b);
                returnTriangles.emplace_back(std::make_pair(triangle.a, triangle.b), tris);
            }
            {
                auto tris = GetTrianglesWithEdge(triangle.b, triangle.c);
                returnTriangles.emplace_back(std::make_pair(triangle.b, triangle.c), tris);
            }
            {
                auto tris = GetTrianglesWithEdge(triangle.c, triangle.a);
                returnTriangles.emplace_back(std::make_pair(triangle.c, triangle.a), tris);
            }

            return returnTriangles;
        }

        /*
            Check whether the mesh is manifold: every edge must be shared by
            exactly 2 triangles. Returns true if no boundary or non-manifold
            edges exist.
        */
        bool IsManifold()
        {
            std::vector<std::pair<size_t, size_t>> contours;

            for (auto &[pair, count] : segmentCounts)
            {
                if (count != 2)
                {
                    contours.push_back(pair);
                }
            }

            return contours.empty();
        }

        /*
            Return the boundary (contour) edges grouped by plane id.
            A contour edge is one that appears exactly once in a given plane's
            segment set -- i.e. it borders only one triangle and thus lies on
            the mesh boundary.
        */
        std::map<size_t, std::vector<std::pair<size_t, size_t>>> GetContourSegments()
        {
            std::map<size_t, std::vector<std::pair<size_t, size_t>>> contours;

            for (auto &[plane, segmentCounts] : planeSegmentCounts)
            {
                for (auto &[pair, count] : segmentCounts)
                {
                    if (count == 1)
                    {
                        contours[plane].push_back(pair);
                    }
                }
            }

            return contours;
        }
    };

    struct SharedPosition
    {
        std::vector<Point> points;
        std::vector<Plane> planes;

        SegmentSet A;
        SegmentSet B;

        // TODO: design flaw
        const Geometry *_linkedA;
        const Geometry *_linkedB;

        Geometry relevantA;
        Geometry relevantB;

        BVH relevantBVHA;
        BVH relevantBVHB;

        std::unordered_map<std::tuple<int64_t, int64_t, int64_t>, std::vector<size_t>> pointGrid;

        /*
            Given a base triangle and a set of candidate triangle ids (all sharing
            an edge with the base and correctly oriented), select the triangle whose
            surface normal is most "upward" relative to the base. This is used during
            winding-order propagation to pick the best neighbour to orient next.
            Returns the id of the selected triangle.
        */
        size_t FindUppermostTriangleId(Triangle &base, const std::vector<size_t> &triangleIds)
        {
            if (triangleIds.size() == 1)
            {
                return triangleIds[0];
            }

            auto baseNorm = GetNormal(base);

            size_t maxDotId = -1;
            double maxDot = -DBL_MAX;

            for (auto id : triangleIds)
            {
                if (id == base.id)
                    continue;

                auto triNorm = GetNormal(A.triangles[id]);

                auto other = A.triangles[id].GetNotShared(base);

                auto above = CalcTriPt(base, other);

                auto dot = (glm::dot(baseNorm, triNorm) + 1.0) / 2.0; // range dot from 0 to 1, positive for above, negative for below

                if (above == TriangleVsPoint::BELOW)
                {
                    dot = dot - 1.0;
                }
                else
                {
                    dot = 1.0 - dot;
                }

                if (dot > maxDot)
                {
                    maxDot = dot;
                    maxDotId = id;
                }
            }

            return maxDotId;
        }

        size_t GetPointWithMaxY()
        {
            double max = -DBL_MAX;
            size_t pointID = 0;

            for (auto &p : points)
            {
                if (p.location3D.y > max)
                {
                    max = p.location3D.y;
                    pointID = p.id;
                }
            }

            return pointID;
        }

        enum class TriangleVsPoint
        {
            ABOVE,
            BELOW,
            ON
        };

        /*
            Classify the position of a point relative to a triangle's plane.
            Returns ABOVE if the point is on the positive-normal side, BELOW if
            on the negative side, or ON if the point lies (within tolerance)
            on the triangle's plane.
        */
        TriangleVsPoint CalcTriPt(Triangle &T, size_t point)
        {
            auto norm = GetNormal(T);

            auto dpt3d = points[point].location3D - points[T.a].location3D;
            auto dot = glm::dot(norm, dpt3d);

            if (std::fabs(dot) < EPS_BIG)
                return TriangleVsPoint::ON;
            if (dot > 0.0)
                return TriangleVsPoint::ABOVE;
            return TriangleVsPoint::BELOW;
        }

        /*
            Determine whether a neighbouring triangle's winding order must be
            flipped to be consistent with triangle T. Uses the relative positions
            of the non-shared vertices: if E (neighbour's unique vertex) is above T,
            then A (T's unique vertex) should be above the neighbour, and vice versa.
            Returns true if the neighbour needs to be flipped.
        */
        bool ShouldFlip(Triangle &T, Triangle &neighbour)
        {
            /*
                        for each n in N, orient n by formula:
                        if triangle T and N are vertices ABCDEF with BCDE as the shared edge, consider three cases:

                        E is above T, then A is above n
                        if dot(normal(T), E) > 0 => dot(normal(n), A) > 0, else flip n
                        E is below T, then A is below n
                        if dot(normal(T), E) < 0 => dot(normal(n), A) < 0, else flip n
                        E is on T, then normal of n and T are equal
                        if dot(normal(T), E) == 0 => dot(normal(n), normal(T)) == 1, else flip n
            */
            auto normT = GetNormal(T);
            auto normNB = GetNormal(neighbour);

            if (T.SamePoints(neighbour))
            {
                // same tri, different winding, flip if not the same normal
                return glm::dot(normT, normNB) < 1.0 - EPS_BIG;
            }

            auto A = T.GetNotShared(neighbour);
            auto E = neighbour.GetNotShared(T);

            TriangleVsPoint EvsT = CalcTriPt(T, E);
            TriangleVsPoint AvsN = CalcTriPt(neighbour, E);

            if (EvsT == TriangleVsPoint::ABOVE && AvsN != TriangleVsPoint::ABOVE)
            {
                return true;
            }
            else if (EvsT == TriangleVsPoint::BELOW && AvsN != TriangleVsPoint::BELOW)
            {
                return true;
            }
            else if (EvsT == TriangleVsPoint::ON && AvsN != TriangleVsPoint::ON)
            {
                return true;
            }
            else
            {
                // neighbour already good!
                return false;
            }
        }

        /*
            Compute the unit normal of a triangle using its three vertex positions.
            Falls back to a default direction if the triangle is degenerate.
        */
        Vec3 GetNormal(Triangle &tri)
        {
            Vec3 temp(-1.0, -1.0, -1.0);
            Vec3 norm = glm::normalize(temp);
            computeSafeNormal(points[tri.a].location3D, points[tri.b].location3D, points[tri.c].location3D, norm, EPS_SMALL);
            return norm;
        }

        SegmentSet &GetSegSetA()
        {
            return A;
        }

        /*
            Add a 3D point to the shared point set with spatial-hash deduplication.
            A uniform grid (cell size = toleranceVectorEquality) is used so that
            only the 27 neighbouring cells need to be checked for existing matches.
            If an existing point within tolerance is found, its id is returned;
            otherwise a new point is created and inserted into the grid.
        */
        size_t AddPoint(const Vec3& newPoint)
        {
            // 1. Compute the grid cell for the query point
            const double cellSize = toleranceVectorEquality;   // same tolerance you already use for ==

            auto getKey = [&](const Vec3& p) -> std::tuple<int64_t, int64_t, int64_t> {
                return {
                    static_cast<int64_t>(std::floor(p.x / cellSize)),
                    static_cast<int64_t>(std::floor(p.y / cellSize)),
                    static_cast<int64_t>(std::floor(p.z / cellSize))
                };
                };

            const auto centerKey = getKey(newPoint);

            // 2. Check the point against all 27 neighbouring cells (guaranteed to contain any point
            //    that is within toleranceVectorEquality because |delta_x|,|delta_y|,|delta_z| < tolerance)
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        const std::tuple<int64_t, int64_t, int64_t> neighbourKey = {
                            std::get<0>(centerKey) + dx,
                            std::get<1>(centerKey) + dy,
                            std::get<2>(centerKey) + dz
                        };

                        auto it = pointGrid.find(neighbourKey);
                        if (it != pointGrid.end()) {
                            for (size_t existingId : it->second) {
                                if (points[existingId] == newPoint) {   // re-uses your existing tolerance check
                                    return existingId;
                                }
                            }
                        }
                    }
                }
            }

            // 3. Point does not exist -> insert it
            Point p;
            p.id = points.size();
            p.location3D = newPoint;
            points.push_back(p);

            // Insert into spatial grid (same cell as the point itself)
            pointGrid[getKey(newPoint)].push_back(p.id);

            return p.id;
        }

        /*
            Add a plane (defined by its normal and signed distance from the origin)
            to the shared plane set. If a plane with the same reference id or
            equivalent (normal, distance) already exists, its id is returned.
            Otherwise a new plane is created with the normal normalized.
        */
        size_t AddPlane(const Vec3 &normal, double d, uint32_t refId)
        {
            for (auto &plane : planes)
            {
                if (plane.refPlane == refId || plane.IsEqualTo(normal, d))
                {
                    return plane.id;
                }
            }

            Plane p;
            p.id = planes.size();
            p.refPlane = refId;
            p.normal = glm::normalize(normal);
            p.distance = d;
            planes.push_back(p);

            return p.id;
        }

        /*
            Main entry point: populate the SharedPosition from two input geometries.
            Each geometry's faces are filtered by AABB overlap with the other and
            their planes, points, and segments are inserted into the shared set.
            isUnion controls how irrelevant (non-overlapping) faces are handled later.
        */
        void Construct(const Geometry &A, const Geometry &B, bool isUnion, const BooleanBudget& budget)
        {
            budget.CheckDeadline("SharedPosition Construct start");

            auto boxA = A.GetAABB();
            auto boxB = B.GetAABB();

            AddGeometry(A, B, boxB, true, isUnion, 0, budget);
            AddGeometry(B, A, boxA, false, isUnion, A.planes.size(), budget);

            _linkedA = &A;
            _linkedB = &B;

            budget.CheckDeadline("SharedPosition Construct complete");
        }

        /*
            Process one input geometry's faces into the shared position.
            For each face, check AABB overlap with the other geometry's bounding
            box and individual face boxes. Faces that don't overlap are classified
            as "irrelevant" (fully outside) or "irrelevant_toTest" (no face-face
            contact but within the other's bounds). Overlapping faces have their
            planes, points, and segments added to the shared set, and a BVH is
            built over the relevant faces for later inside/outside queries.
        */
        void AddGeometry(const Geometry &geom, const Geometry &secondGeom, const AABB &relevantBounds, bool isA, bool isUnion, uint32_t offsetPlane, const BooleanBudget& budget)
        {
#ifdef CSG_DEBUG_OUTPUT
            Geometry relevant;
#endif

            uint64_t iteration = 0;
            for (size_t i = 0; i < geom.numFaces; i++)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline(isA ? "SharedPosition AddGeometry A" : "SharedPosition AddGeometry B");
                }

                Face f = geom.GetFace(i);

                auto faceBox = geom.GetFaceBox(i);

                if (!faceBox.intersects(relevantBounds))
                {
                    if (isA)
                    {
                        A.irrelevantFaces.push_back(i);
                    }
                    else
                    {
                        B.irrelevantFaces.push_back(i);
                    }

                    continue;
                }

                bool contact = false;

                for (size_t j = 0; j < secondGeom.numFaces; j++)
                {
                    if ((iteration++ & 1023ULL) == 0)
                    {
                        budget.CheckDeadline(isA ? "SharedPosition AddGeometry A overlap" : "SharedPosition AddGeometry B overlap");
                    }

                    auto faceBox2 = secondGeom.GetFaceBox(j);

                    if (faceBox.intersects(faceBox2))
                    {
                        contact = true;
                        break;
                    }
                }

                if (!contact)
                {
                    if (isA)
                    {
                        A.irrelevantFaces_toTest.push_back(i);
                    }
                    else
                    {
                        if (isUnion || geom.numFaces < 2000) // TODO: This condition is wrong but efficient for large models
                        {
                            B.irrelevantFaces_toTest.push_back(i);
                        }
                    }

                    continue;
                }

                if (isA)
                {
#ifdef CSG_DEBUG_OUTPUT
//                    DumpGeometry(geom, L"Initial_A.obj");
#endif
                }
                else
                {
#ifdef CSG_DEBUG_OUTPUT
//                    DumpGeometry(geom, L"Initial_B.obj");
#endif
                }

                auto a = geom.GetPoint(f.i0);
                auto b = geom.GetPoint(f.i1);
                auto c = geom.GetPoint(f.i2);

#ifdef CSG_DEBUG_OUTPUT
//                relevant.AddFace(a, b, c, -1);
#endif

                Vec3 norm;
                if (computeSafeNormal(a, b, c, norm, EPS_SMALL))
                {
                    double rs = glm::dot(geom.planes[f.pId].normal, norm);

                    size_t planeId = -1;

                    if (rs < 0)
                    {
                        planeId = AddPlane(-geom.planes[f.pId].normal, -geom.planes[f.pId].distance, f.pId + offsetPlane);
                    }
                    else
                    {
                        planeId = AddPlane(geom.planes[f.pId].normal, geom.planes[f.pId].distance, f.pId + offsetPlane);
                    }

                    auto ia = AddPoint(a);
                    auto ib = AddPoint(b);
                    auto ic = AddPoint(c);

                    double da = glm::dot(norm, a);
                    double db = glm::dot(norm, b);
                    double dc = glm::dot(norm, c);

                    if (!planes[planeId].IsPointOnPlane(a))
                    {
                        if (messages)
                        {
                            printf("unexpected point on plane in AddGeometry\n");
                        }
                        if (messages)
                        {
                            printf("a = (%12.8f, %12.8f, %12.8f), da = %12.8f, distance = %12.8f\n", a.x, a.y, a.z, da, planes[planeId].distance);
                        }
                    }
                    if (!planes[planeId].IsPointOnPlane(b))
                    {
                        if (messages)
                        {
                            printf("unexpected point on plane in AddGeometry\n");
                        }
                        if (messages)
                        {
                            printf("b = (%12.8f, %12.8f, %12.8f), db = %12.8f, distance = %12.8f\n", b.x, b.y, b.z, db, planes[planeId].distance);
                        }
                    }
                    if (!planes[planeId].IsPointOnPlane(c))
                    {
                        if (messages)
                        {
                            printf("unexpected point on plane in AddGeometry\n");
                        }
                        if (messages)
                        {
                            printf("c = (%12.8f, %12.8f, %12.8f), dc = %12.8f, distance = %12.8f\n", c.x, c.y, c.z, dc, planes[planeId].distance);
                        }
                    }
                    planes[planeId].AddPoint(a);
                    planes[planeId].AddPoint(b);
                    planes[planeId].AddPoint(c);

                    if (isA)
                    {
                        A.AddFace(planeId, ia, ib, ic);
                        relevantA.AddFace(a, b, c, planes[planeId].refPlane);
                    }
                    else
                    {
                        B.AddFace(planeId, ia, ib, ic);
                        relevantB.AddFace(a, b, c, planes[planeId].refPlane);
                    }
                }
                else
                {
                    if (messages)
                    {
                        printf("Degenerate face in AddGeometry\n");
                    }
                }
            }

            if (isA)
            {
                relevantBVHA = MakeBVH(relevantA);
            }
            else
            {
                relevantBVHB = MakeBVH(relevantB);
            }

#ifdef CSG_DEBUG_OUTPUT
            // if (isA)
            // {
            //     DumpGeometry(relevant, L"relevantA.obj");
            // }
            // else
            // {
            //     DumpGeometry(relevant, L"relevantB.obj");
            // }
#endif
        }

        std::vector<size_t> GetPointsOnPlane(Plane &p)
        {
            auto cp = planeToPoints[p.id];
            std::sort(cp.begin(), cp.end());
            cp.erase(std::unique(cp.begin(), cp.end()), cp.end());
            return cp;
        }

        /*
            Given two sorted lists of intersection distances along a shared line
            (one per plane), compute the overlapping segments. The overlap region
            is [max(a_first, b_first), min(a_last, b_last)]. All distances from
            both lists that fall in this range are merged and sorted, then
            consecutive pairs form the output segments.
            Returns a list of (start_distance, end_distance) pairs.
        */
        std::vector<std::pair<double, double>> BuildSegments(const std::vector<double> &a, const std::vector<double> &b) const
        {
            if (a.size() == 0 || b.size() == 0)
            {
                return {};
            }

            // we need to figure out the overlap between the two lists of intersections
            // we can be clever here and try to conclude that the first point must be the start of a segment and the next point would end that segment
            // however, we would really shoot ourselves in the foot as any coplanar results are missing from these two sets
            // let's just make some segments that span both intersection lists, and eat the overhead

            double min = std::max(a[0], b[0]);
            double max = std::min(a[a.size() - 1], b[b.size() - 1]);

            std::vector<double> points;

            for (size_t i = 0; i < a.size(); i++)
            {
                double val = a[i];
                if (val >= min && val <= max)
                {
                    points.push_back(val);
                }
            }

            for (size_t i = 0; i < b.size(); i++)
            {
                double val = b[i];
                if (val >= min && val <= max)
                {
                    points.push_back(val);
                }
            }

            const auto double_less = +[](double left, double right)
            { return left < right; };
            std::sort(points.begin(), points.end(), double_less);

            // Epsilon-dedup: two intersection distances that differ only by
            // floating-point noise must collapse to a single point, otherwise
            // the pair-generation below emits a zero-or-near-zero-length
            // segment that later spawns a micro-triangle cascade in
            // Normalize (the main engine of test61's face blow-up).
            std::vector<double> uniq;
            uniq.reserve(points.size());
            for (double v : points)
            {
                if (uniq.empty() || std::abs(v - uniq.back()) > _TOLERANCE_PLANE_INTERSECTION)
                {
                    uniq.push_back(v);
                }
            }

            std::vector<std::pair<double, double>> result;
            result.reserve(uniq.size());

            for (size_t i = 1; i < uniq.size(); i++)
            {
                // Skip any pair whose two endpoints are within tolerance of
                // each other. Belt-and-braces after the dedup above.
                if (std::abs(uniq[i] - uniq[i - 1]) <= _TOLERANCE_PLANE_INTERSECTION) continue;
                result.emplace_back(uniq[i - 1], uniq[i]);
            }

            return result;
        }

        /*
            Merge a line's segments into a sorted, non-overlapping chain of
            point-id pairs. All segment endpoints are collected, sorted by their
            distance along the line, and then consecutive distinct points form
            the output segments. This fills gaps between segments -- those gaps
            are removed later during inside/outside classification.
        */
        std::vector<std::pair<size_t, size_t>> GetNonIntersectingSegments(Line &l)
        {
            std::vector<std::pair<size_t, double>> pointsInOrder;

            for (auto segment : l.GetSegments())
            {
                if (!l.IsPointOnLine(points[segment.first].location3D))
                {
                    if (messages)
                    {
                        printf("point not on line in GetNonIntersectingSegments\n");
                    }
                    if (messages)
                    {
                        printf("points[segment.first].location3D = (%12.8f, %12.8f, %12.8f)\n", points[segment.first].location3D.x, points[segment.first].location3D.y, points[segment.first].location3D.z);
                    }
                    if (messages)
                    {
                        printf("l.origin                         = (%12.8f, %12.8f, %12.8f)\n", l.origin.x, l.origin.y, l.origin.z);
                    }
                    if (messages)
                    {
                        printf("l.direction                      = (%12.8f, %12.8f, %12.8f)\n", l.direction.x, l.direction.y, l.direction.z);
                    }
                }

                if (!l.IsPointOnLine(points[segment.second].location3D))
                {
                    if (messages)
                    {
                        printf("point not on line in GetNonIntersectingSegments\n");
                    }
                    if (messages)
                    {
                        printf("points[segment.second].location3D = (%12.8f, %12.8f, %12.8f)\n", points[segment.second].location3D.x, points[segment.second].location3D.y, points[segment.second].location3D.z);
                    }
                    if (messages)
                    {
                        printf("l.origin                          = (%12.8f, %12.8f, %12.8f)\n", l.origin.x, l.origin.y, l.origin.z);
                    }
                    if (messages)
                    {
                        printf("l.direction                       = (%12.8f, %12.8f, %12.8f)\n", l.direction.x, l.direction.y, l.direction.z);
                    }
                }

                pointsInOrder.emplace_back(segment.first, l.GetPosOnLine(points[segment.first].location3D));
                pointsInOrder.emplace_back(segment.second, l.GetPosOnLine(points[segment.second].location3D));
            }

            std::sort(
                pointsInOrder.begin(), pointsInOrder.end(), [&](const std::pair<size_t, double> &left, const std::pair<size_t, double> &right)
                { return left.second > right.second; });

            std::vector<std::pair<size_t, size_t>> segmentsWithoutIntersections;

            if (pointsInOrder.empty())
                return {};

            size_t cur = pointsInOrder[0].first;
            for (size_t i = 1; i < pointsInOrder.size(); i++)
            {
                // this fills gaps, but we don't care about that in this stage
                // gaps filled will be removed during inside/outside checking
                size_t next = pointsInOrder[i].first;
                if (cur != next)
                {
                    segmentsWithoutIntersections.emplace_back(cur, next);
                    cur = next;
                }
            }
            return segmentsWithoutIntersections;
        }

        /*
            Re-triangulate a plane using Constrained Delaunay Triangulation (CDT).
            All points on the plane are projected onto its 2D basis, and line
            segments become constrained edges. After triangulation, each candidate
            triangle is tested against both meshes' BVHs (isInsideMesh) and the
            projected boundary (isInsideBoundary) to determine whether it belongs
            to the intersection region. Only triangles that lie on the boundary
            of at least one mesh are added to the output geometry.
        */
        void TriangulatePlane(Geometry &geom, Plane &p, const BooleanBudget& budget)
        {
            budget.CheckDeadline("Normalize TriangulatePlane start");

            // grab all points on the plane
            auto pointsOnPlane = GetPointsOnPlane(p);

            // temporarily project all points for triangulation
            // NOTE: these points should not be used as output, only as placeholder for triangulation!

            auto basis = p.MakeBasis();

            std::unordered_map<size_t, size_t> pointToProjectedPoint;
            pointToProjectedPoint.reserve(pointsOnPlane.size());

            std::unordered_map<size_t, size_t> projectedPointToPoint;
            projectedPointToPoint.reserve(pointsOnPlane.size());

            std::vector<glm::dvec2> projectedPoints;
            projectedPoints.reserve(pointsOnPlane.size());

            uint64_t iteration = 0;
            for (auto &pointId : pointsOnPlane)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize TriangulatePlane project points");
                }

                pointToProjectedPoint[pointId] = projectedPoints.size();
                projectedPointToPoint[projectedPoints.size()] = pointId;
                projectedPoints.push_back(basis.project(points[pointId].location3D));
            }

            std::set<std::pair<size_t, size_t>> edges;
            std::set<std::pair<size_t, size_t>> defaultEdges;
            static int i = 0;
            i++;

            for (auto &line : p.lines)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize TriangulatePlane lines");
                }

                // these segments might intersect internally, lets resolve that so we get a valid chain
                auto segments = GetNonIntersectingSegments(line);

                for (auto &segment : segments)
                {
                    if ((iteration++ & 1023ULL) == 0)
                    {
                        budget.CheckDeadline("Normalize TriangulatePlane segments");
                    }

                    if (pointToProjectedPoint.count(segment.first) == 0)
                    {
                        bool expectedOnPlane = p.IsPointOnPlane(points[segment.first].location3D);
                        if (messages)
                        {
                            printf("unknown point in list, repairing in TriangulateLine\n");
                        }

                        pointToProjectedPoint[segment.first] = projectedPoints.size();
                        projectedPointToPoint[projectedPoints.size()] = segment.first;
                        projectedPoints.push_back(basis.project(points[segment.first].location3D));
                    }
                    if (pointToProjectedPoint.count(segment.second) == 0)
                    {
                        bool expectedOnPlane = p.IsPointOnPlane(points[segment.second].location3D);
                        if (messages)
                        {
                            printf("unknown point in list, repairing in TriangulateLine\n");
                        }

                        pointToProjectedPoint[segment.second] = projectedPoints.size();
                        projectedPointToPoint[projectedPoints.size()] = segment.second;
                        projectedPoints.push_back(basis.project(points[segment.second].location3D));
                    }

                    auto projectedIndexA = pointToProjectedPoint[segment.first];
                    auto projectedIndexB = pointToProjectedPoint[segment.second];

                    if (projectedIndexA != projectedIndexB)
                    {
                        defaultEdges.insert(segment);
                        edges.insert(std::make_pair(projectedIndexA, projectedIndexB));
                    }
                }
            }

#ifdef CSG_DEBUG_OUTPUT
            // std::vector<std::vector<glm::dvec2>> edgesPrinted;

            // for (auto &e : edges)
            // {
            //     edgesPrinted.push_back({projectedPoints[e.first], projectedPoints[e.second]});
            // }

            // DumpSVGLines(edgesPrinted, L"poly.html");
#endif

            CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::AsProvided);
            std::vector<CDT::Edge> cdt_edges;
            cdt_edges.reserve(edges.size());

            std::vector<CDT::V2d<double>> cdt_verts;
            cdt_verts.reserve(projectedPoints.size());

            for (auto &point : projectedPoints)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize TriangulatePlane CDT vertices");
                }

                cdt_verts.emplace_back(CDT::V2d<double>::make(point.x, point.y));
            }

            for (auto &edge : edges)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize TriangulatePlane CDT edges");
                }

                cdt_edges.emplace_back((uint32_t)edge.first, (uint32_t)edge.second);
            }

            budget.CheckDeadline("Normalize TriangulatePlane CDT remap");
            auto mapping = CDT::RemoveDuplicatesAndRemapEdges(cdt_verts, cdt_edges).mapping;
            // mapping[] is OLD->NEW (CDT.h line 64). After the call cdt_verts
            // is compacted and cdt_edges is already remapped to new indices,
            // so CDT works in the NEW index space and tri.vertices[i] is a
            // NEW index. The surrounding code needs OLD indices to look into
            // projectedPoints / projectedPointToPoint (those arrays were
            // built before the dedup), so we build the inverse mapping once
            // here. When several old vertices collapse to the same new one
            // (duplicates), any of their old indices represents the same
            // location; we keep the first.
            std::vector<size_t> newToOld(cdt_verts.size(), SIZE_MAX);
            for (size_t oldIdx = 0; oldIdx < mapping.size(); ++oldIdx)
            {
                size_t newIdx = mapping[oldIdx];
                if (newIdx < newToOld.size() && newToOld[newIdx] == SIZE_MAX)
                {
                    newToOld[newIdx] = oldIdx;
                }
            }

            // CDT's spatial locate infinite-loops on non-finite coordinates, and being a
            // synchronous third-party call no budget deadline check can break it. Non-finite
            // projected coordinates come from a degenerate plane: a zero-length / non-finite
            // normal makes MakeBasis produce a NaN basis, and a non-finite plane distance makes
            // the projection origin NaN -- either way every projected point becomes NaN. Such a
            // plane has no real area and yields no valid triangles, so skip its triangulation
            // rather than hang the kernel. Observed on near-coplanar hollow-core-slab cuts
            // (DLGKF2-WNCSE-Bovenbouw.ifc, element 433835).
            for (const auto& v : cdt_verts)
            {
                if (!std::isfinite(v.x) || !std::isfinite(v.y))
                {
                    if (getenv("CXDIAG_CSG"))
                    {
                        // Attribute the non-finite coordinates to their source: a degenerate plane
                        // (basis/origin NaN poisons every projection) or non-finite 3D points.
                        size_t badPoints = 0;
                        for (auto& pointId : pointsOnPlane)
                        {
                            const glm::dvec3& loc = points[pointId].location3D;
                            if (!std::isfinite(loc.x) || !std::isfinite(loc.y) || !std::isfinite(loc.z)) badPoints++;
                        }
                        std::cerr << "DIAG: TP NaN-skip plane=" << p.id << " refPlane=" << p.refPlane
                                  << " n=(" << p.normal.x << "," << p.normal.y << "," << p.normal.z << ") d=" << p.distance
                                  << " pointsOnPlane=" << pointsOnPlane.size() << " nonFinite3DPoints=" << badPoints
                                  << " verts=" << cdt_verts.size() << std::endl;
                    }
                    return;
                }
            }

            cdt.insertVertices(cdt_verts);
            budget.CheckDeadline("Normalize TriangulatePlane CDT insert vertices");
            cdt.insertEdges(cdt_edges);
            budget.CheckDeadline("Normalize TriangulatePlane CDT insert edges");

            cdt.eraseSuperTriangle();
            budget.CheckDeadline("Normalize TriangulatePlane CDT erase super triangle");

            auto triangles = cdt.triangles;

            // auto contourLoop = FindLargestEdgeLoop(projectedPoints, edges);

#ifdef CSG_DEBUG_OUTPUT
            // std::vector<std::vector<glm::dvec2>> edges3DTriangles;
            // std::set<std::pair<size_t, size_t>> edgesTriangles;
            // std::set<std::pair<size_t, size_t>> finalEdgesTriangles;
#endif

            for (auto &tri : triangles)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize TriangulatePlane triangles");
                    budget.CheckFaceCount(geom.numFaces, "Normalize TriangulatePlane triangles");
                }

#ifdef CSG_DEBUG_OUTPUT
                // edgesTriangles.insert(std::make_pair(tri.vertices[0], tri.vertices[1]));
                // edgesTriangles.insert(std::make_pair(tri.vertices[1], tri.vertices[2]));
                // edgesTriangles.insert(std::make_pair(tri.vertices[0], tri.vertices[2]));
#endif

                size_t pointIdA = projectedPointToPoint[newToOld[tri.vertices[0]]];
                size_t pointIdB = projectedPointToPoint[newToOld[tri.vertices[1]]];
                size_t pointIdC = projectedPointToPoint[newToOld[tri.vertices[2]]];

                auto ptA = points[pointIdA].location3D;
                auto ptB = points[pointIdB].location3D;
                auto ptC = points[pointIdC].location3D;

                glm::dvec3 v1 = glm::normalize(ptA - ptB);
                glm::dvec3 v2 = glm::normalize(ptA - ptC);
                glm::dvec3 v3 = glm::normalize(ptB - ptC);
                double rs1 = glm::dot(v1, v2);
                double rs2 = glm::dot(v2, v3);
                double rs3 = glm::dot(v1, v3);

                if (std::abs(rs1) > 1 - toleranceThinTriangle ||
                    std::abs(rs2) > 1 - toleranceThinTriangle ||
                    std::abs(rs3) > 1 - toleranceThinTriangle)
                {
                    // Thin CDT sliver. Rejecting it outright abandons surface area the CDT
                    // assigned to this plane, leaving its neighbours' edges unpaired
                    // (test53d: open-edge fans + double-cover along the wall bottom edge).
                    // The ray gates below are unreliable for slivers (their winding normal is
                    // numerically unstable), so decide purely in 2D: emit the sliver when it
                    // lies inside the face's boundary polygon, otherwise it is fill outside
                    // the real contour (super-triangle artifacts) and stays rejected.
                    static const bool noThinKeep = std::getenv("CX_CSG_NO_THINKEEP") != nullptr;
                    glm::dvec2 s1 = projectedPoints[newToOld[tri.vertices[0]]];
                    glm::dvec2 s2 = projectedPoints[newToOld[tri.vertices[1]]];
                    glm::dvec2 s3 = projectedPoints[newToOld[tri.vertices[2]]];
                    if (!noThinKeep && isInsideBoundary(s1, s2, s3, edges, projectedPoints))
                    {
#if 1 // CLIP_DIAG
                        if (messages){
                            if (std::fabs(p.normal.z) > 0.99) printf("[TPDIAG] thin-keep z=%.1f c=(%.1f,%.1f)\n", ptA.z, ((ptA+ptB+ptC)/3.0).x, ((ptA+ptB+ptC)/3.0).y);
                        }
#endif
                        geom.AddFace(ptB, ptA, ptC, p.refPlane);
                        budget.CheckFaceCount(geom.numFaces, "Normalize TriangulatePlane output");
                    }
#if 1 // CLIP_DIAG
                    else if (messages){
                        if (std::fabs(p.normal.z) > 0.99) printf("[TPDIAG] thin-reject z=%.1f c=(%.1f,%.1f)\n", ptA.z, ((ptA+ptB+ptC)/3.0).x, ((ptA+ptB+ptC)/3.0).y);
                    }
#endif
                    continue;
                }

                auto pt2DA = projectedPoints[newToOld[tri.vertices[0]]];
                auto pt2DB = projectedPoints[newToOld[tri.vertices[1]]];
                auto pt2DC = projectedPoints[newToOld[tri.vertices[2]]];

                auto triCenter = (ptA + ptB + ptC) / 3.0;

                Vec raydir = computeNormal(ptA, ptB, ptC);

                auto posA = isInsideMesh(triCenter, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                auto posB = isInsideMesh(triCenter, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                if (posA.loc != MeshLocation::BOUNDARY && posB.loc != MeshLocation::BOUNDARY)
                {
#if 1 // CLIP_DIAG
                    if (messages) {
                        if (std::fabs(p.normal.z) > 0.99) {
                            printf("[TPDIAG] gate1-reject z=%.1f c=(%.1f,%.1f) posA=%d posB=%d\n", triCenter.z, triCenter.x, triCenter.y, (int)posA.loc, (int)posB.loc);
                        }
                    }
#endif
                    continue;
                }

                // If the 2D triangle is not inside the boundaries of the projected boundary of the face it requires further verification
                // It can't be discarded because inside/outside could fail when boundaries have internal partitions
                // Therefore new tests are required to verify that the triangle is on the boundary of A or B

                glm::dvec2 t1 = projectedPoints[newToOld[tri.vertices[0]]];
                glm::dvec2 t2 = projectedPoints[newToOld[tri.vertices[1]]];
                glm::dvec2 t3 = projectedPoints[newToOld[tri.vertices[2]]];

                bool inside2d = isInsideBoundary(t1, t2, t3, edges, projectedPoints);

                if (!inside2d)
                {

                    auto postA = isInsideMesh(triCenter, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                    auto postB = isInsideMesh(triCenter, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                    if (postA.loc != MeshLocation::BOUNDARY && postB.loc != MeshLocation::BOUNDARY)
                    {
                        continue;
                    }

                    auto ptt = glm::mix(triCenter, ptA, triangleEvaluationFactor);

                    postA = isInsideMesh(ptt, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                    postB = isInsideMesh(ptt, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                    if (postA.loc != MeshLocation::BOUNDARY && postB.loc != MeshLocation::BOUNDARY)
                    {
                        continue;
                    }

                    ptt = glm::mix(triCenter, ptB, triangleEvaluationFactor);

                    postA = isInsideMesh(ptt, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                    postB = isInsideMesh(ptt, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                    if (postA.loc != MeshLocation::BOUNDARY && postB.loc != MeshLocation::BOUNDARY)
                    {
                        continue;
                    }

                    ptt = glm::mix(triCenter, ptC, triangleEvaluationFactor);

                    postA = isInsideMesh(ptt, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                    postB = isInsideMesh(ptt, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                    if (postA.loc != MeshLocation::BOUNDARY && postB.loc != MeshLocation::BOUNDARY)
                    {
#if 1 // CLIP_DIAG
                        if (messages) {
                            if (std::fabs(p.normal.z) > 0.99) printf("[TPDIAG] gate2-reject z=%.1f c=(%.1f,%.1f)\n", triCenter.z, triCenter.x, triCenter.y);
                        }
#endif
                        continue;
                    }
                }

                // although CDT is great, it spits out too many or too little tris, we fix it manually
                // if (!IsPointInsideLoop(projectedPoints, contourLoop, triCenter))
                //{
                //    printf("removing point outside loop\n");
                //    continue;
                //}

                // Emission winding: the (ptB, ptA, ptC) swap compensates the uniformly left-handed
                // plane basis (right = cross(left, up) makes left x right = -up), so emitted
                // fragments follow the shared plane's normal.
                // Session 8 measured-and-reverted: orienting each fragment by the boundary-probe
                // normal (posA/posB) instead fixed the coplanar backfacing-strip class at its
                // source (test61 slm 100 -> 39, suite nm 10 -> 2) but the single-probe normal is
                // too noisy as an oracle -- 8 other files regressed (test51 slm 11 -> 29, test62
                // opens 0 -> 79, suite pass 25 -> 23). A reliable per-fragment surface reference
                // (e.g. the operand face id from the BVH boundary hit instead of one jittered ray)
                // is the exact next lever for this class.
                geom.AddFace(ptB, ptA, ptC, p.refPlane);
                budget.CheckFaceCount(geom.numFaces, "Normalize TriangulatePlane output");

#ifdef CSG_DEBUG_OUTPUT
                // edges3DTriangles.push_back({ glm::dvec2(ptA.z+ ptA.x/2, ptA.y+ ptA.x/2), glm::dvec2(ptB.z+ ptB.x/2, ptB.y+ ptB.x/2) });
                // edges3DTriangles.push_back({ glm::dvec2(ptA.z+ ptA.x/2, ptA.y+ ptA.x/2), glm::dvec2(ptC.z+ ptC.x/2, ptC.y+ ptC.x/2) });
                // edges3DTriangles.push_back({ glm::dvec2(ptB.z+ ptB.x/2, ptB.y+ ptB.x/2), glm::dvec2(ptC.z+ ptC.x/2, ptC.y+ ptC.x/2) });
                // DumpSVGLines(edges3DTriangles, L"edges_tri.html");
#endif

#ifdef CSG_DEBUG_OUTPUT
                // finalEdgesTriangles.insert(std::make_pair(tri.vertices[0], tri.vertices[1]));
                // finalEdgesTriangles.insert(std::make_pair(tri.vertices[1], tri.vertices[2]));
                // finalEdgesTriangles.insert(std::make_pair(tri.vertices[0], tri.vertices[2]));
#endif
            }

#ifdef CSG_DEBUG_OUTPUT
            // std::vector<std::vector<glm::dvec2>> edgesPrinted2;

            // for (auto& e : edgesTriangles)
            // {
            //     edgesPrinted2.push_back({ projectedPoints[e.first], projectedPoints[e.second] });
            // }

            // DumpSVGLines(edgesPrinted2, L"poly_triangulation.html");

            // std::vector<std::vector<glm::dvec2>> finalEdgesPrinted;

            // for (auto& e : finalEdgesTriangles)
            // {
            //     finalEdgesPrinted.push_back({ projectedPoints[e.first], projectedPoints[e.second] });
            // }

            // DumpSVGLines(finalEdgesPrinted, L"final_poly_triangulation.html");
#endif
        }
        std::unordered_map<size_t, std::vector<size_t>> planeToLines;
        std::unordered_map<size_t, std::vector<size_t>> planeToPoints;

        void AddRefPlaneToPoint(size_t point, size_t plane)
        {
            ReferencePlane ref;
            ref.pointID = point;
            ref.planeID = plane;
            points[point].planes.push_back(ref);
            planeToPoints[plane].push_back(point);
        }
    };

    /*
        Insert intersection segments into a plane's line structure.
        For each (start, end) distance pair on templine, the corresponding 3D
        points are added to the shared position and registered on the plane's
        matching (or newly created) line. Note: templine and the plane's line
        may differ because AddLine can return an equivalent but non-identical line.
    */
    inline void AddSegments(Plane &p, SharedPosition &sp, Line &templine, const std::vector<std::pair<double, double>> &segments, const BooleanBudget& budget)
    {
        budget.CheckDeadline("Normalize AddSegments start");

        // NOTE: this is a design flaw, the addline may return a line that is
        // EQUIVALENT BUT NOT IDENTICAL
        // hence line distances mentioned in "segments" DO apply to templine
        // but not necessary (but possibly) to isectLineId
        auto isectLineId = p.AddLine(templine.origin, templine.direction);

        auto &isectLine = p.lines[isectLineId.first];

        if (!p.IsPointOnPlane(isectLine.origin) || !p.IsPointOnPlane(isectLine.origin + isectLine.direction * 100.))
        {
            if (messages)
            {
                printf("Bad isect line in AddSegments\n");
            }
        }

        uint64_t iteration = 0;
        for (auto &seg : segments)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize AddSegments");
            }

            auto pos = templine.GetPosOnLine(seg.first);

            if (!p.aabb.contains(pos))
            {
                if (messages)
                {
                    printf("making pos outside in AddSegments]\n");
                }
            }

            size_t ptA = sp.AddPoint(pos);
            size_t ptB = sp.AddPoint(templine.GetPosOnLine(seg.second));

            if (!p.aabb.contains(sp.points[ptA].location3D))
            {
                if (messages)
                {
                    printf("bad points in AddSegments\n");
                }
            }
            if (!p.aabb.contains(sp.points[ptB].location3D))
            {
                if (messages)
                {
                    printf("bad points in AddS.segments\n");
                }
            }

            // if (ptA != ptB)
            {
                isectLine.AddPointToLine(isectLine.GetPosOnLine(sp.points[ptA].location3D), ptA);
                isectLine.AddPointToLine(isectLine.GetPosOnLine(sp.points[ptB].location3D), ptB);
            }

            if (!p.IsPointOnPlane(sp.points[ptA].location3D))
            {
                if (messages)
                {
                    printf("bad point in AddSegments\n");
                }
            }
            if (!p.IsPointOnPlane(sp.points[ptB].location3D))
            {
                if (messages)
                {
                    printf("bad point in AddSegments\n");
                }
            }

            sp.AddRefPlaneToPoint(ptA, p.id);
            sp.AddRefPlaneToPoint(ptB, p.id);
        }
    }

    /*
        Compute all intersection distances of a candidate line (lineA) with the
        existing segments on a plane. A very long ray along lineA is tested against
        every non-collinear segment on the plane. The resulting intersection
        distances (projected onto lineA) are sorted and deduplicated.
        Returns a sorted vector of distances along lineA.
    */
    inline std::vector<double> ComputeInitialIntersections(Plane &p, SharedPosition &sp, const Line &lineA, const BooleanBudget& budget, double maxPointRadius)
    {
        budget.CheckDeadline("Normalize ComputeInitialIntersections start");

        // The ray along lineA only has to span the whole geometry. The historical code scanned
        // ALL sp.points on EVERY call to find the farthest point from lineA.origin -- an
        // O(points) pass per executed plane pair that dominated Normalize on high-plane-count
        // operands (hollow-core slab: the whole 18.5 s budget burned in this scan). The caller
        // now passes the point-cloud radius (computed once per Normalize); together with
        // |lineA.origin| it upper-bounds the farthest point from the ray origin. The 1e4 floor
        // reproduces the historical minimum (the old scan seeded its max with 1e8 distance^2),
        // so models below that scale get the exact same ray span as before.
        const double size = std::max(1.0E+04, 2.0 * (maxPointRadius + glm::length(lineA.origin)));

        uint64_t iteration = 0;

        auto Astart = lineA.origin + lineA.direction * (size * 2);
        auto Aend = lineA.origin - lineA.direction * (size * 2);

        std::vector<double> distances;
        distances.reserve(p.lines.size());

        // line B is expected to have the segments already filled, line A is not
        for (auto &line : p.lines)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize ComputeInitialIntersections lines");
            }

            // skip collinear
            if (lineA.IsCollinear(line))
                continue;

            for (const auto &seg : line.GetSegments())
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize ComputeInitialIntersections segments");
                }

                auto result = LineLineIntersection(
                    Astart,
                    Aend,
                    sp.points[seg.first].location3D,
                    sp.points[seg.second].location3D);

                if (result.distance < _TOLERANCE_PLANE_INTERSECTION)
                {
                    if (!p.aabb.contains(sp.points[seg.first].location3D))
                    {
                        if (messages)
                        {
                            printf("bad points in ComputeInitialIntersections\n");
                        }
                    }
                    if (!p.aabb.contains(sp.points[seg.second].location3D))
                    {
                        if (messages)
                        {
                            printf("bad points in ComputeInitialIntersections\n");
                        }
                    }

                    // intersection, mark index of line B and distance on line A
                    distances.emplace_back(lineA.GetPosOnLine(result.point2));
                    auto pt = lineA.GetPosOnLine(distances[distances.size() - 1]);

                    if (!p.aabb.contains(result.point2))
                    {
                        if (messages)
                        {
                            printf("bad points in ComputeInitialIntersections\n");
                        }
                    }

                    if (!equals(pt, result.point2, _TOLERANCE_PLANE_INTERSECTION))
                    {
                        if (messages)
                        {
                            printf("BAD POINT in ComputeInitialIntersections\n");
                        }
                    }
                }
            }
        }

        const auto double_less = +[](double left, double right)
        { return left < right; };
        std::sort(distances.begin(), distances.end(), double_less);

        // Epsilon-dedup distances along the line. std::unique on raw doubles
        // misses near-identical intersection positions that differ only by
        // floating-point noise, and the resulting micro-segments cascade
        // into huge face counts downstream (test61 op #8204981 blowup).
        {
            std::vector<double> uniq;
            uniq.reserve(distances.size());
            for (double d : distances)
            {
                if (uniq.empty() || std::abs(d - uniq.back()) > _TOLERANCE_PLANE_INTERSECTION)
                {
                    uniq.push_back(d);
                }
            }
            distances.swap(uniq);
        }

        return distances;
    }

    /*
        Find intersection points between all segment pairs of two lines on a plane.
        For each pair of non-overlapping segments (one from lineA, one from lineB),
        compute their closest-approach point. If the distance is below tolerance,
        the intersection point is added to the shared position and registered on
        both lines and their associated planes.
    */
    inline void AddLineLineIntersections(Plane &p, SharedPosition &sp, Line &lineA, Line &lineB, const BooleanBudget& budget)
    {
        uint64_t iteration = 0;
        for (auto segA : lineA.GetSegments())
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize AddLineLineIntersections A");
            }

            for (auto segB : lineB.GetSegments())
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize AddLineLineIntersections B");
                }

                if (!p.HasOverlap(segA, segB))
                {
                    auto result = LineLineIntersection(
                        sp.points[segA.first].location3D,
                        sp.points[segA.second].location3D,
                        sp.points[segB.first].location3D,
                        sp.points[segB.second].location3D);

                    if (result.distance < SCALED_EPS_BIG)
                    {
                        if (!p.aabb.contains(result.point1))
                        {
                            if (messages)
                                printf("bad points in AddLineLineIntersections\n");
                            continue;
                        }

                        size_t point = sp.AddPoint((result.point1));
                        lineA.AddPointToLine(lineA.GetPosOnLine(sp.points[point].location3D), point);
                        lineB.AddPointToLine(lineB.GetPosOnLine(sp.points[point].location3D), point);

                        {
                            ReferenceLine ref;
                            ref.pointID = point;
                            ref.lineID = lineA.id;
                            ref.location = lineA.GetPosOnLine(result.point1);
                            sp.points[point].lines.push_back(ref);
                            for (auto &plane : lineA.planes)
                                sp.AddRefPlaneToPoint(point, plane.planeID);
                        }
                        {
                            ReferenceLine ref;
                            ref.pointID = point;
                            ref.lineID = lineB.id;
                            ref.location = lineB.GetPosOnLine(result.point2);
                            sp.points[point].lines.push_back(ref);
                            for (auto &plane : lineB.planes)
                                sp.AddRefPlaneToPoint(point, plane.planeID);
                        }
                    }
                }
            }
        }
    }

    /*
        Compute all pairwise line-line intersections on a plane.
        An AABB is pre-computed for each line's extent so that pairs whose
        bounding boxes don't overlap are skipped, avoiding unnecessary
        segment-vs-segment tests.
    */
    inline void AddLineLineIsects(Plane &p, SharedPosition &sp, const BooleanBudget& budget)
    {
        const size_t total = p.lines.size();

        // Pre-compute line bounding boxes for spatial pre-filter.
        // Lines whose AABBs don't overlap cannot have intersecting segments.
        std::vector<AABB> lineAABBs(total);
        uint64_t iteration = 0;
        for (size_t i = 0; i < total; i++) {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize AddLineLineIsects AABBs");
            }

            auto &line = p.lines[i];
            if (line.points.size() >= 2) {
                Vec3 dir = glm::normalize(line.direction);
                lineAABBs[i].merge(line.origin + line.points.front().first * dir);
                lineAABBs[i].merge(line.origin + line.points.back().first * dir);
            }
        }

        for (size_t lineAIndex = 0; lineAIndex < total; lineAIndex++)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize AddLineLineIsects outer");
            }

            for (size_t lineBIndex = lineAIndex + 1; lineBIndex < total; lineBIndex++)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize AddLineLineIsects inner");
                }

                if (!lineAABBs[lineAIndex].intersects(lineAABBs[lineBIndex])) continue;
                AddLineLineIntersections(p, sp, p.lines[lineAIndex], p.lines[lineBIndex], budget);
            }
        }
    }

    /*
        Top-level CSG merging algorithm. Given two input geometries already loaded
        into a SharedPosition (via Construct), this function:
          1. Extracts boundary contour edges from both meshes and adds them as
             lines on their respective planes.
          2. Associates all points with the planes and lines they lie on.
          3. Computes all line-line intersections within each plane.
          4. Intersects every pair of non-parallel, overlapping planes to produce
             new intersection lines and shared segments.
          5. Recomputes line-line intersections after the new lines are added.
          6. Re-triangulates each plane via Constrained Delaunay Triangulation,
             keeping only triangles that lie on a mesh boundary.
          7. Re-adds irrelevant (non-overlapping) faces from both inputs.
        Returns the merged Geometry ready for inside/outside classification.
    */
    inline Geometry Normalize(const Geometry &A, const Geometry &B, SharedPosition &sp, bool UNION, const BooleanBudget& budget)
    {
        budget.CheckDeadline("Normalize start");

        // construct all contours, derive lines
        auto contoursA = sp.A.GetContourSegments();

        uint64_t iteration = 0;
        for (auto &[planeId, contours] : contoursA)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize contours A");
            }

            std::vector<std::vector<glm::dvec2>> edges;

            Plane &p = sp.planes[planeId];

#ifdef CSG_DEBUG_OUTPUT
            // auto basis = p.MakeBasis();

            // for (auto& segment : contours)
            // {
            //     edges.push_back({ basis.project(sp.points[segment.first].location3D), basis.project(sp.points[segment.second].location3D) });
            // }
            // DumpSVGLines(edges, L"contour_A.html");
#endif

            for (auto &segment : contours)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize contour segments A");
                }

                auto lineId = sp.planes[planeId].AddLine(sp.points[segment.first], sp.points[segment.second]);
            }
        }

        auto contoursB = sp.B.GetContourSegments();

        for (auto &[planeId, contours] : contoursB)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize contours B");
            }

            std::vector<std::vector<glm::dvec2>> edges;

            Plane &p = sp.planes[planeId];

#ifdef CSG_DEBUG_OUTPUT
            // auto basis = p.MakeBasis();

            // for (auto& segment : contours)
            // {
            //     edges.push_back({ basis.project(sp.points[segment.first].location3D), basis.project(sp.points[segment.second].location3D) });
            // }
            // DumpSVGLines(edges, L"contour_B.html");
#endif

            for (auto &segment : contours)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize contour segments B");
                }

                auto lineId = sp.planes[planeId].AddLine(sp.points[segment.first], sp.points[segment.second]);
            }
        }

        // put all points on lines/planes
        for (auto &p : sp.points)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize point-plane refs");
            }

            for (auto &plane : sp.planes)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize point-plane refs inner");
                }

                if (plane.IsPointOnPlane(p.location3D))
                {
                    sp.AddRefPlaneToPoint(p.id, plane.id);
                    plane.PutPointOnLines(p);
                }
            }
        }

        for (auto &plane : sp.planes)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize initial line intersections");
            }

            AddLineLineIsects(plane, sp, budget);
        }

        // Gated diagnostics (CXDIAG_CSG): phase stopwatch + pair-loop counters. The budget
        // timeout only names the phase whose deadline check happened to fire; these numbers
        // attribute the real cost.
        const bool normDiag = getenv("CXDIAG_CSG") != nullptr;
        const auto tPairs0 = std::chrono::steady_clock::now();
        size_t dPairsAabbSkip = 0, dPairsParallelSkip = 0, dPairsExec = 0, dPairsSegments = 0;

        // Point-cloud radius for ComputeInitialIntersections ray spans, computed once instead
        // of the historical per-call scan over all shared points (see the comment there).
        // Points created later (segment intersections) lie within the existing cloud up to
        // tolerance; the 2x margin in the per-call bound covers them.
        double maxPointRadius2 = 0.0;
        for (auto &point : sp.points)
        {
            maxPointRadius2 = std::max(maxPointRadius2, glm::length2(point.location3D));
        }
        const double maxPointRadius = std::sqrt(maxPointRadius2);

        // intersect planes
        // Inner loop starts at planeAIndex+1: each unordered pair {A,B} is processed
        // exactly once. AddSegments is called on both planes inside every iteration,
        // so processing (B,A) after (A,B) was purely redundant. This halves the
        // number of plane pairs from N*(N-1) to N*(N-1)/2.
        for (size_t planeAIndex = 0; planeAIndex < sp.planes.size(); planeAIndex++)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize plane-pair outer");
            }

            for (size_t planeBIndex = planeAIndex + 1; planeBIndex < sp.planes.size(); planeBIndex++)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize plane-pair inner");
                }

                auto &planeA = sp.planes[planeAIndex];
                auto &planeB = sp.planes[planeBIndex];

                if (!planeA.aabb.intersects(planeB.aabb))
                {
                    dPairsAabbSkip++;
                    continue;
                }

                // plane intersect results in new lines
                // new lines result in new line intersects
                // new line intersects result in new points

                if (std::fabs(glm::dot(planeA.normal, planeB.normal)) > 1.0 - EPS_BIG)
                {
                    // parallel planes, don't care
                    dPairsParallelSkip++;
                    continue;
                }
                dPairsExec++;

                // calculate plane intersection line
                auto result = PlanePlaneIsect(planeA.normal, planeA.distance, planeB.normal, planeB.distance);

                // Skip degenerate intersections (near-parallel planes that
                // slipped through the dot-product check)
                if (glm::length(result.dir) < EPS_SMALL)
                {
                    continue;
                }

                // TODO: invalid temp line object
                Line intersectionLine;
                intersectionLine.origin = result.pos;
                intersectionLine.direction = result.dir;

                if (!planeA.IsPointOnPlane(intersectionLine.origin) || !planeA.IsPointOnPlane(intersectionLine.origin + intersectionLine.direction * 1000.))
                {
                    if (messages)
                    {
                        printf("Bad isect line in Normalize\n");
                    }
                }
                if (!planeB.IsPointOnPlane(intersectionLine.origin) || !planeB.IsPointOnPlane(intersectionLine.origin + intersectionLine.direction * 1000.))
                {
                    if (messages)
                    {
                        printf("Bad isect line in Normalize\n");
                    }
                }

                // get all intersection points with the shared line and both planes
                auto isectA = ComputeInitialIntersections(planeA, sp, intersectionLine, budget, maxPointRadius);
                auto isectB = ComputeInitialIntersections(planeB, sp, intersectionLine, budget, maxPointRadius);

                // from these, figure out the shared segments on the current line produced by these two planes
                auto segments = sp.BuildSegments(isectA, isectB);

                if (segments.empty())
                {
                    // nothing resulted from this plane-plane intersection
                    continue;
                }
                else
                {
                    dPairsSegments++;
                    AddSegments(planeA, sp, intersectionLine, segments, budget);
                    AddSegments(planeB, sp, intersectionLine, segments, budget);
                }
            }
        }
        const auto tPairs1 = std::chrono::steady_clock::now();

        for (auto &plane : sp.planes)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize final line intersections");
            }

            AddLineLineIsects(plane, sp, budget);
        }
        const auto tLineIsects1 = std::chrono::steady_clock::now();

        for (auto &p : sp.points)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize final point-plane refs");
            }

            for (auto &plane : sp.planes)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize final point-plane refs inner");
                }

                if (plane.IsPointOnPlane(p.location3D))
                {
                    sp.AddRefPlaneToPoint(p.id, plane.id);
                }
            }
        }
        if (normDiag)
        {
            const auto tRefs1 = std::chrono::steady_clock::now();
            auto ms = [](std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b)
            { return std::chrono::duration_cast<std::chrono::milliseconds>(b - a).count(); };
            size_t totalLines = 0, maxLines = 0, totalSegs = 0;
            for (auto& pl : sp.planes)
            {
                totalLines += pl.lines.size();
                maxLines = std::max(maxLines, pl.lines.size());
                for (auto& ln : pl.lines) totalSegs += ln.GetSegments().size();
            }
            std::cerr << "DIAG: Normalize planes=" << sp.planes.size() << " points=" << sp.points.size()
                      << " pairsExec=" << dPairsExec << " pairsSeg=" << dPairsSegments
                      << " aabbSkip=" << dPairsAabbSkip << " parSkip=" << dPairsParallelSkip
                      << " lines(total/max)=" << totalLines << "/" << maxLines << " segs=" << totalSegs
                      << " pairMs=" << ms(tPairs0, tPairs1) << " lineIsectMs=" << ms(tPairs1, tLineIsects1)
                      << " refsMs=" << ms(tLineIsects1, tRefs1) << std::endl;
        }

        // from the inserted geometries, all lines planes and points are now merged into a single set of shared planes lines and points
        // from this starting point, we can triangulate all planes and obtain the triangulation of the intersected set of geometries
        // this mesh itself is not a boolean result, but rather a merging of all operands

        Geometry geom;
        for (auto &plane : sp.planes)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize triangulate planes");
                budget.CheckFaceCount(geom.numFaces, "Normalize triangulate planes");
            }

#ifdef CSG_DEBUG_OUTPUT
            // std::vector<std::vector<glm::dvec2>> edges;

            // auto basis = plane.MakeBasis();

            // for (auto& line : plane.lines) {
            //     // Get line parameters
            //     auto origin = line.origin;      // 3D point (glm::dvec3)
            //     auto direction = line.direction; // 3D vector (glm::dvec3)

            //     // Convert each distance to a 3D point along the line
            //     std::vector<glm::dvec2> lineSegments;
            //     for (auto distance : line.points) {
            //         glm::dvec3 point3D = origin + glm::dvec3(direction.x * distance.first, direction.y * distance.first, direction.z * distance.first);  // 3D calculation
            //         glm::dvec2 point2D = basis.project(point3D);        // Project to 2D
            //         lineSegments.push_back(point2D);
            //     }

            //     // Create edges between consecutive points
            //     for (size_t i = 0; i < lineSegments.size() - 1; ++i) {
            //         edges.push_back({
            //             lineSegments[i],
            //             lineSegments[i + 1]
            //         });
            //     }
            // }

            // DumpSVGLines(edges, L"contour.html");
#endif

            sp.TriangulatePlane(geom, plane, budget);
            budget.CheckFaceCount(geom.numFaces, "Normalize triangulate plane result");

#ifdef CSG_DEBUG_OUTPUT
            // DumpGeometry(geom, L"triangulated.obj");
#endif
        }

        for (auto &plane : A.planes)
        {
            SimplePlane p;
            p.normal = plane.normal;
            p.distance = plane.distance;
            geom.planes.push_back(p);
            geom.hasPlanes = true;
        }

        for (auto &plane : B.planes)
        {
            SimplePlane p;
            p.normal = plane.normal;
            p.distance = plane.distance;
            geom.planes.push_back(p);
            geom.hasPlanes = true;
        }

        uint32_t offsetA = A.planes.size();

        // re-add irrelevant faces that should be tested
        for (auto &faceIndex : sp.A.irrelevantFaces_toTest)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize re-add A test faces");
                budget.CheckFaceCount(geom.numFaces, "Normalize re-add A test faces");
            }

            const Face &f = sp._linkedA->GetFace(faceIndex);

            auto a = sp._linkedA->GetPoint(f.i0);
            auto b = sp._linkedA->GetPoint(f.i1);
            auto c = sp._linkedA->GetPoint(f.i2);

            geom.AddFace(a, b, c, f.pId);
            budget.CheckFaceCount(geom.numFaces, "Normalize re-add A test faces");
        }

        for (auto &faceIndex : sp.B.irrelevantFaces_toTest)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize re-add B test faces");
                budget.CheckFaceCount(geom.numFaces, "Normalize re-add B test faces");
            }

            const Face &f = sp._linkedB->GetFace(faceIndex);

            auto a = sp._linkedB->GetPoint(f.i0);
            auto b = sp._linkedB->GetPoint(f.i1);
            auto c = sp._linkedB->GetPoint(f.i2);

            geom.AddFace(a, b, c, f.pId + offsetA);
            budget.CheckFaceCount(geom.numFaces, "Normalize re-add B test faces");
        }

        geom.data = geom.numFaces;

        // re-add irrelevant faces
        for (auto &faceIndex : sp.A.irrelevantFaces)
        {
            if ((iteration++ & 1023ULL) == 0)
            {
                budget.CheckDeadline("Normalize re-add A irrelevant faces");
                budget.CheckFaceCount(geom.numFaces, "Normalize re-add A irrelevant faces");
            }

            const Face &f = sp._linkedA->GetFace(faceIndex);

            auto a = sp._linkedA->GetPoint(f.i0);
            auto b = sp._linkedA->GetPoint(f.i1);
            auto c = sp._linkedA->GetPoint(f.i2);

            geom.AddFace(a, b, c, f.pId);
            budget.CheckFaceCount(geom.numFaces, "Normalize re-add A irrelevant faces");
        }

        if (UNION)
        {
            for (auto &faceIndex : sp.B.irrelevantFaces)
            {
                if ((iteration++ & 1023ULL) == 0)
                {
                    budget.CheckDeadline("Normalize re-add B irrelevant faces");
                    budget.CheckFaceCount(geom.numFaces, "Normalize re-add B irrelevant faces");
                }

                const Face &f = sp._linkedB->GetFace(faceIndex);

                auto a = sp._linkedB->GetPoint(f.i0);
                auto b = sp._linkedB->GetPoint(f.i1);
                auto c = sp._linkedB->GetPoint(f.i2);

                geom.AddFace(a, b, c, f.pId + offsetA);
                budget.CheckFaceCount(geom.numFaces, "Normalize re-add B irrelevant faces");
            }
        }

        budget.CheckDeadline("Normalize complete");
        budget.CheckFaceCount(geom.numFaces, "Normalize complete");
        return geom;
    }
}
