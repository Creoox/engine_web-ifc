#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <optional>

#include <glm/glm.hpp>

#pragma warning(push)
#pragma warning(disable : 4267)
#include <CDT.h>
#pragma warning(pop)

#include "geometry.h"
#include "aabb.h"
#include "bvh.h"
#include "util.h"
#include "svg.h"
#include "math.h"
#include "loop-finder.h"
#include "obj-exporter.h"
#include "is-inside-mesh.h"
#include "is-inside-boundary.h"

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
                        direction is always stored normalised (set via Plane::AddLine or
                        PlanePlaneIsect, both of which normalise before assignment).
            */
            const Vec3 &d = direction;   // already unit-length - no normalize needed
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
            // direction is always unit-length (see AddLine / PlanePlaneIsect)
            return glm::dot(pos - origin, direction);
        }

        Vec3 GetPosOnLine(const double dist) const
        {
            // direction is always unit-length (see AddLine / PlanePlaneIsect)
            return origin + dist * direction;
        }

        bool IsCollinear(const Line &other) const
        {
            // Both directions are guaranteed unit-length (see AddLine / PlanePlaneIsect).
            const Vec3 &unitDirection      = direction;
            const Vec3 &unitOtherDirection = other.direction;
            return (equals(unitOtherDirection, unitDirection, toleranceCollinear) || equals(unitOtherDirection, -unitDirection, toleranceCollinear));
        }

        /*
                The original version of IsEqualTo compared dir and direction.  This version compares
                normalised copies of dir and direction.
        */
        bool IsEqualTo(const Vec3 &pos, const Vec3 &dir) const
        {
            // check dir
            // 'dir' is an external parameter - normalise it.
            // 'direction' is already unit-length (set via AddLine / PlanePlaneIsect).
            Vec3 unitDir = glm::normalize(dir);
            const Vec3 &unitDirection = direction;
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

        void AddPointToLine(double dist, size_t id)
        {
            // check existing
            for (auto &p : points)
            {
                if (p.second == id)
                    return;
            }

            // Maintain sorted order with binary-search insertion (points is always kept sorted
            // by distance, so std::lower_bound finds the right position in O(log N)).
            // This replaces the previous push_back + full re-sort (O(N log N)) with an
            // O(log N) search + O(N) insert, which is significantly faster for small-to-
            // medium lists and avoids redundant comparisons.
            auto it = std::lower_bound(
                points.begin(), points.end(), dist,
                [](const std::pair<double, size_t> &p, double val) { return p.first < val; });
            points.insert(it, std::make_pair(dist, id));
        }

        void AddCoveredInterval(double start, double end)
        {
            if (end < start)
            {
                std::swap(start, end);
            }

            coveredIntervals.emplace_back(start, end);
            std::sort(coveredIntervals.begin(), coveredIntervals.end());

            std::vector<std::pair<double, double>> merged;
            merged.reserve(coveredIntervals.size());

            for (const auto &interval : coveredIntervals)
            {
                if (merged.empty() || interval.first > merged.back().second + TOLERANCE_SCALAR_EQUALITY)
                {
                    merged.push_back(interval);
                }
                else
                {
                    merged.back().second = std::max(merged.back().second, interval.second);
                }
            }

            coveredIntervals.swap(merged);
        }

        bool ContainsParameter(double param, double eps = TOLERANCE_SCALAR_EQUALITY) const
        {
            for (const auto &interval : coveredIntervals)
            {
                if (param >= interval.first - eps && param <= interval.second + eps)
                {
                    return true;
                }
            }

            return false;
        }

        std::optional<std::pair<size_t, size_t>> FindSegmentAt(double param, double eps = TOLERANCE_SCALAR_EQUALITY) const
        {
            for (const auto &interval : coveredIntervals)
            {
                if (param < interval.first - eps || param > interval.second + eps)
                {
                    continue;
                }

                bool haveStart = false;
                size_t startPoint = SIZE_MAX;
                double startParam = 0.0;

                for (const auto &point : points)
                {
                    if (point.first < interval.first - eps || point.first > interval.second + eps)
                    {
                        continue;
                    }

                    if (!haveStart)
                    {
                        haveStart = true;
                        startPoint = point.second;
                        startParam = point.first;
                        continue;
                    }

                    if (param >= startParam - eps && param <= point.first + eps)
                    {
                        return std::make_pair(startPoint, point.second);
                    }

                    startPoint = point.second;
                    startParam = point.first;
                }
            }

            return std::nullopt;
        }

        std::vector<std::pair<size_t, size_t>> GetSegments(double eps = TOLERANCE_SCALAR_EQUALITY) const
        {
            std::vector<std::pair<size_t, size_t>> segments;

            for (const auto &interval : coveredIntervals)
            {
                bool haveStart = false;
                size_t startPoint = SIZE_MAX;

                for (const auto &point : points)
                {
                    if (point.first < interval.first - eps || point.first > interval.second + eps)
                    {
                        continue;
                    }

                    if (!haveStart)
                    {
                        haveStart = true;
                        startPoint = point.second;
                        continue;
                    }

                    if (startPoint != point.second)
                    {
                        segments.emplace_back(startPoint, point.second);
                    }

                    startPoint = point.second;
                }
            }

            return segments;
        }

        std::vector<std::pair<double, size_t>> points;
        std::vector<std::pair<double, double>> coveredIntervals;

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
        size_t lineLineIsectCheckedUpTo = 0; // tracks how many lines were checked by AddLineLineIsects
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

            const double distA = lines[lineId.first].GetPosOnLine(a.location3D);
            const double distB = lines[lineId.first].GetPosOnLine(b.location3D);
            lines[lineId.first].AddPointToLine(distA, a.id);
            lines[lineId.first].AddPointToLine(distB, b.id);
            lines[lineId.first].AddCoveredInterval(distA, distB);

            return lineId;
        }

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
            return equals(distance, d, toleranceVectorEquality);
        }

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

    struct PairHash {
        size_t operator()(const std::pair<size_t, size_t> &p) const {
            size_t h = std::hash<size_t>()(p.first);
            h ^= std::hash<size_t>()(p.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
            return h;
        }
    };

    struct SegmentSet
    {
        std::vector<std::pair<size_t, size_t>> segments;
        std::vector<Triangle> triangles;
        std::unordered_map<std::pair<size_t, size_t>, size_t, PairHash> segmentCounts;
        std::vector<size_t> irrelevantFaces;
        std::vector<size_t> irrelevantFaces_toTest;

        std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> planeSegments;
        std::unordered_map<size_t, std::unordered_map<std::pair<size_t, size_t>, size_t, PairHash>> planeSegmentCounts;

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

        std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> GetContourSegments()
        {
            std::unordered_map<size_t, std::vector<std::pair<size_t, size_t>>> contours;

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

        // Precomputed conservative bounding radius used by ComputeInitialIntersections
        // to avoid an O(N_points) scan for every plane-pair.  Set once in Normalize()
        // before entering the plane×plane loop.
        double boundingRadius = 0.0;

        // TODO: design flaw
        const Geometry *_linkedA;
        const Geometry *_linkedB;

        Geometry relevantA;
        Geometry relevantB;

        BVH relevantBVHA;
        BVH relevantBVHB;

        std::unordered_map<std::tuple<int64_t, int64_t, int64_t>, std::vector<size_t>> pointGrid;

        //      assumes all triangleIds are connected to base with an edge and are flipped correctly
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

        // simplify, see FindUppermostTriangleId
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

        size_t AddPoint(const Vec3& newPoint)
        {
            // 1. Compute the grid cell for the query point
            const double cellSize = toleranceVectorEquality;

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
                                if (points[existingId] == newPoint) {
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

        void Construct(const Geometry &A, const Geometry &B, bool isUnion)
        {
            auto boxA = A.GetAABB();
            auto boxB = B.GetAABB();

            AddGeometry(A, B, boxB, true, isUnion, 0);
            AddGeometry(B, A, boxA, false, isUnion, A.planes.size());

            _linkedA = &A;
            _linkedB = &B;
        }

        void AddGeometry(const Geometry &geom, const Geometry &secondGeom, const AABB &relevantBounds, bool isA, bool isUnion, uint32_t offsetPlane)
        {
#ifdef CSG_DEBUG_OUTPUT
            Geometry relevant;
#endif

            for (size_t i = 0; i < geom.numFaces; i++)
            {
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

        // pair of lineID, distance
        std::vector<std::pair<double, double>> BuildSegments(const std::vector<double> &a, const std::vector<double> &b) const
        {
            auto buildIntervals = [](const std::vector<double> &values)
            {
                if (values.size() < 2)
                {
                    return std::vector<std::pair<double, double>>{};
                }

                std::vector<double> uniqueValues;
                uniqueValues.reserve(values.size());

                for (double value : values)
                {
                    if (uniqueValues.empty() || std::fabs(uniqueValues.back() - value) > _TOLERANCE_PLANE_INTERSECTION)
                    {
                        uniqueValues.push_back(value);
                    }
                }

                std::vector<std::pair<double, double>> intervals;
                intervals.reserve(uniqueValues.size() / 2);

                for (size_t i = 1; i < uniqueValues.size(); i += 2)
                {
                    const double start = uniqueValues[i - 1];
                    const double end = uniqueValues[i];
                    if (end > start + _TOLERANCE_PLANE_INTERSECTION)
                    {
                        intervals.emplace_back(start, end);
                    }
                }

                return intervals;
            };

            auto intervalsA = buildIntervals(a);
            auto intervalsB = buildIntervals(b);

            if (intervalsA.empty() || intervalsB.empty())
            {
                return {};
            }

            std::vector<std::pair<double, double>> result;
            size_t idxA = 0;
            size_t idxB = 0;

            while (idxA < intervalsA.size() && idxB < intervalsB.size())
            {
                const double start = std::max(intervalsA[idxA].first, intervalsB[idxB].first);
                const double end = std::min(intervalsA[idxA].second, intervalsB[idxB].second);

                if (end > start + _TOLERANCE_PLANE_INTERSECTION)
                {
                    result.emplace_back(start, end);
                }

                if (intervalsA[idxA].second < intervalsB[idxB].second - _TOLERANCE_PLANE_INTERSECTION)
                {
                    idxA++;
                }
                else if (intervalsB[idxB].second < intervalsA[idxA].second - _TOLERANCE_PLANE_INTERSECTION)
                {
                    idxB++;
                }
                else
                {
                    idxA++;
                    idxB++;
                }
            }

            return result;
        }

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

        void TriangulatePlane(Geometry &geom, Plane &p)
        {

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

            for (auto &pointId : pointsOnPlane)
            {
                pointToProjectedPoint[pointId] = projectedPoints.size();
                projectedPointToPoint[projectedPoints.size()] = pointId;
                projectedPoints.push_back(basis.project(points[pointId].location3D));
            }

            std::set<std::pair<size_t, size_t>> edges;

            for (auto &line : p.lines)
            {
                // these segments might intersect internally, lets resolve that so we get a valid chain
                auto segments = GetNonIntersectingSegments(line);

                for (auto &segment : segments)
                {
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
                cdt_verts.emplace_back(CDT::V2d<double>::make(point.x, point.y));
            }

            for (auto &edge : edges)
            {
                cdt_edges.emplace_back((uint32_t)edge.first, (uint32_t)edge.second);
            }

            auto mapping = CDT::RemoveDuplicatesAndRemapEdges(cdt_verts, cdt_edges).mapping;

            cdt.insertVertices(cdt_verts);
            cdt.insertEdges(cdt_edges);

            cdt.eraseSuperTriangle();

            auto triangles = cdt.triangles;

            // Copy constraint edges into a contiguous vector so isInsideBoundary()
            // iterates sequentially (cache-friendly) instead of chasing BST pointers.
            // The crossing count in isInsideBoundary is order-independent (each edge
            // contributes independently), so the result is identical to using the set.
            // The set itself is unchanged and was already used for CDT edge insertion.
            std::vector<std::pair<size_t, size_t>> edgesVec(edges.begin(), edges.end());

            // auto contourLoop = FindLargestEdgeLoop(projectedPoints, edges);

#ifdef CSG_DEBUG_OUTPUT
            // std::vector<std::vector<glm::dvec2>> edges3DTriangles;
            // std::set<std::pair<size_t, size_t>> edgesTriangles;
            // std::set<std::pair<size_t, size_t>> finalEdgesTriangles;
#endif

            for (auto &tri : triangles)
            {
#ifdef CSG_DEBUG_OUTPUT
                // edgesTriangles.insert(std::make_pair(tri.vertices[0], tri.vertices[1]));
                // edgesTriangles.insert(std::make_pair(tri.vertices[1], tri.vertices[2]));
                // edgesTriangles.insert(std::make_pair(tri.vertices[0], tri.vertices[2]));
#endif

                size_t pointIdA = projectedPointToPoint[mapping[tri.vertices[0]]];
                size_t pointIdB = projectedPointToPoint[mapping[tri.vertices[1]]];
                size_t pointIdC = projectedPointToPoint[mapping[tri.vertices[2]]];

                auto ptA = points[pointIdA].location3D;
                auto ptB = points[pointIdB].location3D;
                auto ptC = points[pointIdC].location3D;

                // Skip degenerate/sliver triangles: area check is cheaper than
                // 3 normalizations + 3 dot products (the old toleranceThinTriangle
                // approach). Cross product magnitude = 2*area; threshold matches
                // the old angular tolerance geometrically.
                {
                    Vec crossTri = glm::cross(ptB - ptA, ptC - ptA);
                    double crossLen2 = glm::dot(crossTri, crossTri);
                    double maxEdge2 = std::max({glm::distance2(ptA, ptB),
                                                glm::distance2(ptB, ptC),
                                                glm::distance2(ptC, ptA)});
                    // Minimum altitude = 2*area / maxEdge = |cross| / maxEdge
                    // Skip if altitude^2 / maxEdge < tol^2  =>  |cross|^2 < tol^2 * maxEdge^2
                    if (maxEdge2 > 0 && crossLen2 < toleranceThinTriangle * toleranceThinTriangle * maxEdge2)
                        continue;
                }

                auto triCenter = (ptA + ptB + ptC) / 3.0;

                Vec raydir = computeNormal(ptA, ptB, ptC);

                auto posA = isInsideMesh(triCenter, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                auto posB = isInsideMesh(triCenter, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                if (posA.loc != MeshLocation::BOUNDARY && posB.loc != MeshLocation::BOUNDARY)
                {
                    continue;
                }

                // If the 2D triangle is not inside the boundaries of the projected boundary of the face it requires further verification
                // It can't be discarded because inside/outside could fail when boundaries have internal partitions
                // Therefore new tests are required to verify that the triangle is on the boundary of A or B

                glm::dvec2 t1 = projectedPoints[tri.vertices[0]];
                glm::dvec2 t2 = projectedPoints[tri.vertices[1]];
                glm::dvec2 t3 = projectedPoints[tri.vertices[2]];

                bool inside2d = isInsideBoundary(t1, t2, t3, edgesVec, projectedPoints);

                if (!inside2d)
                {
                    // posA/posB already tested triCenter above — skip redundant re-check
                    // and go straight to the vertex-interpolated probe points.

                    auto ptt = glm::mix(triCenter, ptA, triangleEvaluationFactor);

                    auto probeA = isInsideMesh(ptt, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                    auto probeB = isInsideMesh(ptt, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                    if (probeA.loc != MeshLocation::BOUNDARY && probeB.loc != MeshLocation::BOUNDARY)
                    {
                        continue;
                    }

                    ptt = glm::mix(triCenter, ptB, triangleEvaluationFactor);

                    probeA = isInsideMesh(ptt, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                    probeB = isInsideMesh(ptt, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                    if (probeA.loc != MeshLocation::BOUNDARY && probeB.loc != MeshLocation::BOUNDARY)
                    {
                        continue;
                    }

                    ptt = glm::mix(triCenter, ptC, triangleEvaluationFactor);

                    probeA = isInsideMesh(ptt, glm::dvec3(0), relevantA, relevantBVHA, raydir);
                    probeB = isInsideMesh(ptt, glm::dvec3(0), relevantB, relevantBVHB, raydir);

                    if (probeA.loc != MeshLocation::BOUNDARY && probeB.loc != MeshLocation::BOUNDARY)
                    {
                        continue;
                    }
                }

                // although CDT is great, it spits out too many or too little tris, we fix it manually
                // if (!IsPointInsideLoop(projectedPoints, contourLoop, triCenter))
                //{
                //    printf("removing point outside loop\n");
                //    continue;
                //}

                // TODO: why is this swapped? winding doesnt matter much, but still
                geom.AddFace(ptB, ptA, ptC, p.refPlane);

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
            // Check for duplicate: a point is typically on 1-3 planes,
            // so the linear scan of the small vector is faster than a set.
            auto &refs = points[point].planes;
            for (const auto &r : refs)
                if (r.planeID == plane) return;

            ReferencePlane ref;
            ref.pointID = point;
            ref.planeID = plane;
            refs.push_back(ref);
            planeToPoints[plane].push_back(point);
        }
    };

    inline void AddSegments(Plane &p, SharedPosition &sp, Line &templine, const std::vector<std::pair<double, double>> &segments)
    {
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

        for (auto &seg : segments)
        {
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
                const double distA = isectLine.GetPosOnLine(sp.points[ptA].location3D);
                const double distB = isectLine.GetPosOnLine(sp.points[ptB].location3D);
                isectLine.AddPointToLine(distA, ptA);
                isectLine.AddPointToLine(distB, ptB);
                isectLine.AddCoveredInterval(distA, distB);
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

    inline std::vector<double> ComputeInitialIntersections(Plane &p, SharedPosition &sp, const Line &lineA)
    {
        // -----------------------------------------------------------------------
        // O(L) algorithm: one infinite-line/line intersection per plane-line,
        // instead of one segment/line intersection per segment (O(L×S)).
        //
        // Two non-parallel, non-collinear lines in 3D have at most ONE point of
        // nearest approach.  Because all lines here lie in the same plane and
        // lineA also lies in that plane, collinear lines are filtered out by
        // IsCollinear; the remaining lines each intersect lineA exactly once.
        //
        // The classic "shortest distance between two lines" formula gives us the
        // parameter t along lineA and u along the candidate line in O(1).  We
        // then verify the distance is within tolerance and that u falls on one
        // of the line's real covered intervals, avoiding false positives from
        // unrelated collinear segments that happen to share the same infinite line.
        // -----------------------------------------------------------------------
        std::vector<double> distances;
        distances.reserve(p.lines.size());

        const Vec3 &da = lineA.direction;   // unit-length (guaranteed by AddLine / PlanePlaneIsect)

        for (auto &line : p.lines)
        {
            if (lineA.IsCollinear(line)) continue;
            if (line.points.size() < 2)  continue;   // no segments yet

            const Vec3 &db = line.direction;           // unit-length
            const Vec3  dc = line.origin - lineA.origin;

            // denom = 1 − (da·db)²  (= |da × db|² since |da|=|db|=1)
            const double daDb  = glm::dot(da, db);
            const double denom = 1.0 - daDb * daDb;

            // Guard against nearly-parallel lines (collinear already handled above).
            if (std::fabs(denom) < toleranceParallelLineDenominator) continue;

            const double dcDa = glm::dot(dc, da);
            const double dcDb = glm::dot(dc, db);

            // t : parameter on lineA at the nearest-approach point
            // u : parameter on 'line' at the nearest-approach point
            const double t = (dcDa - dcDb * daDb) / denom;
            const double u = (dcDa * daDb - dcDb) / denom;

            // Distance between the nearest-approach points must be within tolerance.
            const Vec3 pointOnA = lineA.origin + t * da;
            const Vec3 pointOnB = line.origin  + u * db;
            if (glm::distance(pointOnA, pointOnB) >= _TOLERANCE_PLANE_INTERSECTION) continue;

            // Check that the intersection lies on a real covered segment, not just
            // somewhere inside the line's overall parameter span.
            if (!line.ContainsParameter(u, _TOLERANCE_PLANE_INTERSECTION)) continue;

            distances.emplace_back(t);
        }

        const auto double_less = +[](double left, double right)
        { return left < right; };
        std::sort(distances.begin(), distances.end(), double_less);
        distances.erase(std::unique(distances.begin(), distances.end()), distances.end());

        return distances;
    }

    inline void AddLineLineIntersections(Plane &p, SharedPosition &sp, Line &lineA, Line &lineB)
    {
        // Two non-parallel, non-collinear lines in the same plane intersect at most ONCE.
        // Use the infinite-line nearest-approach formula (same as ComputeInitialIntersections)
        // to find that single point in O(1), then verify it falls within an actual segment
        // of each line in O(Sa + Sb).

        if (lineA.points.size() < 2 || lineB.points.size() < 2) return;
        if (lineA.IsCollinear(lineB)) return;

        const Vec3 &da = lineA.direction;          // unit-length (guaranteed by AddLine / PlanePlaneIsect)
        const Vec3 &db = lineB.direction;          // unit-length
        const Vec3  dc = lineB.origin - lineA.origin;

        const double daDb  = glm::dot(da, db);
        const double denom = 1.0 - daDb * daDb;   // |da × db|²
        if (std::fabs(denom) < toleranceParallelLineDenominator) return; // nearly parallel

        const double dcDa = glm::dot(dc, da);
        const double dcDb = glm::dot(dc, db);
        const double t    = (dcDa - dcDb * daDb) / denom;   // parameter on lineA at nearest approach
        const double u    = (dcDa * daDb - dcDb) / denom;   // parameter on lineB at nearest approach

        // Verify the nearest-approach points are within geometric tolerance of each other.
        const Vec3 pointOnA = lineA.origin + t * da;
        const Vec3 pointOnB = lineB.origin + u * db;
        if (glm::distance(pointOnA, pointOnB) >= SCALED_EPS_BIG) return;

        if (!lineA.ContainsParameter(t, SCALED_EPS_BIG)) return;
        if (!lineB.ContainsParameter(u, SCALED_EPS_BIG)) return;

        const auto segA = lineA.FindSegmentAt(t, SCALED_EPS_BIG);
        if (!segA) return;
        const auto segB = lineB.FindSegmentAt(u, SCALED_EPS_BIG);
        if (!segB) return;

        // Skip if the containing segments already share an endpoint - the intersection
        // point is already part of both lines.
        if (p.HasOverlap(*segA, *segB)) return;

        if (!p.aabb.contains(pointOnA))
        {
            if (messages) printf("bad points in AddLineLineIntersections\n");
            return;
        }

        const size_t point = sp.AddPoint(pointOnA);

        lineA.AddPointToLine(lineA.GetPosOnLine(sp.points[point].location3D), point);
        lineB.AddPointToLine(lineB.GetPosOnLine(sp.points[point].location3D), point);

        // Register the new point on the lines and their planes.
        {
            ReferenceLine ref;
            ref.pointID  = point;
            ref.lineID   = lineA.id;
            ref.location = lineA.GetPosOnLine(pointOnA);
            sp.points[point].lines.push_back(ref);

            for (auto &plane : lineA.planes)
                sp.AddRefPlaneToPoint(point, plane.planeID);
        }
        {
            ReferenceLine ref;
            ref.pointID  = point;
            ref.lineID   = lineB.id;
            ref.location = lineB.GetPosOnLine(pointOnB);
            sp.points[point].lines.push_back(ref);

            for (auto &plane : lineB.planes)
                sp.AddRefPlaneToPoint(point, plane.planeID);
        }
    }

    inline void AddLineLineIsects(Plane &p, SharedPosition &sp)
    {
        // Only check pairs involving at least one NEW line (added since last call).
        // Previously-checked pairs (both indices < checkedUpTo) are skipped.
        const size_t prevChecked = p.lineLineIsectCheckedUpTo;
        const size_t total = p.lines.size();
        for (size_t lineAIndex = 0; lineAIndex < total; lineAIndex++)
        {
            // If lineA is old, only pair with new lines (lineBIndex >= prevChecked).
            // If lineA is new, pair with all subsequent lines.
            size_t startB = (lineAIndex < prevChecked) ? std::max(lineAIndex + 1, prevChecked) : lineAIndex + 1;
            for (size_t lineBIndex = startB; lineBIndex < total; lineBIndex++)
            {
                AddLineLineIntersections(p, sp, p.lines[lineAIndex], p.lines[lineBIndex]);
            }
        }
        p.lineLineIsectCheckedUpTo = total;
    }

    inline Geometry Normalize(const Geometry &A, const Geometry &B, SharedPosition &sp, bool UNION)
    {
        // construct all contours, derive lines
        auto contoursA = sp.A.GetContourSegments();
        for (auto &[planeId, contours] : contoursA)
        {
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
                auto lineId = sp.planes[planeId].AddLine(sp.points[segment.first], sp.points[segment.second]);
            }
        }

        auto contoursB = sp.B.GetContourSegments();
        for (auto &[planeId, contours] : contoursB)
        {
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
                auto lineId = sp.planes[planeId].AddLine(sp.points[segment.first], sp.points[segment.second]);
            }
        }

        // put all points on lines/planes
        // AABB pre-filter: each plane's aabb covers only its finite geometry.
        // For typical models each point lies near 1-3 planes, reducing this
        // from O(N_points × P_planes) to O(N_points × ~3).
        for (auto &p : sp.points)
        {
            for (auto &plane : sp.planes)
            {
                if (!plane.aabb.contains(p.location3D)) continue;
                if (plane.IsPointOnPlane(p.location3D))
                {
                    sp.AddRefPlaneToPoint(p.id, plane.id);
                    plane.PutPointOnLines(p);
                }
            }
        }

        for (auto &plane : sp.planes)
        {
            AddLineLineIsects(plane, sp);
        }

        // Precompute a conservative bounding radius once so ComputeInitialIntersections
        // does not have to scan all sp.points (O(N_points)) for every plane pair.
        //
        // Strategy: build the AABB of all current points, then use
        //   radius = (half-diagonal of AABB) + (distance of AABB centre from origin)
        // This guarantees the extension covers the whole geometry even when the
        // plane-plane intersection origin lies outside the AABB.
        {
            AABB globalBounds;
            for (const auto &pt : sp.points)
                globalBounds.merge(pt.location3D);

            const Vec3 center      = (globalBounds.min + globalBounds.max) * 0.5;
            const double halfDiag  = glm::length(globalBounds.max - center);
            const double centerDist = glm::length(center);

            // ×4 safety factor covers intersection-line origins that may lie slightly
            // outside the bounding box (e.g. from nearly-parallel plane pairs).
            sp.boundingRadius = std::max((halfDiag + centerDist) * 4.0, 1.0E+04);
        }

        // Remember how many points exist before the plane×plane loop so that the
        // post-loop point->plane reference sweep (below) only needs to handle the
        // NEW points.  The original points already had their plane references added
        // in the "put all points on lines/planes" loop above.
        const size_t pointCountBeforePlanePlane = sp.points.size();

        // intersect planes
        // NOTE: the inner loop starts at planeAIndex+1 (not 0) so each unordered pair
        // {A,B} is processed exactly once.  The previous loop over all (A,B) and (B,A)
        // called AddSegments on both planes inside every iteration, so processing (B,A)
        // after (A,B) was purely redundant (all points/lines are de-duplicated on insert).
        // Fixing this halves the number of plane pairs: N*(N-1)/2 instead of N*(N-1).
        for (size_t planeAIndex = 0; planeAIndex < sp.planes.size(); planeAIndex++)
        {
            for (size_t planeBIndex = planeAIndex + 1; planeBIndex < sp.planes.size(); planeBIndex++)
            {
                auto &planeA = sp.planes[planeAIndex];
                auto &planeB = sp.planes[planeBIndex];

                if (!planeA.aabb.intersects(planeB.aabb))
                {
                    continue;
                }

                // plane intersect results in new lines
                // new lines result in new line intersects
                // new line intersects result in new points

                if (std::fabs(glm::dot(planeA.normal, planeB.normal)) > 1.0 - toleranceParallelPlaneDot)
                {
                    // parallel planes, don't care
                    continue;
                }

                // calculate plane intersection line
                auto result = PlanePlaneIsect(planeA.normal, planeA.distance, planeB.normal, planeB.distance);

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
                auto isectA = ComputeInitialIntersections(planeA, sp, intersectionLine);
                auto isectB = ComputeInitialIntersections(planeB, sp, intersectionLine);

                // from these, figure out the shared segments on the current line produced by these two planes
                auto segments = sp.BuildSegments(isectA, isectB);

                if (segments.empty())
                {
                    // nothing resulted from this plane-plane intersection
                    continue;
                }
                else
                {
                    AddSegments(planeA, sp, intersectionLine, segments);
                    AddSegments(planeB, sp, intersectionLine, segments);
                }
            }
        }

        for (auto &plane : sp.planes)
        {
            AddLineLineIsects(plane, sp);
        }

        // Only check NEWLY added points (indices >= pointCountBeforePlanePlane).
        // Points that existed before the plane×plane loop already received their
        // plane references in the "put all points on lines/planes" loop above.
        // AddSegments also calls AddRefPlaneToPoint for the two planes it touches,
        // so the remaining gap is other planes these new points may coincide with.
        for (size_t ptIdx = pointCountBeforePlanePlane; ptIdx < sp.points.size(); ptIdx++)
        {
            const auto &p = sp.points[ptIdx];
            for (auto &plane : sp.planes)
            {
                if (!plane.aabb.contains(p.location3D)) continue;
                if (plane.IsPointOnPlane(p.location3D))
                {
                    sp.AddRefPlaneToPoint(p.id, plane.id);
                }
            }
        }

        // from the inserted geometries, all lines planes and points are now merged into a single set of shared planes lines and points
        // from this starting point, we can triangulate all planes and obtain the triangulation of the intersected set of geometries
        // this mesh itself is not a boolean result, but rather a merging of all operands

        Geometry geom;
        for (auto &plane : sp.planes)
        {

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

            sp.TriangulatePlane(geom, plane);

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
            const Face &f = sp._linkedA->GetFace(faceIndex);

            auto a = sp._linkedA->GetPoint(f.i0);
            auto b = sp._linkedA->GetPoint(f.i1);
            auto c = sp._linkedA->GetPoint(f.i2);

            geom.AddFace(a, b, c, f.pId);
        }

        for (auto &faceIndex : sp.B.irrelevantFaces_toTest)
        {
            const Face &f = sp._linkedB->GetFace(faceIndex);

            auto a = sp._linkedB->GetPoint(f.i0);
            auto b = sp._linkedB->GetPoint(f.i1);
            auto c = sp._linkedB->GetPoint(f.i2);

            geom.AddFace(a, b, c, f.pId + offsetA);
        }

        geom.data = geom.numFaces;

        // re-add irrelevant faces
        for (auto &faceIndex : sp.A.irrelevantFaces)
        {
            const Face &f = sp._linkedA->GetFace(faceIndex);

            auto a = sp._linkedA->GetPoint(f.i0);
            auto b = sp._linkedA->GetPoint(f.i1);
            auto c = sp._linkedA->GetPoint(f.i2);

            geom.AddFace(a, b, c, f.pId);
        }

        if (UNION)
        {
            for (auto &faceIndex : sp.B.irrelevantFaces)
            {
                const Face &f = sp._linkedB->GetFace(faceIndex);

                auto a = sp._linkedB->GetPoint(f.i0);
                auto b = sp._linkedB->GetPoint(f.i1);
                auto c = sp._linkedB->GetPoint(f.i2);

                geom.AddFace(a, b, c, f.pId + offsetA);
            }
        }

        return geom;
    }
}
