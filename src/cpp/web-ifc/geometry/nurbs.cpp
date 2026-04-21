#include <spdlog/spdlog.h>
#include "representation/IfcGeometry.h"
#include "representation/geometry.h"
#include "nurbs.h"
#include "operations/geometryutils.h"
#include <tinynurbs/tinynurbs.h>
#include <CDT.h>
#include <numeric>
#include <string>

namespace webifc::geometry{

	void Nurbs::fill_geometry(){
		if (!_initialized) {
			return;
		}

		// Model-unit tolerance used both for loop densification and interior
		// subdivision. Keeping them on the same scale means the CDT polygon and
		// the final triangulation agree on how smooth is smooth.
		double const subdivTol = 0.005 / this->scaling;
		spdlog::debug("[Nurbs::fill_geometry] scaling={} subdivTol={}", this->scaling, subdivTol);

		// Step 1: project each bound loop to UV space, preserving per-loop topology
		// so the CDT below can use each loop as a constraint polygon.
		std::vector<std::vector<glm::dvec2>> uvLoops;
		uvLoops.reserve(this->bounds.size());
		for (auto const& bound : this->bounds) {
			std::vector<glm::dvec2> loop;
			loop.reserve(bound.curve.points.size());
			for (auto const& point : bound.curve.points) {
				auto uv {this->inverse_evaluation(point)};
				auto back {tinynurbs::surfacePoint(*this->nurbs, uv.x, uv.y)};
				double residual = glm::distance(glm::dvec3(back), point);
				spdlog::debug("[NurbsBound] in=({:.4f},{:.4f},{:.4f}) uv=({:.6f},{:.6f}) back=({:.4f},{:.4f},{:.4f}) residual={:.6e} maxErr={:.6e}",
					point.x, point.y, point.z, uv.x, uv.y, back.x, back.y, back.z, residual, this->maxError);
				if (residual > this->maxError) {
					spdlog::warn("[NurbsBound] projection failed, residual={:.6e} > maxError={:.6e} at uv=({:.6f},{:.6f}) for pt=({:.4f},{:.4f},{:.4f})",
						residual, this->maxError, uv.x, uv.y, point.x, point.y, point.z);
				}
				if (!loop.empty() && glm::distance(loop.back(), glm::dvec2(uv.x, uv.y)) < 1e-8) continue;
				loop.emplace_back(uv.x, uv.y);
			}
			if (loop.size() >= 2 && glm::distance(loop.front(), loop.back()) < 1e-8) loop.pop_back();
			if (loop.size() < 3) continue;
			uvLoops.push_back(std::move(loop));
		}
		if (uvLoops.empty()) return;

		// Step 1b: densify each UV loop until chord-to-surface deviation drops
		// below subdivTol. Without this pass the CDT polygon is just the chord
		// polygon of the original IFC edge samples, so curved bounds render as
		// visible facets no matter how finely we subdivide the interior.
		for (auto& loop : uvLoops) {
			for (int pass = 0; pass < 5; ++pass) {
				std::vector<glm::dvec2> next;
				next.reserve(loop.size() * 2);
				bool inserted = false;
				for (size_t i = 0; i < loop.size(); ++i) {
					size_t const j = (i + 1) % loop.size();
					glm::dvec2 const& a = loop[i];
					glm::dvec2 const& b = loop[j];
					glm::dvec2 const m = (a + b) * 0.5;
					glm::dvec3 const pa = tinynurbs::surfacePoint(*this->nurbs, a.x, a.y);
					glm::dvec3 const pb = tinynurbs::surfacePoint(*this->nurbs, b.x, b.y);
					glm::dvec3 const pm = tinynurbs::surfacePoint(*this->nurbs, m.x, m.y);
					next.push_back(a);
					if (glm::distance(pm, (pa + pb) * 0.5) > subdivTol) {
						next.push_back(m);
						inserted = true;
					}
				}
				loop = std::move(next);
				if (!inserted) break;
			}
		}

		// Step 2: constrained Delaunay triangulation with bound loops as edge constraints.
		std::vector<CDT::V2d<double>> cdtVerts;
		std::vector<CDT::Edge> cdtEdges;
		for (auto const& loop : uvLoops) {
			size_t const baseIdx = cdtVerts.size();
			for (auto const& uv : loop) {
				cdtVerts.push_back(CDT::V2d<double>::make(uv.x, uv.y));
			}
			size_t const n = loop.size();
			for (size_t i = 0; i < n; ++i) {
				cdtEdges.push_back(CDT::Edge(
					static_cast<uint32_t>(baseIdx + i),
					static_cast<uint32_t>(baseIdx + (i + 1) % n)));
			}
		}
		CDT::RemoveDuplicatesAndRemapEdges(cdtVerts, cdtEdges);
		cdtEdges.erase(
			std::remove_if(cdtEdges.begin(), cdtEdges.end(),
				[](CDT::Edge const& e) { return e.v1() == e.v2(); }),
			cdtEdges.end());
		if (cdtVerts.size() < 3) return;

		CDT::Triangulation<double> cdt(CDT::VertexInsertionOrder::AsProvided);
		try {
			cdt.insertVertices(cdtVerts);
			cdt.insertEdges(cdtEdges);
			cdt.eraseSuperTriangle();
		}
		catch (...) {
			spdlog::warn("[Nurbs::fill_geometry] CDT failed; skipping surface");
			return;
		}

		// Step 3: clip to the bound polygon interior with an even-odd fill test
		// on the original UV loops. This matches the approach used in
		// TriangulateRevolution and is more robust than eraseOuterTrianglesAndHoles
		// when input loops have numerical imprecision.
		auto insideUVPoly = [&](double pu, double pv) -> bool {
			int crossings = 0;
			for (auto const& loop : uvLoops) {
				size_t const n = loop.size();
				for (size_t i = 0, j = n - 1; i < n; j = i++) {
					double yi = loop[i].y, yj = loop[j].y;
					if ((yi > pv) != (yj > pv)) {
						double xCross = loop[j].x + (pv - yj) / (yi - yj) * (loop[i].x - loop[j].x);
						if (pu < xCross) crossings++;
					}
				}
			}
			return (crossings & 1) != 0;
		};

		// Step 4: adaptive subdivision driven by 3D deviation. Each triangle is
		// split into 4 only if the true surface sags from the flat triangle by
		// more than `subdivTol` (in model units). Replaces the previous 3 fixed
		// midpoint passes which multiplied the triangle count by 64 regardless
		// of curvature.
		struct UVTri { glm::dvec2 a, b, c; int depth; };
		int const maxDepth = 5;
		std::vector<UVTri> work;
		work.reserve(cdt.triangles.size());
		for (auto const& tri : cdt.triangles) {
			auto const& v0 = cdt.vertices[tri.vertices[0]];
			auto const& v1 = cdt.vertices[tri.vertices[1]];
			auto const& v2 = cdt.vertices[tri.vertices[2]];
			double const cu = (v0.x + v1.x + v2.x) / 3.0;
			double const cv = (v0.y + v1.y + v2.y) / 3.0;
			if (!insideUVPoly(cu, cv)) continue;
			work.push_back({ glm::dvec2(v0.x, v0.y), glm::dvec2(v1.x, v1.y), glm::dvec2(v2.x, v2.y), 0 });
		}

		while (!work.empty()) {
			UVTri const t = work.back();
			work.pop_back();
			glm::dvec3 const p0 = tinynurbs::surfacePoint(*this->nurbs, t.a.x, t.a.y);
			glm::dvec3 const p1 = tinynurbs::surfacePoint(*this->nurbs, t.b.x, t.b.y);
			glm::dvec3 const p2 = tinynurbs::surfacePoint(*this->nurbs, t.c.x, t.c.y);
			if (t.depth < maxDepth) {
				glm::dvec2 const mab = (t.a + t.b) * 0.5;
				glm::dvec2 const mbc = (t.b + t.c) * 0.5;
				glm::dvec2 const mca = (t.c + t.a) * 0.5;
				glm::dvec3 const pab = tinynurbs::surfacePoint(*this->nurbs, mab.x, mab.y);
				glm::dvec3 const pbc = tinynurbs::surfacePoint(*this->nurbs, mbc.x, mbc.y);
				glm::dvec3 const pca = tinynurbs::surfacePoint(*this->nurbs, mca.x, mca.y);
				double const dev = std::max({
					glm::distance(pab, (p0 + p1) * 0.5),
					glm::distance(pbc, (p1 + p2) * 0.5),
					glm::distance(pca, (p2 + p0) * 0.5)
				});
				if (dev > subdivTol) {
					work.push_back({ t.a, mab, mca, t.depth + 1 });
					work.push_back({ mab, t.b, mbc, t.depth + 1 });
					work.push_back({ mca, mbc, t.c,  t.depth + 1 });
					work.push_back({ mab, mbc, mca, t.depth + 1 });
					continue;
				}
			}
			geometry.AddFace(p0, p1, p2);
		}
	}

	Nurbs::Nurbs(IfcGeometry& geometry, std::vector<IfcBound3D>const & bounds, IfcSurface const& surface, double const scaling)
		: geometry{geometry},
			bounds{bounds},
			bspline_surface{surface.BSplineSurface},
			num_u{this->bspline_surface.ControlPoints.size()},
			num_v{this->bspline_surface.ControlPoints.front().size()},
			scaling{scaling} {
			this->init();
	}

	void Nurbs::init() {
		// Check that the control point grid has sufficient dimensions.
		// We need at least (degree + 1) control points in each direction.
		if (this->num_u < static_cast<size_t>(this->bspline_surface.UDegree) + 1) {
			spdlog::error("Insufficient control point rows: num_u = {} but UDegree = {} requires at least {} rows",
				this->num_u, this->bspline_surface.UDegree, this->bspline_surface.UDegree + 1);
			return; // Or throw an exception.
		}
		if (this->num_v < static_cast<size_t>(this->bspline_surface.VDegree) + 1) {
			spdlog::error("Insufficient control point columns: num_v = {} but VDegree = {} requires at least {} columns",
				this->num_v, this->bspline_surface.VDegree, this->bspline_surface.VDegree + 1);
			return; // Or throw an exception.
		}

		// Validate degrees.
		if (this->bspline_surface.UDegree < 0 || this->bspline_surface.VDegree < 0) {
			spdlog::error("Invalid degree values: UDegree={}, VDegree={}",
				this->bspline_surface.UDegree, this->bspline_surface.VDegree);
			return;
		}

		// Validate control points count.
		size_t expectedCount = this->num_u * this->num_v;
		auto controlPoints = this->get_control_points();
		if (controlPoints.size() != expectedCount) {
			spdlog::error("Control points count mismatch: expected {}, got {}", expectedCount, controlPoints.size());
			return;
		}
		auto weights = this->get_weights();
		if (weights.size() != expectedCount) {
			spdlog::error("Weights count mismatch: expected {}, got {}", expectedCount, weights.size());
			return;
		}

		// Create the NURBS surface.
		this->nurbs = std::make_shared<tinynurbs::RationalSurface3d>(
			static_cast<int>(this->bspline_surface.UDegree),
			static_cast<int>(this->bspline_surface.VDegree),
			this->get_knots(this->bspline_surface.UKnots, this->bspline_surface.UMultiplicity),
			this->get_knots(this->bspline_surface.VKnots, this->bspline_surface.VMultiplicity),
			tinynurbs::array2<glm::dvec3>{this->num_u, this->num_v, controlPoints},
			tinynurbs::array2<double>{this->num_u, this->num_v, weights}
		);

		// Check that the knot vectors are large enough.
		if (this->nurbs->knots_u.size() < static_cast<size_t>(this->nurbs->degree_u + 1)) {
			spdlog::error("Invalid knots_u: size {} is less than degree_u+1 ({})",
				this->nurbs->knots_u.size(), this->nurbs->degree_u + 1);
			return;
		}
		if (this->nurbs->knots_v.size() < static_cast<size_t>(this->nurbs->degree_v + 1)) {
			spdlog::error("Invalid knots_v: size {} is less than degree_v+1 ({})",
				this->nurbs->knots_v.size(), this->nurbs->degree_v + 1);
			return;
		}

		// Helper lambda to check if a knot vector is monotonic increasing.
		auto check_monotonic = [](const std::vector<double>& knots, const std::string& name) -> bool {
			for (size_t i = 1; i < knots.size(); i++) {
				if (knots[i] < knots[i - 1]) {
					spdlog::error("{} is not monotonic increasing at index {} ({} < {})", name, i, knots[i], knots[i - 1]);
					return false;
				}
			}
			return true;
			};

		if (!check_monotonic(this->nurbs->knots_u, "knots_u")) return;
		if (!check_monotonic(this->nurbs->knots_v, "knots_v")) return;

		// Ensure that we have enough knots to set the range.
		if (this->nurbs->knots_u.size() <= this->nurbs->degree_u ||
			this->nurbs->knots_u.size() <= this->nurbs->degree_u + 1) {
			spdlog::error("Not enough knots in knots_u to determine range, size={}, degree_u={}",
				this->nurbs->knots_u.size(), this->nurbs->degree_u);
			return;
		}
		if (this->nurbs->knots_v.size() <= this->nurbs->degree_v ||
			this->nurbs->knots_v.size() <= this->nurbs->degree_v + 1) {
			spdlog::error("Not enough knots in knots_v to determine range, size={}, degree_v={}",
				this->nurbs->knots_v.size(), this->nurbs->degree_v);
			return;
		}
		this->range_knots_u = {
			this->nurbs->knots_u[this->nurbs->degree_u],
			this->nurbs->knots_u[this->nurbs->knots_u.size() - this->nurbs->degree_u - 1]
		};
		this->range_knots_v = {
			this->nurbs->knots_v[this->nurbs->degree_v],
			this->nurbs->knots_v[this->nurbs->knots_v.size() - this->nurbs->degree_v - 1]
		};
		spdlog::debug("[Nurbs::init] degree=({},{}) grid={}x{} u_range=[{},{}] v_range=[{},{}]",
			this->bspline_surface.UDegree, this->bspline_surface.VDegree,
			this->num_u, this->num_v,
			this->range_knots_u.x, this->range_knots_u.y,
			this->range_knots_v.x, this->range_knots_v.y);

		// Compute sample surface points.
		this->ptc = tinynurbs::surfacePoint(*this->nurbs, EPS_TINY, EPS_TINY);
		this->pth = tinynurbs::surfacePoint(*this->nurbs, 1.0, EPS_TINY);
		this->ptv = tinynurbs::surfacePoint(*this->nurbs, EPS_TINY, 1.0);

		if (!std::isfinite(this->ptc.x)) this->ptc.x = 0.0;
		if (!std::isfinite(this->ptc.y)) this->ptc.y = 0.0;
		if (!std::isfinite(this->ptc.z)) this->ptc.z = 0.0;

		if (!std::isfinite(this->pth.x)) this->pth.x = 0.0;
		if (!std::isfinite(this->pth.y)) this->pth.y = 0.0;
		if (!std::isfinite(this->pth.z)) this->pth.z = 0.0;

		if (!std::isfinite(this->ptv.x)) this->ptv.x = 0.0;
		if (!std::isfinite(this->ptv.y)) this->ptv.y = 0.0;
		if (!std::isfinite(this->ptv.z)) this->ptv.z = 0.0;

		// Compute distances for further use.
		this->dh = glm::distance(ptc, pth);
		this->dv = glm::distance(ptc, ptv);
		this->pr = (dh + 1) / (dv + 1);
		if (!std::isfinite(this->pr)) this->pr = 1.0;

		// Scale error tolerances.
		this->minError /= this->scaling;
		this->maxError /= this->scaling;
		_initialized = true;
	}

	std::vector<double> Nurbs::get_weights() const{
		std::vector<double> result(this->num_u * this->num_v);
		std::fill(result.begin(), result.end(), 1.0);
		size_t flat = 0;
		for (auto const& row : this->bspline_surface.Weights)
			for (size_t i = 0; i < row.size(); ++i)
				result[flat++] = row[i];
		return result;
	}

	std::vector<double> Nurbs::get_knots(std::vector<double>const & bs_knots, std::vector<uint32_t> const & bs_mults) const{
		std::vector<double> result;
		auto knots_no_expanded {this->check_knots(bs_knots)};
		auto const num_srf_knots {std::accumulate(bs_mults.begin(), bs_mults.end(), 0.0)};
		result.reserve(num_srf_knots);
		for(size_t knot_i{0}; knot_i < bs_knots.size(); ++knot_i){
			auto const knot {knots_no_expanded[knot_i]};
			auto const knot_mult {bs_mults[knot_i]};
			for(size_t i{0}; i < knot_mult; ++i) result.push_back(knot);
		}
		return result;
	}

	std::vector<glm::dvec3> Nurbs::get_control_points() const{
		std::vector<glm::dvec3> result;
		size_t num_points{0};
		for(auto const& row : this->bspline_surface.ControlPoints) num_points += row.size();
		result.reserve(num_points);
		for(auto const& row : this->bspline_surface.ControlPoints)  std::copy(row.begin(), row.end(), std::back_inserter(result));
		return result;
	}

	auto Nurbs::get_approximation(glm::dvec3 const& pt, uv_point_t const& range_u, uv_point_t const& range_v) const{
		double fU{0.0};
		double fV{0.0};
    int const grid_size {10};
    auto min_distance = std::numeric_limits<double>::max();
		auto const portion_u {std::abs((range_u.y - range_u.x) / grid_size)};
		auto const portion_v {std::abs((range_v.y - range_v.x) / grid_size)};
		auto new_range_u {range_u};
		auto new_range_v {range_v};
		for (int i = 0; i < grid_size; ++i) {
				auto const step_u {portion_u * i};
				auto const u {range_u.x + (step_u != 0.0 ? step_u : std::numeric_limits<double>::epsilon())};
        for (int j = 0; j < grid_size; ++j) {
					auto const step_v {portion_v * j};
					auto const v {range_v.x + (step_v != 0.0 ? step_v : std::numeric_limits<double>::epsilon())};
					auto const pt_grid {tinynurbs::surfacePoint(*this->nurbs, u, v)};
					auto const dist {glm::distance(pt_grid, pt)};
					if (dist < min_distance) {
							min_distance = dist;
							fU = u;
							fV = v;
							new_range_u = {range_u.x + portion_u * i, range_u.x + portion_u * (i + 1)};
							new_range_v = {range_v.x + portion_v * j, range_v.x + portion_v * (j + 1)};
					}
        }
    }
		return std::tuple {min_distance, fU, fV, new_range_u, new_range_v};
	}

	Nurbs::uv_point_t Nurbs::inverse_evaluation(glm::dvec3 const& pt) const
	{
		spdlog::debug("[BSplineInverseEvaluation()]");
		return inverse_method(pt);
	}

	Nurbs::uv_point_t Nurbs::inverse_method(glm::dvec3 const& pt) const
	{
		spdlog::debug("[InverseMethod()]");
		glm::highp_dvec3 pt00{};
		// Seed the search with a 10x10 grid-sampled global nearest point instead of
		// the UV midpoint. Starting at the midpoint caused the rotational hill-climb
		// below to settle in local minima for curved or off-center surfaces, producing
		// bound projections on the wrong side of the patch and mesh tears after CDT.
		auto const seed {this->get_approximation(pt, this->range_knots_u, this->range_knots_v)};
		double fU {std::get<1>(seed)};
		double fV {std::get<2>(seed)};
		auto max_distance {std::get<0>(seed)};
		if (max_distance <= maxError) return {fU, fV};

		double divisor {100.0};
		while (max_distance > maxError && divisor < 10000)
		{
			for (double r = 1; r < 5; r++)
			{
				int round = 0;
				auto mul_divisor {r * r * divisor};
				while (max_distance > minError && round < 3)
				{
					for (double i = 0; i < rotations; i++)
					{
						double rads = (i / rotations) * pi2;
						double incU = glm::sin(rads) / mul_divisor;
						double incV = glm::cos(rads) / mul_divisor;
						if (pr > 1) incV *= pr;
						else incU /= pr;
						while (true)
						{
							double ffU = fU + incU;
							double ffV = fV + incV;
							if (ffU < range_knots_u.x)ffU = range_knots_u.y - (range_knots_u.x - ffU);
							else if (ffU > range_knots_u.y) ffU = range_knots_u.x + (ffU - range_knots_u.y);
							if (ffV < range_knots_v.x) ffV = range_knots_v.y - (range_knots_v.x - ffV);
							else if (ffV > range_knots_v.y) ffV = range_knots_v.x + (ffV - range_knots_v.y);
							pt00 = tinynurbs::surfacePoint(*this->nurbs, ffU, ffV);
							auto const di {glm::distance(pt00, pt)};
							if (di < max_distance)
							{
								max_distance = di;
								fU = ffU;
								fV = ffV;
							}
							else
							{
								break;
							}
						}
					}
					round++;
				}
			}
			divisor *= 3;
		}
		return {fU, fV};
	}

	std::vector<double> Nurbs::get_zscores(std::vector<double> const& knots) const{
		std::vector<double> result(knots.size());
		double mean = std::accumulate(knots.begin(), knots.end(), 0.0) / knots.size();
		double sq_sum = std::inner_product(knots.begin(), knots.end(), knots.begin(), 0.0);
		double variance = sq_sum / knots.size() - mean * mean;
		double stdev = variance > 0.0 ? std::sqrt(variance) : 1.0;
		for (size_t i = 0; i < knots.size(); ++i) {
				result[i] = (knots[i] - mean) / stdev;
		}
		return result;
	}

	std::vector<double> Nurbs::check_knots(std::vector<double> const& knots) const{
		std::vector<double> result(knots.size());
		auto const num_knots {knots.size()};
		if(num_knots == 2){
			result[0] = 0.0;
			result[1] = 1.0;
			return result;
		}
		auto threshold {3.0};
		auto zscores {get_zscores(knots)};
		for(size_t i{0}; i < num_knots; ++i){
			if(std::abs(zscores[i]) > threshold){
				if(i == 0)									result[i] = knots[i+1];
				else if (i == num_knots -1)	result[i] = knots[i-1];
				else 												result[i] = (knots[i-1]+knots[i+1]) / 2.0;
			}
			else result[i] = knots[i];
		}
		return result;
	}
}
