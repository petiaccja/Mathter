// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Decompositions/DecomposeQR.hpp"
#include "../Transforms/IdentityBuilder.hpp"
#include "Hyperplane.hpp"
#include "Line.hpp"
#include "LineSegment.hpp"
#include "Ray.hpp"
#include "Triangle.hpp"

#include <cmath>
#include <optional>


namespace mathter {

namespace impl {

	template <class T1, class T2, int Dim>
	auto IntersectHyperplaneLine(const Hyperplane<T1, Dim>& plane, const Line<T2, Dim>& line)
		-> std::optional<common_arithmetic_type_t<T1, T2>> {
		using T = common_arithmetic_type_t<T1, T2>;

		const auto bp = plane.Base();
		const auto& bl = line.base;

		const auto bpMinusBl = bp - bl;
		const auto nDotD = Dot(plane.normal, line.direction);

		if (nDotD == T(0)) {
			return std::nullopt;
		}

		const auto t = Dot(plane.normal, bpMinusBl) / nDotD;
		return t;
	}


	// Assumes the point is on the line.
	template <class T1, class T2, int Dim>
	auto IsWithin(const LineSegment<T1, Dim>& lineSegment, const Vector<T2, Dim>& point) {
		using T = common_arithmetic_type_t<T1, T2>;
		const auto end = lineSegment.point2 - lineSegment.point1;
		const auto absEnd = Abs(end);
		const auto offset = point - lineSegment.point1;
		const auto axisIdx = std::max_element(absEnd.begin(), absEnd.end()) - absEnd.begin();
		const auto limit = end[axisIdx];
		const auto coord = std::copysign(T1(1), limit) * offset[axisIdx];
		return T(0) <= T(coord) && T(coord) <= T(std::abs(limit));
	}


	// Assumes the point is on the line.
	template <class T1, class T2, int Dim>
	auto IsWithin(const Ray<T1, Dim>& ray, const Vector<T2, Dim>& point) {
		using T = common_arithmetic_type_t<T1, T2>;
		const auto absDir = Abs(ray.direction);
		const auto offset = point - ray.base;
		const auto axisIdx = std::max_element(absDir.begin(), absDir.end()) - absDir.begin();
		const auto limit = ray.direction[axisIdx];
		const auto coord = std::copysign(T1(1), limit) * offset[axisIdx];
		return T(0) <= T(coord);
	}

} // namespace impl


//------------------------------------------------------------------------------
// Plane x Line
//------------------------------------------------------------------------------

/// <summary> Find the intersection of a line and a hyperplane. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line is parallel to the plane. </returns>
template <class T1, class T2, int Dim>
auto Intersect(const Hyperplane<T1, Dim>& plane, const Line<T2, Dim>& line)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, Dim, false>> {
	using T = common_arithmetic_type_t<T1, T2>;
	const auto t = impl::IntersectHyperplaneLine(plane, line);
	if (t) {
		return Line<T, Dim>(line).PointAt(*t);
	}
	else {
		return std::nullopt;
	}
}


/// <summary> Find the intersection of a line and a hyperplane. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line is parallel to the plane. </returns>
template <class T1, class T2, int Dim>
auto Intersect(const Line<T1, Dim>& line, const Hyperplane<T2, Dim>& plane)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, Dim, false>> {
	return Intersect(plane, line);
}


//------------------------------------------------------------------------------
// Plane x Line segment
//------------------------------------------------------------------------------

/// <summary> Find the intersection of a line segment and a hyperplane. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segment is parallel to the plane. </returns>
template <class T1, class T2, int Dim>
auto Intersect(const Hyperplane<T1, Dim>& plane, const LineSegment<T2, Dim>& lineSegment)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, Dim, false>> {
	using T = common_arithmetic_type_t<T1, T2>;
	const auto t = impl::IntersectHyperplaneLine(plane, lineSegment.Line());
	if (0 <= t && t < Distance(lineSegment.point1, lineSegment.point2)) {
		return Line<T, Dim>(lineSegment.Line()).PointAt(*t);
	}
	return std::nullopt;
}


/// <summary> Find the intersection of a line segment and a hyperplane. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segment is parallel to the plane. </returns>
template <class T1, class T2, int Dim>
auto Intersect(const LineSegment<T1, Dim>& lineSegment, const Hyperplane<T2, Dim>& plane)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, Dim, false>> {
	return Intersect(plane, lineSegment);
}


//------------------------------------------------------------------------------
// Plane x Ray
//------------------------------------------------------------------------------

/// <summary> Find the intersection of a line segment and a hyperplane. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segment is parallel to the plane. </returns>
template <class T1, class T2, int Dim>
auto Intersect(const Hyperplane<T1, Dim>& plane, const Ray<T2, Dim>& ray)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, Dim, false>> {
	using T = common_arithmetic_type_t<T1, T2>;
	const auto t = impl::IntersectHyperplaneLine(plane, ray.Line());
	if (0 <= t) {
		return Line<T, Dim>(ray.Line()).PointAt(*t);
	}
	return std::nullopt;
}


/// <summary> Find the intersection of a line segment and a hyperplane. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segment is parallel to the plane. </returns>
template <class T1, class T2, int Dim>
auto Intersect(const Ray<T1, Dim>& ray, const Hyperplane<T2, Dim>& plane)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, Dim, false>> {
	return Intersect(plane, ray);
}


//------------------------------------------------------------------------------
// Line<2> x Line<2>
//------------------------------------------------------------------------------

/// <summary> Find the intersection of two lines in 2D. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the lines are parallel. </returns>
template <class T1, class T2>
auto Intersect(const Line<T1, 2>& lhs, const Line<T2, 2>& rhs)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, 2, false>> {
	return Intersect(Hyperplane(lhs), rhs);
}


//------------------------------------------------------------------------------
// LineSegment<2> x LineSegment<2>
//------------------------------------------------------------------------------

/// <summary> Find the intersection of two line segments in 2D. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segments are parallel. </returns>
template <class T1, class T2>
auto Intersect(const LineSegment<T1, 2>& lhs, const LineSegment<T2, 2>& rhs)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, 2, false>> {
	const auto intersection = Intersect(lhs.Line(), rhs.Line());
	if (intersection && impl::IsWithin(lhs, *intersection) && impl::IsWithin(rhs, *intersection)) {
		return intersection;
	}
	return std::nullopt;
}


//------------------------------------------------------------------------------
// Line<2> x LineSegment<2>
//------------------------------------------------------------------------------

/// <summary> Find the intersection of a line and a line segment in 2D. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segments are parallel. </returns>
template <class T1, class T2>
auto Intersect(const Line<T1, 2>& lhs, const LineSegment<T2, 2>& rhs)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, 2, false>> {
	return Intersect(Hyperplane(lhs), rhs);
}


/// <summary> Find the intersection of a line and a line segment in 2D. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segments are parallel. </returns>
template <class T1, class T2>
auto Intersect(const LineSegment<T1, 2>& lhs, const Line<T2, 2>& rhs)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, 2, false>> {
	return Intersect(rhs, lhs);
}


//------------------------------------------------------------------------------
// Ray<3> x Triangle<3>
//------------------------------------------------------------------------------

/// <summary> Find the intersection of a line and a line segment in 2D. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segments are parallel. </returns>
template <class T1, class T2>
auto Intersect(const Ray<T1, 3>& ray, const Triangle<T2, 3>& triangle)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, 3, false>> {
	using T = common_arithmetic_type_t<T1, T2>;

	constexpr T epsilon = std::numeric_limits<T>::epsilon();

	const auto edge1 = triangle.corners[1] - triangle.corners[0];
	const auto edge2 = triangle.corners[2] - triangle.corners[0];

	const auto h = Cross(ray.Direction(), edge2);
	const auto a = Dot(edge1, h);

	if (std::abs(a) < epsilon) {
		return std::nullopt;
	}

	const auto f = T(1) / a;
	const auto s = ray.Base() - triangle.corners[0];

	const auto u = f * Dot(s, h);
	if (u < T(0) || u > T(1)) {
		return std::nullopt;
	}

	const auto q = Cross(s, edge1);
	const auto v = f * Dot(ray.Direction(), q);

	if (v < 0.0 || u + v > 1.0) {
		return std::nullopt;
	}

	const auto t = f * Dot(edge2, q);
	if (t < epsilon) {
		return std::nullopt;
	}
	const auto point = Ray<T, 3>(ray).PointAt(t);
	return point;
}


/// <summary> Find the intersection of a line and a line segment in 2D. </summary>
/// <returns> The coordinates of the intersection point, or nullopt if the line segments are parallel. </returns>
template <class T1, class T2>
auto Intersect(const Triangle<T1, 3>& triangle, const Ray<T2, 3>& ray)
	-> std::optional<Vector<common_arithmetic_type_t<T1, T2>, 3, false>> {
	return Intersect(ray, triangle);
}

} // namespace mathter