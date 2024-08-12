// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Math.hpp"
#include "../Vector/Vector.hpp"
#include "Line.hpp"


namespace mathter {

template <class T, int Dim>
class LineSegment {
	using Vec = Vector<T, Dim>;

public:
	/// <summary> Does not initialize the object. </summary>
	LineSegment() = default;

	/// <summary> Converts from a line segment with different scalar type. </summary>
	template <class TOther>
	LineSegment(const LineSegment<TOther, Dim>& other) : point1(other.point1), point2(other.point2) {}

	/// <summary> Construct a line segment a starting point, a direction and a length. </summary>
	LineSegment(const Vec& base, const Vec& direction, T length)
		: point1(base), point2(base + direction * length) {}

	/// <summary> Constructs a line segment between two points. </summary>
	LineSegment(const Vec& point1, const Vec& point2)
		: point1(point1), point2(point2) {}

	/// <summary> Return the length of the line segment. </summary>
	T Length() const {
		return mathter::Length(point2 - point1);
	}

	/// <summary> Returns the direction of the line segment. </summary>
	Vec Direction() const {
		return NormalizePrecise(point2 - point1);
	}

	/// <summary> Returns the starting point of the line segment. </summary>
	Vec Start() const {
		return point1;
	}

	/// <summary> Returns the end point of the line segment. </summary>
	Vec End() const {
		return point2;
	}

	/// <summary> Interpolates between the start and the end points of the line segment. </summary>
	Vec Interpol(const T& t) const {
		return t * point2 + (T(1) - t) * point1;
	}

	/// <summary> Calculates the interpolation paramter that gives <paramref name="point"/>. </summary>
	/// <remarks> If the point if not actually on the line, it returns the interpolation
	///		parameter to get the closest point to it. </remarks>
	T InterpolOf(const Vec& point) const {
		const auto offset = point - point1;
		const auto direction = point2 - point1;
		const auto length = LengthPrecise(direction);
		return Dot(offset, direction / length) / length;
	}

	/// <summary> Returns a line colinear to the line segment. </summary>
	mathter::Line<T, Dim> Line() const {
		return mathter::Line<T, Dim>{ point1, Direction() };
	}

public:
	Vec point1;
	Vec point2;
};

} // namespace mathter