// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Math.hpp"
#include "../Vector/Vector.hpp"
#include "Line.hpp"

#include <cmath>


namespace mathter {

template <class T, int Dim>
class LineSegment {
	using Vec = Vector<T, Dim>;

public:
	/// <summary> Does not initialize the object. </summary>
	LineSegment() = default;

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

	/// <summary> Returns a line colinear to the line segment. </summary>
	Line<T, Dim> Line() const {
		return mathter::Line<T, Dim>{ point1, Direction() };
	}

public:
	Vec point1, point2;
};

} // namespace mathter