// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Math.hpp"
#include "../Vector/Vector.hpp"

#include <cmath>


namespace mathter {

template <class T, int Dim>
class Line {
	using Vec = Vector<T, Dim>;
	static_assert(!is_complex_v<T>, "Must use real numbers to describe line.");
	static_assert(!std::is_integral_v<T>, "Must use real numbers to describe line.");

public:
	/// <summary> Does not zero-initialize members. </summary>
	Line() = default;

	/// <summary> Converts from a line with different scalar type. </summary>
	template <class TOther>
	Line(const Line<TOther, Dim>& other) : direction(other.direction), base(other.base) {}

	/// <summary> Construct a line through <paramref name="base"/> in given <paramref name="direction"/>. </summary>
	/// <param name="base"> Any point in 3D space. </param>
	/// <param name="direction"> Must be normalized. </param>
	Line(const Vec& base, const Vec& direction) : direction(direction), base(base) {
		assert(std::abs(T(1) - Length(direction)) < T(0.001));
	}

	/// <summary> Constructs a line through both points. </summary>
	/// <param name="point1"> Base of the line. </param>
	/// <param name="point2"> Specifies direction only. </param>
	static Line Through(const Vec& point1, const Vec& point2) {
		return Line(point1, NormalizePrecise(point2 - point1));
	}

	/// <summary> Return the signed direction of the line (as given in constructor). </summary>
	const Vec& Direction() const {
		return direction;
	}

	/// <summary> Returns the base point or point1 as given in constructor. </summary>
	const Vec& Base() const {
		return base;
	}

	/// <summary> Returns the point at <paramref name="param"/> distance from the base point along direction. </summary>
	Vec PointAt(const T& param) const {
		return base + param * direction;
	}

public:
	Vec direction;
	Vec base;
};

} // namespace mathter