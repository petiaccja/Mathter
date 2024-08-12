// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Math.hpp"
#include "../Vector/Vector.hpp"

#include <array>


namespace mathter {

template <class T, int Dim>
class Triangle {
	using Vec = Vector<T, Dim>;

public:
	/// <summary> Does not initialize the object. </summary>
	Triangle() = default;

	/// <summary> Construct a triangle from its corners. </summary>
	Triangle(const Vec& a, const Vec& b, const Vec& c)
		: corners{ a, b, c } {}

	/// <summary> Construct a triangle from its corners. </summary>
	explicit Triangle(const std::array<Vec, 3>& corners)
		: corners{ corners } {}

	/// <summary> Convert from a triangle with different scalar type. </summary>
	template <class TOther>
	Triangle(const Triangle<TOther, Dim>& other)
		: corners{ Vec(other.corners[0]), Vec(other.corners[1]), Vec(other.corners[2]) } {}

	/// <summary> Returns the area of the triangle. </summary>
	T Area() const {
		const auto l0 = LengthPrecise(corners[0] - corners[1]);
		const auto l1 = LengthPrecise(corners[0] - corners[2]);
		const auto l2 = LengthPrecise(corners[1] - corners[2]);
		const auto s = T(0.5) * (l0 + l1 + l2);
		return std::sqrt(s * (s - l0) * (s - l1) * (s - l2));
	}

	/// <summary> Returns the centroid of the triangle. </summary>
	/// <remarks> This is also the center of mass for triangles of uniform density. </remarks>
	Vec Centroid() const {
		return (corners[0] + corners[1] + corners[2]) / T(3);
	}

public:
	std::array<Vec, 3> corners;
};

} // namespace mathter