// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Vector.hpp"


namespace mathter {

template <class T, int Dim, int Order>
class BezierCurve {
	static_assert(Order >= 1, "Bezier curve must have order n>=1.");

public:
	using Vec = Vector<T, Dim>;

	/// <summary> Does not initialize the object. </summary>
	BezierCurve() = default;

	/// <summary> Construct from array of control points. </summary>
	template <class TOther>
	BezierCurve(const BezierCurve<TOther, Dim, Order>& other) {
		auto it = controlPoints.begin();
		for (auto& cp : other.controlPoints) {
			*it++ = Vec(cp);
		}
	}

	/// <summary> Construct from array of control points. </summary>
	explicit BezierCurve(const std::array<Vec, Order + 1>& controlPoints) : controlPoints(controlPoints) {}

	/// <summary> Construct from control points. </summary>
	template <class... Vectors, class = std::enable_if_t<sizeof...(Vectors) == Order + 1>>
	explicit BezierCurve(const Vectors&... controlPoints) : controlPoints{ controlPoints... } {}

	/// <summary> Interpolates the Bezier curve. </summary>
	Vec operator()(const T& t) const {
		return PointAt(t);
	}

	/// <summary> Interpolates the Bezier curve. </summary>
	Vec PointAt(const T& t) const;

public:
	std::array<Vec, Order + 1> controlPoints;
};


template <class Vec, size_t Size>
BezierCurve(const std::array<Vec, Size>&) -> BezierCurve<scalar_type_t<Vec>, dimension_v<Vec>, Size - 1>;


template <class... Vectors, class = std::enable_if_t<sizeof...(Vectors) >= 2>>
BezierCurve(const Vectors&...) -> BezierCurve<scalar_type_t<common_arithmetic_type_t<Vectors...>>, dimension_v<common_arithmetic_type_t<Vectors...>>, sizeof...(Vectors) - 1>;


template <class T, int Dim, int Order>
auto BezierCurve<T, Dim, Order>::PointAt(const T& t) const -> Vec {
	std::array<Vec, Order + 1> reduction = controlPoints;

	const T u = T(1) - t;

	for (int i = Order; i >= 1; --i) {
		for (int j = 1; j <= i; ++j) {
			reduction[j - 1] = u * reduction[j - 1] + t * reduction[j];
		}
	}

	return reduction[0];
}

} // namespace mathter