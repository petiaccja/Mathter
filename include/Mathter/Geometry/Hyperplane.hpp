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
class Hyperplane {
	using Vec = Vector<T, Dim>;

public:
	/// <summary> Does not initialize the object. </summary>
	Hyperplane() = default;

	/// <summary> Converts from a hyperplane with different scalar type. </summary>
	template <class TOther>
	Hyperplane(const Hyperplane<TOther, Dim>& other) : normal(other.normal), scalar(static_cast<T>(other.scalar)) {}

	/// <summary> Construct a plane through a point and a vector normal to the plane. </summary>
	Hyperplane(const Vec& base, const Vec& normal) : normal(normal), scalar(Dot(normal, base)) {
		assert(std::abs(T(1) - Length(normal)) < 0.0001f);
	}

	/// <summary> Construct a plane by its algebraic equation. </summary>
	Hyperplane(const Vec& normal, const T& scalar) : normal(normal), scalar(scalar) {
		assert(std::abs(T(1) - Length(normal)) < 0.0001f);
	}

	/// <summary> Convert a 2D line to a 2D hyperplane, as they are the same objects. </summary>
	explicit Hyperplane(const Line<T, 2>& line)
		: normal(Cross(line.Direction())), scalar(Dot(normal, line.Base())) {
		static_assert(Dim == 2, "Plane dimension must be two, which is a line.");
	}

	/// <summary> Convert a 2D hyperplane to a 2D line, as they are the same objects. </summary>
	explicit operator Line<T, 2>() const {
		static_assert(Dim == 2, "Plane dimension must be two, which is a line.");
		return Line<T, 2>(Base(), Cross(normal));
	}

	/// <summary> Returns the hyperplane's normal vector. </summary>
	const Vec& Normal() const {
		return normal;
	}

	/// <summary> Returns the scalar part of the hyperplane's equation. </summary>
	T Scalar() const {
		return scalar;
	}

	/// <summary> Returns a point on the hyperplane. </summary>
	Vec Base() const {
		return normal * scalar;
	}

	/// <summary> Return the distance of a point from the hyperplane. </summary>
	T Distance(const Vector<T, Dim>& point) const {
		return Dot(point, normal) - scalar;
	}

public:
	Vec normal;
	T scalar;
};


template <class T>
Hyperplane(const Line<T, 2>&) -> Hyperplane<T, 2>;


template <class T>
Line(const Hyperplane<T, 2>&) -> Line<T, 2>;

} // namespace mathter