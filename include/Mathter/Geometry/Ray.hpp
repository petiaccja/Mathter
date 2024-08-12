// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "Line.hpp"


namespace mathter {

/// <summary> An line with a starting point and an open end. </summary>
template <class T, int Dim>
class Ray : protected Line<T, Dim> {
	// The implementation is the same as for lines, hence the inheritance.
	// However, rays are conceptually different (i.e. for intersections), so
	// the inheritance is protected.
public:
	using Line<T, Dim>::Line;
	using Line<T, Dim>::Through;
	using Line<T, Dim>::Direction;
	using Line<T, Dim>::Base;
	using Line<T, Dim>::PointAt;
	using Line<T, Dim>::direction;
	using Line<T, Dim>::base;

	/// <summary> Converts from a line with different scalar type. </summary>
	template <class TOther>
	Ray(const Ray<TOther, Dim>& other) : Ray(other.base, other.direction) {}

	/// <summary> Returns a line colinear to the ray. </summary>
	mathter::Line<T, Dim> Line() const {
		return static_cast<mathter::Line<T, Dim>>(*this);
	}
};

} // namespace mathter