// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Vector.hpp"
#include "Arithmetic.hpp"
#include "Math.hpp"
#include "Quaternion.hpp"

namespace mathter {

//------------------------------------------------------------------------------
// Quaternion x Vector
//------------------------------------------------------------------------------

template <class T, mathter::eQuaternionLayout Layout, bool Packed>
template <class TOther, bool PackedOther>
auto mathter::Quaternion<T, Layout, Packed>::operator()(const Vector<TOther, 3, PackedOther>& v) const {
	// Sandwich product :)
	const auto vq = Quaternion<T, Layout, Packed>(v);
	const auto rotated = *this * vq * Inverse(*this);
	return Vector(rotated.vector);
}


/// <summary> Rotate the vector by the quaterion. </summary>
/// <remarks> This method is deprecated because its meaning is a bit ambiguous.
///		Use the operator() instead, which is pretty straightforward. </remarks>
template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, bool Packed2>
[[deprecated("use Quaternion::operator()(Vector) instead")]]
auto operator*(const Quaternion<T1, Layout1, Packed1>& q, const Vector<T2, 3, Packed2>& v) {
	return q(v);
}


/// <summary> Rotate the vector by the quaterion. </summary>
/// <remarks> This method is deprecated because its meaning is a bit ambiguous.
///		Use the operator() instead, which is pretty straightforward. </remarks>
template <class T1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
[[deprecated("use Quaternion::operator()(Vector) instead")]]
auto operator*(const Vector<T1, 3, Packed1>& v, const Quaternion<T2, Layout2, Packed2>& q) {
	return q(v);
}


/// <summary> Rotate the vector by the quaterion. </summary>
/// <remarks> This method is deprecated because its meaning is a bit ambiguous.
///		Use the operator() instead, which is pretty straightforward. </remarks>
template <class T1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
[[deprecated("use Quaternion::operator()(Vector) instead")]]
auto& operator*=(Vector<T1, 3, Packed1>& v, const Quaternion<T2, Layout2, Packed2>& q) {
	return v = q(v);
}

} // namespace mathter