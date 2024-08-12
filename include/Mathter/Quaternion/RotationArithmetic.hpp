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

template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, bool Packed2>
auto operator*(const Quaternion<T1, Layout1, Packed1>& q, const Vector<T2, 3, Packed2>& v) {
	return Vector((q * Quaternion<T2, Layout1, Packed1>(T2(0), v) * Inverse(q)).vector);
}


template <class T1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto operator*(const Vector<T1, 3, Packed1>& v, const Quaternion<T2, Layout2, Packed2>& q) {
	return Vector((q * Quaternion<T2, Layout2, Packed2>(T2(0), v) * Inverse(q)).vector);
}


template <class T1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto& operator*=(Vector<T1, 3, Packed1>& v, const Quaternion<T2, Layout2, Packed2>& q) {
	return v = v * q;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, bool Packed2>
auto operator*(const Quaternion<T1, Layout1, Packed1>& q, const Vector<T2, 4, Packed2>& v) {
	return Vector((q * Vector(v.xyz)), v.w);
}


template <class T1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto operator*(const Vector<T1, 4, Packed1>& v, const Quaternion<T2, Layout2, Packed2>& q) {
	return Vector((Vector(v.xyz) * q), v.w);
}


template <class T1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto& operator*=(Vector<T1, 4, Packed1>& v, const Quaternion<T2, Layout2, Packed2>& q) {
	return v = v * q;
}

} // namespace mathter