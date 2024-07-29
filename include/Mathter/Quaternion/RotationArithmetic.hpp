// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Vector.hpp"
#include "Arithmetic.hpp"
#include "Quaternion.hpp"
#include "Math.hpp"

namespace mathter {

//------------------------------------------------------------------------------
// Quaternion x Vector
//------------------------------------------------------------------------------

template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, bool Packed2>
auto operator*(const Quaternion<T1, Layout1, Packed1>& q, const Vector<T2, 3, Packed2>& v) {
	return q * Quaternion<T2, Layout1, Packed1>(T2(0), v) * Inverse(q);
}


template <class T2, bool Packed2,
		  class T1, eQuaternionLayout Layout1, bool Packed1>
auto operator*(const Vector<T2, 3, Packed2>& v, const Quaternion<T1, Layout1, Packed1>& q) {
	return q * Quaternion<T2, Layout1, Packed1>(T2(0), v) * Inverse(q);
}


template <class T2, bool Packed2,
		  class T1, eQuaternionLayout Layout1, bool Packed1>
auto& operator*=(Vector<T2, 3, Packed2>& v, const Quaternion<T1, Layout1, Packed1>& q) {
	return v = v * q;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, bool Packed2>
auto operator*(const Quaternion<T1, Layout1, Packed1>& q, const Vector<T2, 4, Packed2>& v) {
	return (q * Vector(v.xyz)) | v.w;
}


template <class T2, bool Packed2,
		  class T1, eQuaternionLayout Layout1, bool Packed1>
auto operator*(const Vector<T2, 4, Packed2>& v, const Quaternion<T1, Layout1, Packed1>& q) {
	return (Vector(v.xyz) * q) | v.w;
}


template <class T2, bool Packed2,
		  class T1, eQuaternionLayout Layout1, bool Packed1>
auto& operator*=(Vector<T2, 4, Packed2>& v, const Quaternion<T1, Layout1, Packed1>& q) {
	return v = v * q;
}

} // namespace mathter