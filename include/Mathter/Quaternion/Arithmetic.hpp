// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "Quaternion.hpp"

namespace mathter {


namespace impl {

	template <class T1, eQuaternionLayout Layout1, bool Packed1,
			  class T2, eQuaternionLayout Layout2, bool Packed2>
	auto Multiply(const Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
		const Vector lhsv = lhs.canonical;
		const Vector rhsv = rhs.canonical;
		const Vector<T1, 4, Packed1> s1 = { -1, 1, -1, 1 };
		const Vector<T1, 4, Packed1> s2 = { -1, 1, 1, -1 };
		const Vector<T1, 4, Packed1> s3 = { -1, -1, 1, 1 };

		// [ 3, 2, 1, 0 ]
		// [ 0, 3, 2, 1 ]
		const Vector t0 = lhsv.xxxx;
		const Vector t1 = rhsv.xyzw;

		const Vector t2 = lhsv.yyyy;
		const Vector t3 = rhsv.yxwz;

		const Vector t4 = lhsv.zzzz;
		const Vector t5 = rhsv.zwxy;

		const Vector t6 = lhsv.wwww;
		const Vector t7 = rhsv.wzyx;

		const auto m0 = t0 * t1;
		const auto m1 = s1 * t2 * t3;
		const auto m2 = s2 * t4 * t5;
		const auto m3 = s3 * t6 * t7;

		const auto canonical = m0 + m1 + m2 + m3;
		using Canonical = std::decay_t<decltype(canonical)>;
		return Quaternion < scalar_type_t<Canonical>, Layout1, Packed1 && Packed2 > { canonicalArg, canonical };
	}

} // namespace impl


//------------------------------------------------------------------------------
// Quaternion x Quaternion
//------------------------------------------------------------------------------

template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto operator+(const Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	using R = common_arithmetic_type_t<T1, T2>;
	return Quaternion < R, Layout1, Packed1 && Packed2 > (canonicalArg, lhs.canonical + rhs.canonical);
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto operator-(const Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	using R = common_arithmetic_type_t<T1, T2>;
	return Quaternion < R, Layout1, Packed1 && Packed2 > (canonicalArg, lhs.canonical - rhs.canonical);
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto operator*(const Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return impl::Multiply(lhs, rhs);
}


//------------------------------------------------------------------------------
// Quaternion x Quaternion assign
//------------------------------------------------------------------------------

template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto& operator+=(Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return lhs = lhs + rhs;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto& operator-=(Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return lhs = lhs - rhs;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
auto& operator*=(Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return lhs = lhs * rhs;
}


//------------------------------------------------------------------------------
// Quaternion x Scalar
//------------------------------------------------------------------------------

template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto operator+(const Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	return lhs + Quaternion<T2, Layout1, Packed1>(rhs);
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto operator-(const Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	return lhs - Quaternion<T2, Layout1, Packed1>(rhs);
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto operator*(const Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	using R = common_arithmetic_type_t<T1, T2>;
	return Quaternion<R, Layout1, Packed1>(Vector(lhs) * rhs);
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto operator/(const Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	using R = common_arithmetic_type_t<T1, T2>;
	return Quaternion<R, Layout1, Packed1>(Vector(lhs) / rhs);
}

//------------------------------------------------------------------------------
// Scalar x Quaternion
//------------------------------------------------------------------------------

template <class T1,
		  class T2, eQuaternionLayout Layout2, bool Packed2,
		  class = std::enable_if_t<is_scalar_v<T1>>>
auto operator+(const T1& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return Quaternion<T2, Layout2, Packed2>(lhs) + rhs;
}


template <class T1,
		  class T2, eQuaternionLayout Layout2, bool Packed2,
		  class = std::enable_if_t<is_scalar_v<T1>>>
auto operator-(const T1& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return Quaternion<T2, Layout2, Packed2>(lhs) - rhs;
}


template <class T1,
		  class T2, eQuaternionLayout Layout2, bool Packed2,
		  class = std::enable_if_t<is_scalar_v<T1>>>
auto operator*(const T1& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	using R = common_arithmetic_type_t<T1, T2>;
	return Quaternion<R, Layout2, Packed2>(lhs * Vector(rhs));
}


//------------------------------------------------------------------------------
// Quaternion x Scalar assign
//------------------------------------------------------------------------------

template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto& operator+=(Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	return lhs = lhs + rhs;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto& operator-=(Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	return lhs = lhs - rhs;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto& operator*=(Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	return lhs = lhs * rhs;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2,
		  class = std::enable_if_t<is_scalar_v<T2>>>
auto& operator/=(Quaternion<T1, Layout1, Packed1>& lhs, const T2& rhs) {
	return lhs = lhs / rhs;
}

} // namespace mathter