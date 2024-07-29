// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma once

#include "Swizzle.hpp"
#include "Vector.hpp"

#include <algorithm>


namespace mathter {

//--------------------------------------
// Vector x Vector
//--------------------------------------

/// <summary> Exactly compares two vectors. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator==(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs) {
	return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

/// <summary> Exactly compares two vectors. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator!=(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs) {
	return !operator==(lhs, rhs);
}

//--------------------------------------
// Vector x Swizzle
//--------------------------------------

/// <summary> Exactly compares a vector and a swizzle. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Dim2, bool Packed2, int... Indices2>
bool operator==(const Vector<T1, Dim1, Packed1>& lhs, const Swizzle<T2, Dim2, Packed2, Indices2...>& rhs) {
	constexpr auto TargetDim2 = target_dimension_v<std::decay_t<decltype(rhs)>>;
	static_assert(Dim1 == TargetDim2);
	const auto rhsLin = rhs.Linearize();
	return std::equal(lhs.begin(), lhs.end(), rhsLin.begin(), rhsLin.begin() + TargetDim2);
}

/// <summary> Exactly compares a vector and a swizzle. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Dim2, bool Packed2, int... Indices2>
bool operator!=(const Vector<T1, Dim1, Packed1>& lhs, const Swizzle<T2, Dim2, Packed2, Indices2...>& rhs) {
	return !operator==(lhs, rhs);
}

//--------------------------------------
// Swizzle x Vector
//--------------------------------------

/// <summary> Exactly compares a vector and a swizzle. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Dim2, bool Packed2, int... Indices1>
bool operator==(const Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const Vector<T2, Dim2, Packed2>& rhs) {
	constexpr auto TargetDim1 = target_dimension_v<std::decay_t<decltype(lhs)>>;
	static_assert(TargetDim1 == Dim2);
	const auto lhsLin = lhs.Linearize();
	return std::equal(lhsLin.begin(), lhsLin.begin() + TargetDim1, rhs.begin(), rhs.end());
}

/// <summary> Exactly compares a vector and a swizzle. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Dim2, bool Packed2, int... Indices1>
bool operator!=(const Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const Vector<T2, Dim2, Packed2>& rhs) {
	return !operator==(lhs, rhs);
}

//--------------------------------------
// Swizzle x Swizzle
//--------------------------------------

/// <summary> Exactly compares two swizzles. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Dim2, bool Packed2, int... Indices1, int... Indices2>
bool operator==(const Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const Swizzle<T2, Dim2, Packed2, Indices2...>& rhs) {
	constexpr auto TargetDim1 = target_dimension_v<std::decay_t<decltype(lhs)>>;
	constexpr auto TargetDim2 = target_dimension_v<std::decay_t<decltype(rhs)>>;
	static_assert(TargetDim1 == TargetDim2);
	const auto lhsLin = lhs.Linearize();
	const auto rhsLin = rhs.Linearize();
	return std::equal(lhsLin.begin(), lhsLin.begin() + TargetDim1, rhsLin.begin(), rhsLin.begin() + TargetDim2);
}

/// <summary> Exactly compares two swizzles. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Dim2, bool Packed2, int... Indices1, int... Indices2>
bool operator!=(const Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const Swizzle<T2, Dim2, Packed2, Indices2...>& rhs) {
	return !operator==(lhs, rhs);
}

//--------------------------------------
// Swizzle<1> x Scalar
//--------------------------------------

/// <summary> Exactly compares a single-element swizzle with a scalar. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Index>
auto operator==(const Swizzle<T1, Dim1, Packed1, Index>& lhs, const T2& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T2>>, bool> {
	return lhs[0] == rhs;
}

/// <summary> Exactly compares a single-element swizzle with a scalar. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, int Dim1, bool Packed1, class T2, int Index>
auto operator!=(const Swizzle<T1, Dim1, Packed1, Index>& lhs, const T2& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T2>>, bool> {
	return lhs[0] != rhs;
}

template <class T1, int Dim1, bool Packed1, class T2, int Index>
auto operator<(const Swizzle<T1, Dim1, Packed1, Index>& lhs, const T2& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T2>>, bool> {
	return lhs[0] < rhs;
}

template <class T1, int Dim1, bool Packed1, class T2, int Index>
auto operator>(const Swizzle<T1, Dim1, Packed1, Index>& lhs, const T2& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T2>>, bool> {
	return lhs[0] > rhs;
}

template <class T1, int Dim1, bool Packed1, class T2, int Index>
auto operator<=(const Swizzle<T1, Dim1, Packed1, Index>& lhs, const T2& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T2>>, bool> {
	return lhs[0] <= rhs;
}

template <class T1, int Dim1, bool Packed1, class T2, int Index>
auto operator>=(const Swizzle<T1, Dim1, Packed1, Index>& lhs, const T2& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T2>>, bool> {
	return lhs[0] >= rhs;
}

//--------------------------------------
// Scalar x Swizzle<1>
//--------------------------------------

/// <summary> Exactly compares a single-element swizzle with a scalar. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, class T2, int Dim2, bool Packed2, int Index>
auto operator==(const T1& lhs, const Swizzle<T2, Dim2, Packed2, Index>& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T1>>, bool> {
	return lhs == rhs[0];
}

/// <summary> Exactly compares a single-element swizzle with a scalar. </summary>
/// <remarks> &lt;The usual warning about floating point numbers&gt; </remarks>
template <class T1, class T2, int Dim2, bool Packed2, int Index>
auto operator!=(const T1& lhs, const Swizzle<T2, Dim2, Packed2, Index>& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T1>>, bool> {
	return lhs != rhs[0];
}

template <class T1, class T2, int Dim2, bool Packed2, int Index>
auto operator<(const T1& lhs, const Swizzle<T2, Dim2, Packed2, Index>& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T1>>, bool> {
	return lhs < rhs[0];
}

template <class T1, class T2, int Dim2, bool Packed2, int Index>
auto operator>(const T1& lhs, const Swizzle<T2, Dim2, Packed2, Index>& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T1>>, bool> {
	return lhs > rhs[0];
}

template <class T1, class T2, int Dim2, bool Packed2, int Index>
auto operator<=(const T1& lhs, const Swizzle<T2, Dim2, Packed2, Index>& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T1>>, bool> {
	return lhs <= rhs[0];
}

template <class T1, class T2, int Dim2, bool Packed2, int Index>
auto operator>=(const T1& lhs, const Swizzle<T2, Dim2, Packed2, Index>& rhs)
	-> std::enable_if_t<is_scalar_v<std::decay_t<T1>>, bool> {
	return lhs >= rhs[0];
}

} // namespace mathter