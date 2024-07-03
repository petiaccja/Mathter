// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/Functional.hpp"
#include "VectorImpl.hpp"
#include "VectorUtils.hpp"

#include <functional>


namespace mathter {

//------------------------------------------------------------------------------
// Vector arithmetic
//------------------------------------------------------------------------------

/// <summary> Elementwise (Hadamard) vector product. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator*(const Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return DoBinaryOp(lhs, rhs, std::multiplies{});
}

/// <summary> Elementwise vector division. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator/(const Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return DoBinaryOp(lhs, rhs, std::divides{});
}

/// <summary> Elementwise vector addition. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator+(const Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return DoBinaryOp(lhs, rhs, std::plus{});
}

/// <summary> Elementwise vector subtraction. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator-(const Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return DoBinaryOp(lhs, rhs, std::minus{});
}

//------------------------------------------------------------------------------
// Vector assign arithmetic
//------------------------------------------------------------------------------

/// <summary> Elementwise (Hadamard) vector product. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator*=(Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return lhs = lhs * Vector<T1, Dim, Packed>(rhs);
}

/// <summary> Elementwise vector division. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator/=(Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return lhs = lhs / Vector<T1, Dim, Packed>(rhs);
}

/// <summary> Elementwise vector addition. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator+=(Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return lhs = lhs + Vector<T1, Dim, Packed>(rhs);
}

/// <summary> Elementwise vector subtraction. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator-=(Vector<T1, Dim, Packed>& lhs, const Vector<T2, Dim, Packed>& rhs) {
	return lhs = lhs - Vector<T1, Dim, Packed>(rhs);
}


//------------------------------------------------------------------------------
// Scalar arithmetic
//------------------------------------------------------------------------------

/// <summary> Scales the vector by <paramref name="rhs"/>. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator*(const Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs * Vector<T2, Dim, Packed>(rhs);
}

/// <summary> Scales the vector by 1/<paramref name="rhs"/>. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator/(const Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs / Vector<T2, Dim, Packed>(rhs);
}

/// <summary> Adds <paramref name="rhs"/> to each element of the vector. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator+(const Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs + Vector<T2, Dim, Packed>(rhs);
}

/// <summary> Subtracts <paramref name="rhs"/> from each element of the vector. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator-(const Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs - Vector<T2, Dim, Packed>(rhs);
}


/// <summary> Scales vector by <paramref name="lhs"/>. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator*(T1 lhs, const Vector<T2, Dim, Packed>& rhs) {
	return Vector<T2, Dim, Packed>(lhs) * rhs;
}

/// <summary> Adds <paramref name="lhs"/> to all elements of the vector. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator+(T1 lhs, const Vector<T2, Dim, Packed>& rhs) {
	return Vector<T2, Dim, Packed>(lhs) + rhs;
}

/// <summary> Makes a vector with <paramref name="lhs"/> as all elements, then subtracts <paramref name="rhs"> from it. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator-(T1 lhs, const Vector<T2, Dim, Packed>& rhs) {
	return Vector<T2, Dim, Packed>(lhs) - rhs;
}

/// <summary> Makes a vector with <paramref name="lhs"/> as all elements, then divides it by <paramref name="rhs">. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto operator/(T1 lhs, const Vector<T2, Dim, Packed>& rhs) {
	return Vector<T2, Dim, Packed>(lhs) / rhs;
}


//------------------------------------------------------------------------------
// Scalar assign arithmetic
//------------------------------------------------------------------------------

/// <summary> Scales the vector by <paramref name="rhs"/>. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator*=(Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs = lhs * Vector<T1, Dim, Packed>(static_cast<T1>(rhs));
}

/// <summary> Scales the vector by 1/<paramref name="rhs"/>. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator/=(Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs = lhs / Vector<T1, Dim, Packed>(static_cast<T1>(rhs));
}

/// <summary> Adds <paramref name="rhs"/> to each element of the vector. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator+=(Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs = lhs + Vector<T1, Dim, Packed>(static_cast<T1>(rhs));
}

/// <summary> Subtracts <paramref name="rhs"/> from each element of the vector. </summary>
template <class T1, class T2, int Dim, bool Packed>
auto& operator-=(Vector<T1, Dim, Packed>& lhs, T2 rhs) {
	return lhs = lhs - Vector<T1, Dim, Packed>(static_cast<T1>(rhs));
}


//------------------------------------------------------------------------------
// Extra
//------------------------------------------------------------------------------

/// <summary> Return (a*b)+c. Performs MAD or FMA if supported by target architecture. </summary>
template <class T1, class T2, class T3, int Dim, bool Packed>
auto MultiplyAdd(const Vector<T1, Dim, Packed>& a, const Vector<T2, Dim, Packed>& b, const Vector<T3, Dim, Packed>& c) {
	return DoTernaryOp(a, b, c, mathter::fma{});
}

/// <summary> Negates all elements of the vector. </summary>
template <class T, int Dim, bool Packed>
auto operator-(const Vector<T, Dim, Packed>& arg) {
	return arg * T(-1);
}

/// <summary> Optional plus sign, leaves the vector as is. </summary>
template <class T, int Dim, bool Packed>
auto operator+(const Vector<T, Dim, Packed>& arg) {
	return arg;
}

//------------------------------------------------------------------------------
// Swizzle-vector
//------------------------------------------------------------------------------

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator*(const Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v * Vector<T2, Dim, SwPacked>(s);
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator/(const Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v / Vector<T2, Dim, SwPacked>(s);
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator+(const Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v + Vector<T2, Dim, SwPacked>(s);
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator-(const Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v - Vector<T2, Dim, SwPacked>(s);
}


template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator*(const Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return Vector<T2, Dim, SwPacked>(s) * v;
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator/(const Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return Vector<T2, Dim, SwPacked>(s) / v;
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator+(const Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return Vector<T2, Dim, SwPacked>(s) + v;
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator-(const Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return Vector<T2, Dim, SwPacked>(s) - v;
}


template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator*=(Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v = v * Vector<T2, Dim, SwPacked>(s);
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator/=(Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v = v / Vector<T2, Dim, SwPacked>(s);
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator+=(Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v = v + Vector<T2, Dim, SwPacked>(s);
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto operator-=(Vector<T1, Dim, Packed>& v, const Swizzle<T2, SwDim, SwPacked, Indices...>& s) {
	return v = v - Vector<T2, Dim, SwPacked>(s);
}


template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto& operator*=(Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return s = s * v;
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto& operator/=(Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return s = s / v;
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto& operator+=(Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return s = s + v;
}

template <class T1, class T2, int Dim, bool Packed, int SwDim, bool SwPacked, int... Indices>
auto& operator-=(Swizzle<T2, SwDim, SwPacked, Indices...>& s, const Vector<T1, Dim, Packed>& v) {
	return s = s - v;
}


//------------------------------------------------------------------------------
// Swizzle-swizzle
//------------------------------------------------------------------------------


template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto operator*(const Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	static_assert(sizeof...(Indices1) == sizeof...(Indices2));
	constexpr bool Packed = Packed1 && Packed2;
	using V1 = Vector<T1, sizeof...(Indices1), Packed>;
	using V2 = Vector<T2, sizeof...(Indices2), Packed>;
	return V1(s1) * V2(s2);
}

template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto operator/(const Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	static_assert(sizeof...(Indices1) == sizeof...(Indices2));
	constexpr bool Packed = Packed1 && Packed2;
	using V1 = Vector<T1, sizeof...(Indices1), Packed>;
	using V2 = Vector<T2, sizeof...(Indices2), Packed>;
	return V1(s1) / V2(s2);
}

template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto operator+(const Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	static_assert(sizeof...(Indices1) == sizeof...(Indices2));
	constexpr bool Packed = Packed1 && Packed2;
	using V1 = Vector<T1, sizeof...(Indices1), Packed>;
	using V2 = Vector<T2, sizeof...(Indices2), Packed>;
	return V1(s1) + V2(s2);
}

template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto operator-(const Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	static_assert(sizeof...(Indices1) == sizeof...(Indices2));
	constexpr bool Packed = Packed1 && Packed2;
	using V1 = Vector<T1, sizeof...(Indices1), Packed>;
	using V2 = Vector<T2, sizeof...(Indices2), Packed>;
	return V1(s1) - V2(s2);
}


template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto& operator*=(Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	return s1 = s1 * s2;
}

template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto& operator/=(Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	return s1 = s1 / s2;
}

template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto& operator+=(Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	return s1 = s1 + s2;
}

template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>
auto& operator-=(Swizzle<T1, Dim1, Packed1, Indices1...>& s1, const Swizzle<T2, Dim2, Packed2, Indices2...>& s2) {
	return s1 = s1 - s2;
}

//------------------------------------------------------------------------------
// Swizzle-scalar
//------------------------------------------------------------------------------

template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator*(const Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	return VectorT(lhs) * rhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator/(const Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	return VectorT(lhs) / rhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator+(const Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	return VectorT(lhs) + rhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator-(const Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	return VectorT(lhs) - rhs;
}



template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator*(U lhs, const Swizzle<T, Dim, Packed, Indices...>& rhs) {
	return rhs * lhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator/(U lhs, const Swizzle<T, Dim, Packed, Indices...>& rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	return lhs / VectorT(rhs);
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator+(U lhs, const Swizzle<T, Dim, Packed, Indices...>& rhs) {
	return rhs + lhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto operator-(U lhs, const Swizzle<T, Dim, Packed, Indices...>& rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	return lhs - VectorT(rhs);
}



template <class T, int Dim, bool Packed, int... Indices, class U>
auto& operator*=(Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	lhs = VectorT(lhs) * rhs;
	return lhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto& operator/=(Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	lhs = VectorT(lhs) / rhs;
	return lhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto& operator+=(Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	lhs = VectorT(lhs) + rhs;
	return lhs;
}

template <class T, int Dim, bool Packed, int... Indices, class U>
auto& operator-=(Swizzle<T, Dim, Packed, Indices...>& lhs, U rhs) {
	using VectorT = Vector<T, sizeof...(Indices), Packed>;
	lhs = VectorT(lhs) - rhs;
	return lhs;
}



} // namespace mathter