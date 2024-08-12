// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/Functional.hpp"
#include "OperationUtil.hpp"
#include "Vector.hpp"

#include <functional>


namespace mathter {


template <class Vec, class Fun>
auto AvoidDivByZero(const Vec& vec, const Fun&) {
	using VecDecay = std::decay_t<Vec>;
	if constexpr (VecDecay::isBatched && std::is_same_v<std::decay_t<Fun>, std::divides<void>>) {
		return VecDecay(FillMasked<dimension_v<VecDecay>>(vec.elements.Load(), static_cast<scalar_type_t<VecDecay>>(1)));
	}
	return vec;
}


#define MATHTER_ARITHMETIC_VEC_x_VEC(OP, FUNCTOR)                                                \
	template <class T1, class T2, int Dim, bool Packed1, bool Packed2>                           \
	auto operator OP(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs) { \
		return DoBinaryOp(lhs, AvoidDivByZero(rhs, FUNCTOR{}), FUNCTOR{});                       \
	}


#define MATHTER_ARITHMETIC_VEC_x_SCALAR(OP, FUNCTOR)                                                  \
	template <class T1, int Dim1, bool Packed1, class T2, std::enable_if_t<is_scalar_v<T2>, int> = 0> \
	auto operator OP(const Vector<T1, Dim1, Packed1>& lhs, const T2& rhs) {                           \
		const auto rhsv = Vector<T2, Dim1, Packed1>(rhs);                                             \
		return DoBinaryOp(lhs, rhsv, FUNCTOR{});                                                      \
	}


#define MATHTER_ARITHMETIC_SCALAR_x_VEC(OP, FUNCTOR)                                                  \
	template <class T1, class T2, int Dim2, bool Packed2, std::enable_if_t<is_scalar_v<T1>, int> = 0> \
	auto operator OP(const T1& lhs, const Vector<T2, Dim2, Packed2>& rhs) {                           \
		const auto lhsv = Vector<T1, Dim2, Packed2>(lhs);                                             \
		return DoBinaryOp(lhsv, AvoidDivByZero(rhs, FUNCTOR{}), FUNCTOR{});                           \
	}


#define MATHTER_ARITHMETIC_VEC_x_SWIZZLE(OP, FUNCTOR)                                                            \
	template <class T1, int Dim1, bool Packed1, class T2, int Dim2, bool Packed2, int... Indices2>               \
	auto operator OP(const Vector<T1, Dim1, Packed1>& lhs, const Swizzle<T2, Dim2, Packed2, Indices2...>& rhs) { \
		static_assert(Dim1 == sizeof...(Indices2), "vector and swizzle must have the same size");                \
		const auto rhsv = Vector(rhs);                                                                           \
		return DoBinaryOp(lhs, AvoidDivByZero(rhsv, FUNCTOR{}), FUNCTOR{});                                      \
	}


#define MATHTER_ARITHMETIC_SWIZZLE_x_VEC(OP, FUNCTOR)                                                            \
	template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2>               \
	auto operator OP(const Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const Vector<T2, Dim2, Packed2>& rhs) { \
		static_assert(Dim2 == sizeof...(Indices1), "vector and swizzle must have the same size");                \
		const auto lhsv = Vector(lhs);                                                                           \
		return DoBinaryOp(lhsv, AvoidDivByZero(rhs, FUNCTOR{}), FUNCTOR{});                                      \
	}


#define MATHTER_ARITHMETIC_SWIZZLE_x_SWIZZLE(OP, FUNCTOR)                                                                      \
	template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, int Dim2, bool Packed2, int... Indices2>            \
	auto operator OP(const Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const Swizzle<T2, Dim2, Packed2, Indices2...>& rhs) { \
		static_assert(sizeof...(Indices1) == sizeof...(Indices2), "swizzles must have the same size");                         \
		const auto lhsv = Vector(lhs);                                                                                         \
		const auto rhsv = Vector(rhs);                                                                                         \
		const auto result = DoBinaryOp(lhsv, AvoidDivByZero(rhsv, FUNCTOR{}), FUNCTOR{});                                      \
		if constexpr (dimension_v<std::decay_t<decltype(result)>> == 1) {                                                      \
			return result[0];                                                                                                  \
		}                                                                                                                      \
		else {                                                                                                                 \
			return result;                                                                                                     \
		}                                                                                                                      \
	}


#define MATHTER_ARITHMETIC_SWIZZLE_x_SCALAR(OP, FUNCTOR)                                                               \
	template <class T1, int Dim1, bool Packed1, int... Indices1, class T2, std::enable_if_t<is_scalar_v<T2>, int> = 0> \
	auto operator OP(const Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const T2& rhs) {                              \
		const auto lhsv = Vector(lhs);                                                                                 \
		const auto rhsv = Vector<T2, sizeof...(Indices1), Packed1>(rhs);                                               \
		const auto result = DoBinaryOp(lhsv, rhsv, FUNCTOR{});                                                         \
		if constexpr (dimension_v<std::decay_t<decltype(result)>> == 1) {                                              \
			return result[0];                                                                                          \
		}                                                                                                              \
		else {                                                                                                         \
			return result;                                                                                             \
		}                                                                                                              \
	}


#define MATHTER_ARITHMETIC_SCALAR_x_SWIZZLE(OP, FUNCTOR)                                                               \
	template <class T1, class T2, int Dim2, bool Packed2, int... Indices2, std::enable_if_t<is_scalar_v<T1>, int> = 0> \
	auto operator OP(const T1& lhs, const Swizzle<T2, Dim2, Packed2, Indices2...>& rhs) {                              \
		const auto lhsv = Vector<T1, sizeof...(Indices2), Packed2>(lhs);                                               \
		const auto rhsv = Vector(rhs);                                                                                 \
		const auto result = DoBinaryOp(lhsv, AvoidDivByZero(rhsv, FUNCTOR{}), FUNCTOR{});                              \
		if constexpr (dimension_v<std::decay_t<decltype(result)>> == 1) {                                              \
			return result[0];                                                                                          \
		}                                                                                                              \
		else {                                                                                                         \
			return result;                                                                                             \
		}                                                                                                              \
	}


MATHTER_ARITHMETIC_VEC_x_VEC(*, std::multiplies);
MATHTER_ARITHMETIC_VEC_x_VEC(/, std::divides);
MATHTER_ARITHMETIC_VEC_x_VEC(+, std::plus);
MATHTER_ARITHMETIC_VEC_x_VEC(-, std::minus);

MATHTER_ARITHMETIC_VEC_x_SCALAR(*, std::multiplies);
MATHTER_ARITHMETIC_VEC_x_SCALAR(/, std::divides);
MATHTER_ARITHMETIC_VEC_x_SCALAR(+, std::plus);
MATHTER_ARITHMETIC_VEC_x_SCALAR(-, std::minus);

MATHTER_ARITHMETIC_SCALAR_x_VEC(*, std::multiplies);
MATHTER_ARITHMETIC_SCALAR_x_VEC(/, std::divides);
MATHTER_ARITHMETIC_SCALAR_x_VEC(+, std::plus);
MATHTER_ARITHMETIC_SCALAR_x_VEC(-, std::minus);

MATHTER_ARITHMETIC_VEC_x_SWIZZLE(*, std::multiplies);
MATHTER_ARITHMETIC_VEC_x_SWIZZLE(/, std::divides);
MATHTER_ARITHMETIC_VEC_x_SWIZZLE(+, std::plus);
MATHTER_ARITHMETIC_VEC_x_SWIZZLE(-, std::minus);

MATHTER_ARITHMETIC_SWIZZLE_x_VEC(*, std::multiplies);
MATHTER_ARITHMETIC_SWIZZLE_x_VEC(/, std::divides);
MATHTER_ARITHMETIC_SWIZZLE_x_VEC(+, std::plus);
MATHTER_ARITHMETIC_SWIZZLE_x_VEC(-, std::minus);

MATHTER_ARITHMETIC_SWIZZLE_x_SWIZZLE(*, std::multiplies);
MATHTER_ARITHMETIC_SWIZZLE_x_SWIZZLE(/, std::divides);
MATHTER_ARITHMETIC_SWIZZLE_x_SWIZZLE(+, std::plus);
MATHTER_ARITHMETIC_SWIZZLE_x_SWIZZLE(-, std::minus);

MATHTER_ARITHMETIC_SWIZZLE_x_SCALAR(*, std::multiplies);
MATHTER_ARITHMETIC_SWIZZLE_x_SCALAR(/, std::divides);
MATHTER_ARITHMETIC_SWIZZLE_x_SCALAR(+, std::plus);
MATHTER_ARITHMETIC_SWIZZLE_x_SCALAR(-, std::minus);

MATHTER_ARITHMETIC_SCALAR_x_SWIZZLE(*, std::multiplies);
MATHTER_ARITHMETIC_SCALAR_x_SWIZZLE(/, std::divides);
MATHTER_ARITHMETIC_SCALAR_x_SWIZZLE(+, std::plus);
MATHTER_ARITHMETIC_SCALAR_x_SWIZZLE(-, std::minus);


#define MATHTER_ARITHMETIC_VECTOR_ASSIGN(OP)                                  \
	template <class T1, int Dim1, bool Packed1, class T2Any>                  \
	auto& operator OP##=(Vector<T1, Dim1, Packed1>& lhs, const T2Any & rhs) { \
		return lhs = lhs OP rhs;                                              \
	}


#define MATHTER_ARITHMETIC_SWIZZLE_ASSIGN(OP)                                               \
	template <class T1, int Dim1, bool Packed1, int... Indices1, class T2Any>               \
	auto& operator OP##=(Swizzle<T1, Dim1, Packed1, Indices1...>& lhs, const T2Any & rhs) { \
		return lhs = lhs OP rhs;                                                            \
	}


MATHTER_ARITHMETIC_VECTOR_ASSIGN(*);
MATHTER_ARITHMETIC_VECTOR_ASSIGN(/);
MATHTER_ARITHMETIC_VECTOR_ASSIGN(+);
MATHTER_ARITHMETIC_VECTOR_ASSIGN(-);

MATHTER_ARITHMETIC_SWIZZLE_ASSIGN(*);
MATHTER_ARITHMETIC_SWIZZLE_ASSIGN(/);
MATHTER_ARITHMETIC_SWIZZLE_ASSIGN(+);
MATHTER_ARITHMETIC_SWIZZLE_ASSIGN(-);


template <class T1, int Dim1, bool Packed1>
const auto& operator+(const Vector<T1, Dim1, Packed1>& arg) {
	return arg;
}

template <class T1, int Dim1, bool Packed1, int... Indices1>
const auto& operator+(const Swizzle<T1, Dim1, Packed1, Indices1...>& arg) {
	return arg;
}

template <class T1, int Dim1, bool Packed1>
auto operator-(const Vector<T1, Dim1, Packed1>& arg) {
	return arg * -static_cast<T1>(1);
}

template <class T1, int Dim1, bool Packed1, int... Indices1>
auto operator-(const Swizzle<T1, Dim1, Packed1, Indices1...>& arg) {
	return arg * -static_cast<T1>(1);
}


template <class T1, bool Packed1,
		  class T2, bool Packed2,
		  class T3, bool Packed3,
		  int Dim>
auto MultiplyAdd(
	const Vector<T1, Dim, Packed1>& a,
	const Vector<T2, Dim, Packed2>& b,
	const Vector<T3, Dim, Packed3>& c) {
	return DoTernaryOp(a, b, c, ::mathter::fma{});
}

} // namespace mathter