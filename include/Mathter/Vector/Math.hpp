// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#if _MSC_VER && defined(min)
#pragma push_macro("min")
#pragma push_macro("max")
#undef min
#undef max
#define MATHTER_MINMAX
#endif


#include "../Common/Functional.hpp"
#include "../Common/MathUtil.hpp"
#include "Arithmetic.hpp"
#include "OperationUtil.hpp"
#include "Vector.hpp"

#include <numeric>


namespace mathter {


/// <summary> Returns the element-wise minimum of arguments </summary>
template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
auto Min(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs) {
	return DoBinaryOp(lhs, rhs, min{});
}


/// <summary> Returns the element-wise maximum of arguments </summary>
template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
auto Max(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs) {
	return DoBinaryOp(lhs, rhs, max{});
}


/// <summary> Returns the minimum element of the vector. </summary>
template <class T, int Dim, bool Packed>
T Min(const Vector<T, Dim, Packed>& v) {
#if MATHTER_ENABLE_SIMD
	if constexpr (v.isBatched) {
		const auto value = v.elements.Load();
		constexpr auto filler = std::numeric_limits<remove_complex_t<T>>::max();
		const auto filled = FillMasked<Dim>(value, filler);
		return xsimd::reduce_min(filled);
	}
#endif
	return *std::min_element(v.begin(), v.end());
}


/// <summary> Returns the maximum element of the vector. </summary>
template <class T, int Dim, bool Packed>
T Max(const Vector<T, Dim, Packed>& v) {
#if MATHTER_ENABLE_SIMD
	if constexpr (IsBatched<T, Dim, Packed>()) {
		const auto value = v.elements.Load();
		constexpr auto filler = std::numeric_limits<remove_complex_t<T>>::lowest();
		const auto filled = FillMasked<Dim>(value, filler);
		return xsimd::reduce_max(filled);
	}
#endif
	return *std::max_element(v.begin(), v.end());
}


/// <summary> Returns the elementwise absolute value of the vector. </summary>
template <class T, int Dim, bool Packed>
auto Abs(const Vector<T, Dim, Packed>& v) {
	return DoUnaryOp(v, abs{});
}



/// <summary> Returns the elementwise absolute value of the vector. </summary>
template <class T, int Dim, bool Packed>
auto Conj(const Vector<T, Dim, Packed>& v) {
	return DoUnaryOp(v, conj{});
}


/// <summary> Returns the sum of the element of the vector. </summary>
template <class T, int Dim, bool Packed>
T Sum(const Vector<T, Dim, Packed>& v) {
#if MATHTER_ENABLE_SIMD
	if constexpr (v.isBatched) {
		const auto value = v.elements.Load();
		constexpr auto filler = T(0);
		const auto filled = FillMasked<Dim>(value, filler);
		return xsimd::reduce_add(filled);
	}
#endif
	return std::reduce(v.begin(), v.end());
}


/// <summary> Calculates the scalar product (dot product) of the two arguments. </summary>
template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
auto Dot(const Vector<T1, Dim, Packed1>& lhs, const Vector<T2, Dim, Packed2>& rhs) {
	return Sum(lhs * Conj(rhs));
}


/// <summary> Returns the squared length of the vector. </summary>
template <class T, int Dim, bool Packed>
auto LengthSquared(const Vector<T, Dim, Packed>& v) {
	return std::real(Dot(v, v));
}


/// <summary> Returns the length of the vector. </summary>
template <class T, int Dim, bool Packed>
auto Length(const Vector<T, Dim, Packed>& v) {
	return std::sqrt(LengthSquared(v));
}


/// <summary> Returns the length of the vector. </summary>
/// <remarks> Avoids overflow and underflow, so it's more expensive. </remarks>
template <class T, int Dim, bool Packed>
auto LengthPrecise(const Vector<T, Dim, Packed>& v) {
	const auto maxElement = Max(Abs(v));
	using Base = remove_complex_t<T>;
	if (maxElement == std::numeric_limits<Base>::infinity()) {
		return std::numeric_limits<Base>::infinity();
	}
	if (maxElement == T(0)) {
		return static_cast<Base>(0);
	}
	const auto scaled = v / maxElement;
	return static_cast<Base>(Length(scaled) * maxElement);
}


/// <summary> Returns the euclidean distance between to vectors. </summary>
template <class T, class U, int Dim, bool Packed1, bool Packed2>
auto Distance(const Vector<T, Dim, Packed1>& lhs, const Vector<U, Dim, Packed2>& rhs) {
	return Length(lhs - rhs);
}


/// <summary> Returns the euclidean distance between to vectors. </summary>
/// <remarks> Avoids overflow and underflow, so it's more expensive. </remarks>
template <class T, class U, int Dim, bool Packed1, bool Packed2>
auto DistancePrecise(const Vector<T, Dim, Packed1>& lhs, const Vector<U, Dim, Packed2>& rhs) {
	return LengthPrecise(lhs - rhs);
}


/// <summary> Makes a unit vector, but keeps direction. </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Normalize(const Vector<T, Dim, Packed>& v) {
	const auto length = Length(v);
	const auto zero = static_cast<std::decay_t<decltype(length)>>(0);
	assert(length != zero);
	return v / length;
}


/// <summary> Makes a unit vector, but keeps direction. </summary>
/// <param name="degenerate"> Returned if <paramref name="v"/> is a null vector. Should be a unit vector. </param>
/// <remarks> Unlike the regular <see cref="Normalize"/>, this does can handle null vectors and under/overflow. </remarks>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> NormalizePrecise(const Vector<T, Dim, Packed>& v, const Vector<T, Dim, Packed>& degenerate) {
	const auto length = LengthPrecise(v);
	const auto zero = static_cast<std::decay_t<decltype(length)>>(0);
	return length != zero ? v / length : degenerate;
}


/// <summary> Makes a unit vector, but keeps direction. </summary>
/// <remarks> Unlike the regular <see cref="Normalize"/>, this does can handle null vectors and under/overflow. </remarks>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> NormalizePrecise(const Vector<T, Dim, Packed>& v) {
	Vector<T, Dim, Packed> degenerate(static_cast<T>(0));
	degenerate[0] = static_cast<T>(1);
	return NormalizePrecise(v, degenerate);
}


/// <summary> Sets all elements of the vector to the same value. </summary>
template <class T, int Dim, bool Packed, class U, std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
void Fill(Vector<T, Dim, Packed>& lhs, U&& all) {
	lhs = Vector<T, Dim, Packed>(static_cast<T>(std::forward<U>(all)));
}


/// <summary> Returns the generalized cross-product in N dimensions. </summary>
/// <remarks> You must supply N-1 arguments of type Vector&lt;N&gt;.
/// The function returns the generalized cross product as defined by
/// https://en.wikipedia.org/wiki/Cross_product#Multilinear_algebra. </remarks>
template <class T, int Dim, bool Packed, class... Args>
auto Cross(const Vector<T, Dim, Packed>& head, Args&&... args) -> Vector<T, Dim, Packed>;


/// <summary> Returns the generalized cross-product in N dimensions. </summary>
/// <remarks> See https://en.wikipedia.org/wiki/Cross_product#Multilinear_algebra for definition. </remarks>
template <class T, int Dim, bool Packed>
auto Cross(const std::array<const Vector<T, Dim, Packed>*, Dim - 1>& args) -> Vector<T, Dim, Packed>;

/// <summary> Returns the 2-dimensional cross product, which is a vector perpendicular to the argument. </summary>
template <class T, bool Packed>
Vector<T, 2, Packed> Cross(const Vector<T, 2, Packed>& arg) {
	return Vector<T, 2, Packed>(-arg.y,
								arg.x);
}
/// <summary> Returns the 2-dimensional cross product, which is a vector perpendicular to the argument. </summary>
template <class T, bool Packed>
Vector<T, 2, Packed> Cross(const std::array<const Vector<T, 2, Packed>*, 1>& arg) {
	return Cross(*(arg[0]));
}


/// <summary> Returns the 3-dimensional cross-product. </summary>
template <class T, bool Packed>
Vector<T, 3, Packed> Cross(const Vector<T, 3, Packed>& lhs, const Vector<T, 3, Packed>& rhs) {
	return lhs.yzx * rhs.zxy - lhs.zxy * rhs.yzx;
}


/// <summary> Returns the 3-dimensional cross-product. </summary>
template <class T, bool Packed>
Vector<T, 3, Packed> Cross(const std::array<const Vector<T, 3, Packed>*, 2>& args) {
	return Cross(*(args[0]), *(args[1]));
}

} // namespace mathter


/*
// Generalized cross-product unfortunately needs matrix determinant.
#include "../Matrix/MatrixFunction.hpp"

namespace mathter {

template <class T, int Dim, bool Packed>
auto Cross(const std::array<const Vector<T, Dim, Packed>*, Dim - 1>& args) -> Vector<T, Dim, Packed> {
	Vector<T, Dim, Packed> result;
	Matrix<T, Dim - 1, Dim - 1, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false> detCalc;

	// Calculate elements of result on-by-one
	int sign = 2 * (Dim % 2) - 1;
	for (int base = 0; base < result.Dimension(); ++base, sign *= -1) {
		// Fill up sub-matrix the determinant of which yields the coefficient of base-vector.
		for (int j = 0; j < base; ++j) {
			for (int i = 0; i < detCalc.RowCount(); ++i) {
				detCalc(i, j) = (*(args[i]))[j];
			}
		}
		for (int j = base + 1; j < result.Dimension(); ++j) {
			for (int i = 0; i < detCalc.RowCount(); ++i) {
				detCalc(i, j - 1) = (*(args[i]))[j];
			}
		}

		T coefficient = T(sign) * Determinant(detCalc);
		result(base) = coefficient;
	}

	return result;
}


template <class T, int Dim, bool Packed, class... Args>
auto Cross(const Vector<T, Dim, Packed>& head, Args&&... args) -> Vector<T, Dim, Packed> {
	static_assert(1 + sizeof...(args) == Dim - 1, "Number of arguments must be (Dimension - 1).");

	std::array<const Vector<T, Dim, Packed>*, Dim - 1> vectors = { &head, &args... };
	return Cross(vectors);
}

} // namespace mathter
*/


#if defined(MATHTER_MINMAX)
#pragma pop_macro("min")
#pragma pop_macro("max")
#endif