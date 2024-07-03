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
#include "VectorImpl.hpp"
#include "VectorUtils.hpp"

#include <numeric>


namespace mathter {


/// <summary> Returns the element-wise minimum of arguments </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Min(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	return DoBinaryOp(lhs, rhs, min{});
}


/// <summary> Returns the element-wise maximum of arguments </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Max(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	return DoBinaryOp(lhs, rhs, max{});
}


/// <summary> Returns the minimum element of the vector. </summary>
template <class T, int Dim, bool Packed>
T Min(const Vector<T, Dim, Packed>& v) {
	if constexpr (IsBatched<T, Dim, Packed>()) {
		using MyBatch = Batch<T, Dim, Packed>;
		const auto value = MyBatch::load_unaligned(v.data());
		constexpr auto filler = std::numeric_limits<remove_complex_t<T>>::max();
		const auto filled = FillMasked<Dim>(value, filler);
		return xsimd::reduce_min(filled);
	}
	else {
		return *std::min_element(v.begin(), v.end());
	}
}


/// <summary> Returns the maximum element of the vector. </summary>
template <class T, int Dim, bool Packed>
T Max(const Vector<T, Dim, Packed>& v) {
	if constexpr (IsBatched<T, Dim, Packed>()) {
		using MyBatch = Batch<T, Dim, Packed>;
		const auto value = MyBatch::load_unaligned(v.data());
		constexpr auto filler = std::numeric_limits<remove_complex_t<T>>::lowest();
		const auto filled = FillMasked<Dim>(value, filler);
		return xsimd::reduce_max(filled);
	}
	else {
		return *std::max_element(v.begin(), v.end());
	}
}


/// <summary> Returns the maximum element of the vector. </summary>
template <class T, int Dim, bool Packed>
auto Abs(const Vector<T, Dim, Packed>& v) {
	return DoUnaryOp(v, abs{});
}


/// <summary> Returns the sum of the element of the vector. </summary>
template <class T, int Dim, bool Packed>
T Sum(const Vector<T, Dim, Packed>& v) {
	if constexpr (IsBatched<T, Dim, Packed>()) {
		using MyBatch = Batch<T, Dim, Packed>;
		const auto value = MyBatch::load_unaligned(v.data());
		constexpr auto filler = T(0);
		const auto filled = FillMasked<Dim>(value, filler);
		return xsimd::reduce_add(filled);
	}
	else {
		return std::reduce(v.begin(), v.end());
	}
}


/// <summary> Calculates the scalar product (dot product) of the two arguments. </summary>
template <class T, int Dim, bool Packed>
T Dot(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	return Sum(lhs * rhs);
}


/// <summary> Returns the squared length of the vector. </summary>
template <class T, int Dim, bool Packed>
T LengthSquared(const Vector<T, Dim, Packed>& v) {
	return Dot(v, v);
}

/// <summary> Returns the length of the vector. </summary>
template <class T, int Dim, bool Packed>
T Length(const Vector<T, Dim, Packed>& v) {
	return static_cast<T>(std::sqrt(LengthSquared(v)));
}

/// <summary> Returns the length of the vector, avoids overflow and underflow, so it's more expensive. </summary>
template <class T, int Dim, bool Packed>
T LengthPrecise(const Vector<T, Dim, Packed>& v) {
	const auto maxElement = Max(Abs(v));
	if (maxElement == std::numeric_limits<remove_complex_t<T>>::infinity()) {
		return T(std::numeric_limits<remove_complex_t<T>>::infinity());
	}
	if (maxElement == T(0)) {
		return T(0);
	}
	const auto scaled = v / maxElement;
	return Length(scaled) * maxElement;
}


/// <summary> Returns true if the vector's length is too small for precise calculations (i.e. normalization). </summary>
/// <remarks> "Too small" means smaller than the square root of the smallest number representable by the underlying scalar.
///			This value is ~10^-18 for floats and ~10^-154 for doubles. </remarks>
template <class T, int Dim, bool Packed>
bool IsNullvector(const Vector<T, Dim, Packed>& v) {
	static constexpr T epsilon = T(1) / impl::ConstexprExp10<T>(impl::ConstexprAbs(std::numeric_limits<T>::min_exponent10) / 2);
	T length = Length(v);
	return length < epsilon;
}


/// <summary> Returns the euclidean distance between to vectors. </summary>
template <class T, class U, int Dim, bool Packed1, bool Packed2>
auto Distance(const Vector<T, Dim, Packed1>& lhs, const Vector<U, Dim, Packed2>& rhs) {
	return Length(lhs - rhs);
}

/// <summary> Makes a unit vector, but keeps direction. </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Normalize(const Vector<T, Dim, Packed>& v) {
	assert(!IsNullvector(v));
	T l = Length(v);
	return v / l;
}

/// <summary> Checks if the vector is unit vector. There's some tolerance due to floating points. </summary>
template <class T, int Dim, bool Packed>
bool IsNormalized(const Vector<T, Dim, Packed>& v) {
	T n = LengthSquared(v);
	return T(0.9999) <= n && n <= T(1.0001);
}

/// <summary> Makes a unit vector, but keeps direction. Leans towards (1,0,0...) for nullvectors, costs more. </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> SafeNormalize(const Vector<T, Dim, Packed>& v) {
	Vector<T, Dim, Packed> vmod = v;
	vmod(0) = std::abs(v(0)) > std::numeric_limits<T>::denorm_min() ? v(0) : std::numeric_limits<T>::denorm_min();
	T l = LengthPrecise(vmod);
	return vmod / l;
}

/// <summary> Makes a unit vector, but keeps direction. Leans towards <paramref name="degenerate"/> for nullvectors, costs more. </summary>
/// <param name="degenerate"> Must be a unit vector. </param>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> SafeNormalize(const Vector<T, Dim, Packed>& v, const Vector<T, Dim, Packed>& degenerate) {
	assert(IsNormalized(degenerate));
	T length = LengthPrecise(v);
	if (length == 0) {
		return degenerate;
	}
	return v / length;
}

/// <summary> Sets all elements of the vector to the same value. </summary>
template <class T, int Dim, bool Packed, class U, std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
void Fill(Vector<T, Dim, Packed>& lhs, U all) {
	lhs = Vector<T, Dim, Packed>((T)all);
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

/// <summary> Returns the 2-dimensional cross prodct, which is a vector perpendicular to the argument. </summary>
template <class T, bool Packed>
Vector<T, 2, Packed> Cross(const Vector<T, 2, Packed>& arg) {
	return Vector<T, 2, Packed>(-arg.y,
								arg.x);
}
/// <summary> Returns the 2-dimensional cross prodct, which is a vector perpendicular to the argument. </summary>
template <class T, bool Packed>
Vector<T, 2, Packed> Cross(const std::array<const Vector<T, 2, Packed>*, 1>& arg) {
	return Cross(*(arg[0]));
}


/// <summary> Returns the 3-dimensional cross-product. </summary>
template <class T, bool Packed>
Vector<T, 3, Packed> Cross(const Vector<T, 3, Packed>& lhs, const Vector<T, 3, Packed>& rhs) {
	using VecT = Vector<T, 3, Packed>;
	return VecT(lhs.yzx) * VecT(rhs.zxy) - VecT(lhs.zxy) * VecT(rhs.yzx);
}


/// <summary> Returns the 3-dimensional cross-product. </summary>
template <class T, bool Packed>
Vector<T, 3, Packed> Cross(const std::array<const Vector<T, 3, Packed>*, 2>& args) {
	return Cross(*(args[0]), *(args[1]));
}

} // namespace mathter



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


#if defined(MATHTER_MINMAX)
#pragma pop_macro("min")
#pragma pop_macro("max")
#endif