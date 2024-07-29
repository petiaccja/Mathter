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

#include <array>
#include <numeric>
#include <optional>


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
auto Real(const Vector<T, Dim, Packed>& v) {
	return DoUnaryOp(v, real{});
}


/// <summary> Returns the elementwise absolute value of the vector. </summary>
template <class T, int Dim, bool Packed>
auto Imag(const Vector<T, Dim, Packed>& v) {
	return DoUnaryOp(v, imag{});
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


namespace impl {

	template <class IterFirst, class IterLast, class Vec = typename std::iterator_traits<IterFirst>::value_type>
	auto Cross(IterFirst first, IterLast last) -> std::enable_if_t<is_vector_v<Vec>, Vec>;

} // namespace impl


/// <summary> Returns the generalized cross-product in N dimensions. </summary>
/// <remarks> See https://en.wikipedia.org/wiki/Cross_product#Multilinear_algebra for definition. </remarks>
template <class... Vectors, class Vec = common_arithmetic_type_t<std::decay_t<Vectors>...>>
auto Cross(const Vectors&... vectors) -> std::enable_if_t<(... && is_vector_v<Vectors>), Vec> {
	static_assert(sizeof...(Vectors) == dimension_v<Vec> - 1, "ND cross product needs exactly N-1 vectors.");
	const auto container = std::initializer_list<Vec>{ vectors... };
	return impl::Cross(container.begin(), container.end());
}

} // namespace mathter



// Generalized cross-product unfortunately needs matrix determinant.
#include "../Matrix/Math.hpp"

namespace mathter {
namespace impl {

	template <class IterFirst, class IterLast, class Vec = typename std::iterator_traits<IterFirst>::value_type>
	auto CrossND(IterFirst first, IterLast last) -> std::enable_if_t<is_vector_v<Vec>, Vec> {
		using Scalar = scalar_type_t<Vec>;
		constexpr auto Dim = dimension_v<Vec>;
		Vec result;
		Matrix<Scalar, Dim - 1, Dim - 1> detCalc;

		std::array<std::optional<std::reference_wrapper<const Vec>>, Dim - 1> vectors;
		auto [argIt, outIt] = std::tuple(first, vectors.begin());
		for (; argIt != last && outIt != vectors.end(); ++argIt, ++outIt) {
			*outIt = std::ref(*argIt);
		}
		if (outIt != vectors.end()) {
			throw std::invalid_argument("not enough arguments for cross product");
		}

		// Calculate elements of result on-by-one
		int sign = 2 * (Dim % 2) - 1;
		for (size_t idx = 0; idx < result.Dimension(); ++idx, sign *= -1) {
			// Fill up sub-matrix the determinant of which yields the coefficient of base-vector.
			for (int j = 0; j < idx; ++j) {
				for (int i = 0; i < detCalc.RowCount(); ++i) {
					detCalc(i, j) = (*vectors[i]).get()[j];
				}
			}
			for (int j = idx + 1; j < result.Dimension(); ++j) {
				for (int i = 0; i < detCalc.RowCount(); ++i) {
					detCalc(i, j - 1) = (*vectors[i]).get()[j];
				}
			}

			const Scalar coefficient = static_cast<Scalar>(sign) * Determinant(detCalc);
			result(idx) = coefficient;
		}

		return result;
	}

	template <class IterFirst, class IterLast, class Vec>
	auto Cross(IterFirst first, IterLast last) -> std::enable_if_t<is_vector_v<Vec>, Vec> {
		constexpr auto Dim = dimension_v<Vec>;

		if constexpr (dimension_v<Vec> == 2) {
			if (first == last) {
				throw std::invalid_argument("not enough arguments for cross product");
			}
			return Vec(-first->y, first->x);
		}
		if constexpr (dimension_v<Vec> == 3) {
			if (first == last) {
				throw std::invalid_argument("not enough arguments for cross product");
			}
			const auto& a = *first++;
			if (first == last) {
				throw std::invalid_argument("not enough arguments for cross product");
			}
			const auto& b = *first++;
			return Vec(a.yzx * b.zxy - a.zxy * b.yzx);
		}
		else {
			return CrossND(first, last);
		}
	}
} // namespace impl
} // namespace mathter


#if defined(MATHTER_MINMAX)
#pragma pop_macro("min")
#pragma pop_macro("max")
#endif