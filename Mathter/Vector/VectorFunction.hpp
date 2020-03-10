#pragma once

#if _MSC_VER && defined(min)
#pragma push_macro("min")
#pragma push_macro("max")
#undef min
#undef max
#define MATHTER_MINMAX
#endif


#include "VectorImpl.hpp"

namespace mathter {


// Misc
template <class T, int Dim, bool Packed>
inline void Fill(Vector<T, Dim, Packed>& lhs, T all) {
	if constexpr (!traits::HasSimd<Vector<T, Dim, Packed>>::value) {
		for (auto& v : lhs) {
			v = all;
		}
	}
	else {
		using SimdT = decltype(VectorData<T, Dim, Packed>::simd);
		lhs.simd = SimdT::spread(all);
	}
}

template <class T, int Dim, bool Packed>
inline T Dot(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	if constexpr (!traits::HasSimd<Vector<T, Dim, Packed>>::value) {
		T sum = T(0);
		for (int i = 0; i < Dim; ++i) {
			sum += lhs.data[i] * rhs.data[i];
		}
		return sum;
	}
	else {
		using SimdT = decltype(VectorData<T, Dim, Packed>::simd);
		return SimdT::template dot<Dim>(lhs.simd, rhs.simd);
	}
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
	return Vector<T, 3, Packed>(lhs.y * rhs.z - lhs.z * rhs.y,
								lhs.z * rhs.x - lhs.x * rhs.z,
								lhs.x * rhs.y - lhs.y * rhs.x);
}
/// <summary> Returns the 3-dimensional cross-product. </summary>
template <class T, bool Packed>
Vector<T, 3, Packed> Cross(const std::array<const Vector<T, 3, Packed>*, 2>& args) {
	return Cross(*(args[0]), *(args[1]));
}


/// <summary> Returns the element-wise minimum of arguments </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Min(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	Vector<T, Dim, Packed> res;
	for (int i = 0; i < lhs.Dimension(); ++i) {
		res[i] = std::min(lhs[i], rhs[i]);
	}
	return res;
}
/// <summary> Returns the element-wise maximum of arguments </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Max(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	Vector<T, Dim, Packed> res;
	for (int i = 0; i < lhs.Dimension(); ++i) {
		res[i] = std::max(lhs[i], rhs[i]);
	}
	return res;
}


/// <summary> Returns the euclidean distance between to vectors. </summary>
template <class T, class U, int Dim, bool Packed1, bool Packed2>
auto Distance(const Vector<T, Dim, Packed1>& lhs, const Vector<U, Dim, Packed2>& rhs) {
	return (lhs - rhs).Length();
}

/// <summary> Returns the normalized version of <paramref name="arg">. </summary>
template <class T, int Dim, bool Packed>
auto Normalized(const Vector<T, Dim, Packed>& arg) {
	return arg.Normalized();
}

} // namespace mathter



// Generalized cross-product unfortunately needs matrix determinant.
#include "../Matrix.hpp"

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