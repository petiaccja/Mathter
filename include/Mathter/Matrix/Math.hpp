// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Math.hpp"
#include "Arithmetic.hpp"
#include "Matrix.hpp"


namespace mathter {


/// <summary> Returns the maximum element of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Min(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	T minElement = Min(m.stripes[0]);
	for (size_t stripe = 1; stripe < m.stripes.size(); ++stripe) {
		minElement = std::min(minElement, Min(m.stripes[stripe]));
	}
	return minElement;
}


/// <summary> Returns the maximum element of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Max(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	T maxElement = Max(m.stripes[0]);
	for (size_t stripe = 1; stripe < m.stripes.size(); ++stripe) {
		maxElement = std::max(maxElement, Max(m.stripes[stripe]));
	}
	return maxElement;
}


/// <summary> Computes a divisor that scales the matrix so that its largest element is close to 1. </summary>
/// <remarks> This is similar to dividing by the infinity norm, but it's much faster and
///		less accurate for complex numbers. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto ScaleElements(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	static_assert(!std::is_integral_v<T>, "No need to scale integer matrices.");

	remove_complex_t<T> scale(0);
	for (size_t stripe = 0; stripe < m.stripes.size(); ++stripe) {
		scale = std::max(scale, ScaleElements(m.stripes[stripe]));
	}
	return scale;
}


/// <summary> Returns the sum of the elements of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Sum(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	T sum(0);
	for (size_t stripe = 0; stripe < m.stripes.size(); ++stripe) {
		sum += Sum(m.stripes[stripe]);
	}
	return sum;
}


/// <summary> Returns the maximum element of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Abs(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	using R = remove_complex_t<T>;
	Matrix<R, Rows, Columns, Order, Layout, Packed> r;
	for (size_t stripe = 0; stripe < m.stripes.size(); ++stripe) {
		r.stripes[stripe] = Abs(m.stripes[stripe]);
	}
	return r;
}


/// <summary> Returns the real part of a real or complex matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Real(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	using R = remove_complex_t<T>;
	Matrix<R, Rows, Columns, Order, Layout, Packed> r;
	for (size_t stripe = 0; stripe < m.stripes.size(); ++stripe) {
		r.stripes[stripe] = Real(m.stripes[stripe]);
	}
	return r;
}


/// <summary> Returns the imaginary part of a real or complex matrix. </summary>
/// <remarks> Returns a zero matrix for real matrices. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Imag(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	using R = remove_complex_t<T>;
	Matrix<R, Rows, Columns, Order, Layout, Packed> r;
	for (size_t stripe = 0; stripe < m.stripes.size(); ++stripe) {
		r.stripes[stripe] = Imag(m.stripes[stripe]);
	}
	return r;
}


/// <summary> Returns the maximum element of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Conj(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	Matrix<T, Rows, Columns, Order, Layout, Packed> r;
	for (size_t stripe = 0; stripe < m.stripes.size(); ++stripe) {
		r.stripes[stripe] = Conj(m.stripes[stripe]);
	}
	return r;
}


/// <summary> Calculates the square of the Frobenius norm of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto NormSquared(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	return Sum(Real(Hadamard(m, Conj(m))));
}


/// <summary> Calculates the Frobenius norm of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Norm(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	return std::sqrt(NormSquared(m));
}


/// <summary> Returns the Frobenius norm of the matrix, avoids overflow and underflow, so it's more expensive. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto NormPrecise(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	if constexpr (is_complex_v<T>) {
		const auto re = Abs(Real(m));
		const auto im = Abs(Imag(m));
		const auto maxElement = std::max(Max(re), Max(im));
		if (maxElement == static_cast<remove_complex_t<T>>(0)) {
			return static_cast<remove_complex_t<T>>(0);
		}
		const auto reScaled = re / maxElement;
		const auto imScaled = im / maxElement;
		const auto sq = Hadamard(reScaled, reScaled) + Hadamard(imScaled, imScaled);
		return std::sqrt(Sum(sq)) * maxElement;
	}
	else {
		const auto re = Abs(Real(m));
		const auto maxElement = Max(re);
		if (maxElement == static_cast<remove_complex_t<T>>(0)) {
			return static_cast<remove_complex_t<T>>(0);
		}
		const auto reScaled = re / maxElement;
		const auto sq = Hadamard(reScaled, reScaled);
		return std::sqrt(Sum(sq)) * maxElement;
	}
}


/// <summary> Returns the trace (sum of diagonal elements) of the matrix. </summary>
template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Trace(const Matrix<T, Dim, Dim, Order, Layout, Packed>& m) {
	Vector<T, Dim, Packed> diagonal;
	for (size_t i = 0; i < Dim; ++i) {
		diagonal[i] = m(i, i);
	}
	return Sum(diagonal);
}


/// <summary> Returns the determinant of a 2x2 matrix. </summary>
template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Determinant(const Matrix<T, 2, 2, Order, Layout, Packed>& m) {
	return m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
}

/// <summary> Returns the determinant of a 3x3 matrix. </summary>
template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Determinant(const Matrix<T, 3, 3, Order, Layout, Packed>& m) {
	const Vector r0_zyx = m.stripes[0].zyx;
	const Vector r1_xzy = m.stripes[1].xzy;
	const Vector r1_yxz = m.stripes[1].yxz;
	const Vector r2_yxz = m.stripes[2].yxz;
	const Vector r2_xzy = m.stripes[2].xzy;

	T det = Sum(r0_zyx * (r1_xzy * r2_yxz - r1_yxz * r2_xzy));

	return det;
}

/// <summary> Returns the determinant of a 4x4 matrix. </summary>
template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Determinant(const Matrix<T, 4, 4, Order, Layout, Packed>& m) {
	const Vector<T, 4, Packed> evenPair = { 1, -1, -1, 1 };
	const Vector<T, 4, Packed> oddPair = { -1, 1, 1, -1 };

	const auto& r0 = m.stripes[0];
	const auto& r1 = m.stripes[1];
	const auto& r2 = m.stripes[2];
	const auto& r3 = m.stripes[3];

	const Vector r2_zwzw = r2.zwzw;
	const Vector r0_yyxx = r0.yyxx;
	const Vector r1_wwxy = r1.wwxy;
	const Vector r2_xyzz = r2.xyzz;
	const Vector r3_wwww = r3.wwww;
	const Vector r1_zzxy = r1.zzxy;
	const Vector r0_yxyx = r0.yxyx;
	const Vector r3_xxyy = r3.xxyy;
	const Vector r1_wzwz = r1.wzwz;
	const Vector r2_xyww = r2.xyww;
	const Vector r3_zzzz = r3.zzzz;

	const Vector r2_yxz = r2.yxz;
	const Vector r3_xzy = r3.xzy;
	const Vector r2_xzy = r2.xzy;
	const Vector r3_yxz = r3.yxz;
	const Vector r2_yxw = r2.yxw;
	const Vector r1_zyx = r1.zyx;
	const Vector r3_yxw = r3.yxw;
	const Vector r2_xwy = r2.xwy;
	const Vector r3_xwy = r3.xwy;
	const Vector r1_wyx = r1.wyx;
	const T r0_w = r0.w;
	const T r0_z = r0.z;

	const T det = Sum(evenPair * (r0_yyxx * r1_wzwz * r2_zwzw * r3_xxyy))
				  + Sum(oddPair * (r0_yxyx * r1_wwxy * r2_xyww * r3_zzzz))
				  + Sum(evenPair * (r0_yxyx * r1_zzxy * r2_xyzz * r3_wwww))
				  + (r0_w * Sum(r1_zyx * (r2_yxz * r3_xzy - r2_xzy * r3_yxz)))
				  + (r0_z * Sum(r1_wyx * (r2_xwy * r3_yxw - r2_yxw * r3_xwy)));

	return det;
}


namespace impl {

	template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	T BareissAlgorithm(const Matrix<T, Dim, Dim, Order, Layout, Packed>& m) {
		auto inPlace = m;
		bool flipSign = false;
		for (size_t k = 0; k < Dim - 1; ++k) {
			auto minor = k == 0 ? static_cast<T>(1) : inPlace(k - 1, k - 1);
			for (size_t i = k; i < Dim && minor == static_cast<T>(0); ++i) {
				assert(k > 0); // Should be ensured as minor is always 1 for k=0.
				const auto& candidate = inPlace(i, k - 1);
				if (candidate != static_cast<T>(0)) {
					auto rowk1 = inPlace.Row(k - 1);
					auto rowi = inPlace.Row(i);
					inPlace.Row(k - 1, rowi);
					inPlace.Row(i, rowk1);
					minor = candidate;
					flipSign = !flipSign;
				}
			}

			// I'm not 100% sure this is correct, but what else?
			if (minor == static_cast<T>(0)) {
				return static_cast<T>(0);
			}

			for (size_t i = k + 1; i < Dim; ++i) {
				for (size_t j = k + 1; j < Dim; ++j) {
					inPlace(i, j) = (inPlace(i, j) * inPlace(k, k) - inPlace(i, k) * inPlace(k, j)) / minor;
				}
			}
		}

		return flipSign ? -inPlace(Dim - 1, Dim - 1) : inPlace(Dim - 1, Dim - 1);
	}

} // namespace impl


/// <summary> Returns the determinant of the matrix. </summary>
template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T Determinant(const Matrix<T, Dim, Dim, Order, Layout, Packed>& m) {
	return impl::BareissAlgorithm(m);
}

/// <summary> Transposes the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Transpose(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	Matrix<T, Columns, Rows, Order, Layout, Packed> result;
	for (int i = 0; i < m.RowCount(); ++i) {
		for (int j = 0; j < m.ColumnCount(); ++j) {
			result(j, i) = m(i, j);
		}
	}
	return result;
}


/// <summary> Return the conjugate transpose of the matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto ConjTranspose(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	return Conj(Transpose(m));
}


/// <summary> Returns the inverse of a 2x2 matrix. </summary>
template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Inverse(const Matrix<T, 2, 2, Order, Layout, Packed>& m) {
	static_assert(!std::is_integral_v<T>, "Integer matrices cannot be inverted.");

	Matrix<T, 2, 2, Order, Layout, Packed> result;

	const auto& r0 = m.stripes[0];
	const auto& r1 = m.stripes[1];

	result.stripes[0] = { r1.y, -r0.y };
	result.stripes[1] = { -r1.x, r0.x };

	const T det = r0.x * r1.y - r0.y * r1.x;
	return result / det;
}


/// <summary> Returns the inverse of a 3x3 matrix. </summary>
template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Inverse(const Matrix<T, 3, 3, Order, Layout, Packed>& m) {
	static_assert(!std::is_integral_v<T>, "Integer matrices cannot be inverted.");

	// This code below uses notation for row-major matrices' stripes.
	// It, however, "magically" works for column-major layout as well.
	const Vector r0_zxy = m.stripes[0].zxy;
	const Vector r0_yzx = m.stripes[0].yzx;
	const Vector r1_yzx = m.stripes[1].yzx;
	const Vector r1_zxy = m.stripes[1].zxy;
	const Vector r2_zxy = m.stripes[2].zxy;
	const Vector r2_yzx = m.stripes[2].yzx;

	const Vector c0 = r1_yzx * r2_zxy - r1_zxy * r2_yzx;
	const Vector c1 = r0_zxy * r2_yzx - r0_yzx * r2_zxy;
	const Vector c2 = r0_yzx * r1_zxy - r0_zxy * r1_yzx;

	const Vector r0_zyx = m.stripes[0].zyx;
	const Vector r1_xzy = m.stripes[1].xzy;
	const Vector r1_yxz = m.stripes[1].yxz;
	const Vector r2_yxz = m.stripes[2].yxz;
	const Vector r2_xzy = m.stripes[2].xzy;

	const Matrix<T, 3, 3, Order, Layout, Packed> result{
		stripeArg,
		Vector{ c0[0], c1[0], c2[0] },
		Vector{ c0[1], c1[1], c2[1] },
		Vector{ c0[2], c1[2], c2[2] },
	};

	const T det = Sum(r0_zyx * (r1_xzy * r2_yxz - r1_yxz * r2_xzy));
	return result / det;
}


/// <summary> Returns the inverse of a 4x4 matrix. </summary>
template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Inverse(const Matrix<T, 4, 4, Order, Layout, Packed>& m) {
	static_assert(!std::is_integral_v<T>, "Integer matrices cannot be inverted.");

	const Vector<T, 4, false> even = { 1, -1, 1, -1 };
	const Vector<T, 4, false> odd = { -1, 1, -1, 1 };
	const Vector<T, 4, false> evenPair = { 1, -1, -1, 1 };
	const Vector<T, 4, false> oddPair = { -1, 1, 1, -1 };

	const auto& r0 = m.stripes[0];
	const auto& r1 = m.stripes[1];
	const auto& r2 = m.stripes[2];
	const auto& r3 = m.stripes[3];

	const Vector r0_wwwz = r0.wwwz;
	const Vector r0_yxxx = r0.yxxx;
	const Vector r0_zzyy = r0.zzyy;
	const Vector r1_wwwz = r1.wwwz;
	const Vector r1_yxxx = r1.yxxx;
	const Vector r1_zzyy = r1.zzyy;
	const Vector r2_wwwz = r2.wwwz;
	const Vector r2_yxxx = r2.yxxx;
	const Vector r2_zzyy = r2.zzyy;
	const Vector r3_wwwz = r3.wwwz;
	const Vector r3_yxxx = r3.yxxx;
	const Vector r3_zzyy = r3.zzyy;

	const Vector r0_wwwz_r1_yxxx = r0_wwwz * r1_yxxx;
	const Vector r0_wwwz_r1_zzyy = r0_wwwz * r1_zzyy;
	const Vector r0_yxxx_r1_wwwz = r0_yxxx * r1_wwwz;
	const Vector r0_yxxx_r1_zzyy = r0_yxxx * r1_zzyy;
	const Vector r0_zzyy_r1_wwwz = r0_zzyy * r1_wwwz;
	const Vector r0_zzyy_r1_yxxx = r0_zzyy * r1_yxxx;
	const Vector r2_wwwz_r3_yxxx = r2_wwwz * r3_yxxx;
	const Vector r2_wwwz_r3_zzyy = r2_wwwz * r3_zzyy;
	const Vector r2_yxxx_r3_wwwz = r2_yxxx * r3_wwwz;
	const Vector r2_yxxx_r3_zzyy = r2_yxxx * r3_zzyy;
	const Vector r2_zzyy_r3_wwwz = r2_zzyy * r3_wwwz;
	const Vector r2_zzyy_r3_yxxx = r2_zzyy * r3_yxxx;

	const Vector c0 = odd * (r1_wwwz * r2_zzyy_r3_yxxx - r1_zzyy * r2_wwwz_r3_yxxx - r1_wwwz * r2_yxxx_r3_zzyy + r1_yxxx * r2_wwwz_r3_zzyy + r1_zzyy * r2_yxxx_r3_wwwz - r1_yxxx * r2_zzyy_r3_wwwz);
	const Vector c1 = even * (r0_wwwz * r2_zzyy_r3_yxxx - r0_zzyy * r2_wwwz_r3_yxxx - r0_wwwz * r2_yxxx_r3_zzyy + r0_yxxx * r2_wwwz_r3_zzyy + r0_zzyy * r2_yxxx_r3_wwwz - r0_yxxx * r2_zzyy_r3_wwwz);
	const Vector c2 = odd * (r0_wwwz_r1_zzyy * r3_yxxx - r0_zzyy_r1_wwwz * r3_yxxx - r0_wwwz_r1_yxxx * r3_zzyy + r0_yxxx_r1_wwwz * r3_zzyy + r0_zzyy_r1_yxxx * r3_wwwz - r0_yxxx_r1_zzyy * r3_wwwz);
	const Vector c3 = even * (r0_wwwz_r1_zzyy * r2_yxxx - r0_zzyy_r1_wwwz * r2_yxxx - r0_wwwz_r1_yxxx * r2_zzyy + r0_yxxx_r1_wwwz * r2_zzyy + r0_zzyy_r1_yxxx * r2_wwwz - r0_yxxx_r1_zzyy * r2_wwwz);

	const Matrix<T, 4, 4, Order, Layout, Packed> result{
		stripeArg,
		Vector{ c0[0], c1[0], c2[0], c3[0] },
		Vector{ c0[1], c1[1], c2[1], c3[1] },
		Vector{ c0[2], c1[2], c2[2], c3[2] },
		Vector{ c0[3], c1[3], c2[3], c3[3] },
	};

	const Vector r2_zwzw = r2.zwzw;
	const Vector r0_yyxx = r0.yyxx;
	const Vector r1_wwxy = r1.wwxy;
	const Vector r2_xyzz = r2.xyzz;
	const Vector r3_wwww = r3.wwww;
	const Vector r1_zzxy = r1.zzxy;
	const Vector r0_yxyx = r0.yxyx;
	const Vector r3_xxyy = r3.xxyy;
	const Vector r1_wzwz = r1.wzwz;
	const Vector r2_xyww = r2.xyww;
	const Vector r3_zzzz = r3.zzzz;

	const Vector r2_yxz = r2.yxz;
	const Vector r3_xzy = r3.xzy;
	const Vector r2_xzy = r2.xzy;
	const Vector r3_yxz = r3.yxz;
	const Vector r2_yxw = r2.yxw;
	const Vector r1_zyx = r1.zyx;
	const Vector r3_yxw = r3.yxw;
	const Vector r2_xwy = r2.xwy;
	const Vector r3_xwy = r3.xwy;
	const Vector r1_wyx = r1.wyx;
	const T r0_w = r0.w;
	const T r0_z = r0.z;

	const T det = Sum(evenPair * (r0_yyxx * r1_wzwz * r2_zwzw * r3_xxyy))
				  + Sum(oddPair * (r0_yxyx * r1_wwxy * r2_xyww * r3_zzzz))
				  + Sum(evenPair * (r0_yxyx * r1_zzxy * r2_xyzz * r3_wwww))
				  + (r0_w * Sum(r1_zyx * (r2_yxz * r3_xzy - r2_xzy * r3_yxz)))
				  + (r0_z * Sum(r1_wyx * (r2_xwy * r3_yxw - r2_yxw * r3_xwy)));

	return result / det;
}

} // namespace mathter

#include "../Decompositions/DecomposeLU.hpp"

namespace mathter {

/// <summary> Returns the inverse of the matrix. </summary>
template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
Matrix<T, Dim, Dim, Order, Layout, Packed> Inverse(const Matrix<T, Dim, Dim, Order, Layout, Packed>& m) {
	static_assert(!std::is_integral_v<T>, "Integer matrices cannot be inverted.");

	return DecomposeLUP(m).Inverse();
}

} // namespace mathter