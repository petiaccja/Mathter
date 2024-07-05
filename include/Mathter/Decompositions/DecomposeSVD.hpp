// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "DecomposeQR.hpp"

#include <Mathter/Matrix/MatrixCast.hpp>

#include <algorithm>


namespace mathter {


/// <summary> A utility class that can do common operations with the singular value decomposition,
///		i.e. solving equation systems. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class DecompositionSVD {
	template <int Rows_, int Columns_>
	using MatrixT = Matrix<T, Rows_, Columns_, Order, Layout, Packed>;

	static constexpr int Sdim = std::min(Rows, Columns);
	static constexpr int Udim = Rows;
	static constexpr int Vdim = Columns;

public:
	MatrixT<Udim, Sdim> U;
	MatrixT<Sdim, Sdim> S;
	MatrixT<Sdim, Vdim> V;
};


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct DecompositionSVD_V2 {
	static constexpr int PDim = std::min(Rows, Columns);

	Matrix<T, Rows, PDim, Order, Layout, Packed> U;
	Vector<T, PDim> S;
	Matrix<T, PDim, Columns, Order, Layout, Packed> V;
};


namespace impl {

	template <class T, class... Rest>
	T MaximumAbsolute(T first, T second, Rest... rest) {
		const auto partial = std::max(std::abs(first), std::abs(second));
		if constexpr (sizeof...(Rest) != 0) {
			return MaximumAbsolute(partial, rest...);
		}
		return partial;
	}

	/// <summary> The RQ decomposition of a 2x2 matrix. </summary>
	/// <remarks>
	///	Decomposes A such that A = RQ, where R is upper triangular and Q is a rotation matrix.
	///	The layout is R = {{r11, r12}, {0, r22}} and Q = {{c, -s}, {s, c}}.
	/// </remarks>
	template <class T>
	struct DecompositionRQ2x2 {
		T r11;
		T r12;
		T r22;
		T cq;
		T sq;
	};

	/// <summary> The SVD if a 2x2 matrix. </summary>
	/// <remarks>
	/// Decomposes A such that A = USV, where U is a rotation matrix, S is diagonal, and V is also a rotation matrix.
	///	The layout is U = {{cu, -su}, {su, cu}}, S = {{s11, 0}, {0, s22}}, and V = {{cv, -sv}, {sv, cv}}.
	/// </remarks>
	template <class T>
	struct DecompositionSVD2x2 {
		T cu;
		T su;
		T s11;
		T s22;
		T cv;
		T sv;
	};

	template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	DecompositionRQ2x2<T> DecomposeRQ2x2(const Matrix<T, 2, 2, Order, Layout, Packed>& A) {
		const auto a11 = A(0, 0);
		const auto a12 = A(0, 1);
		const auto a21 = A(1, 0);
		const auto a22 = A(1, 1);

		// If a21 is not precisely zero, the scaling will take care of it and the rest works fine.
		if (a21 == T(0)) {
			return { a11, a12, a22, T(1), T(0) };
		}

		// Rescale matrix elements to avoid underflow and overflow.
		const auto scaleNum = MaximumAbsolute(a11, a12, std::numeric_limits<T>::min()); // Avoid NaN.
		const auto scaleDiv = MaximumAbsolute(a21, a22); // a21 is never zero. See if statement above.
		const auto a11s = a11 / scaleNum;
		const auto a12s = a12 / scaleNum;
		const auto a21s = a21 / scaleDiv;
		const auto a22s = a22 / scaleDiv;

		// Compute the cosine and sine of Q and the values of R.
		const auto divisor = std::sqrt(a21s * a21s + a22s * a22s);

		const auto cq = a22s / divisor;
		const auto sq = a21s / divisor;
		const auto r11 = scaleNum * ((a11s * a22s - a12s * a21s) / divisor);
		const auto r12 = scaleNum * ((a11s * a21s + a12s * a22s) / divisor);
		const auto r22 = scaleDiv * divisor;

		return { r11, r12, r22, cq, sq };
	}


	template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	DecompositionSVD2x2<T> DecomposeSVD2x2(const Matrix<T, 2, 2, Order, Layout, Packed>& A) {
		// The algorithm is not trivial, see the description and derivation in the docs.

		// C++17 does not have the math constants yet.
		constexpr auto sqrt2 = T(1.4142135623730950488016887242096980785696718753769480731766797379);
		constexpr auto rsqrt2 = T(0.7071067811865475244008443621048490392848359376884740365883398689);

		// RQ preconditioning.
		const auto rq = DecomposeRQ2x2(A);

		// Get the elements of the upper triangular R matrix.
		const auto a11 = rq.r11;
		const auto a12 = rq.r12;
		const auto a22 = rq.r22;

		// Rescale matrix elements to avoid underflow and overflow.
		const auto scaler = MaximumAbsolute(a11, a12, a22);
		if (scaler == T(0)) {
			return { T(1), T(0), T(0), T(0), T(1), T(0) };
		}
		const auto a11s = a11 / scaler;
		const auto a12s = a12 / scaler;
		const auto a22s = a22 / scaler;

		// Compute the cosine and sine of the rotation matrix V in R = USV.
		const auto z = (a11s - a22s) * (a11s + a22s) - a12s * a12s;
		const auto g = a11s * a12s;
		const auto zabs = std::abs(z);
		const auto zsign = std::copysign(T(1), z);
		const auto p1 = std::hypot(zabs, T(2) * g);
		const auto pab = std::sqrt((zabs + p1) / p1);

		const auto cv = p1 != 0 ? std::clamp(-pab * rsqrt2, -T(1), T(1)) : T(1);
		const auto sv = p1 != 0 ? std::clamp(zsign * sqrt2 * g / (p1 * pab), -T(1), T(1)) : T(0);

		// Compute US = RV^T using the computed cv and sv.
		const auto us11 = cv * a11 - sv * a12;
		const auto us12 = sv * a11 + cv * a12;
		const auto us21 = -sv * a22;
		const auto us22 = cv * a22;

		// Compute the singular values, i.e. the diagonal of matrix S.
		const auto detUs = us11 * us22 - us12 * us21; // det(S) must match det(US), as det(U) == 1.
		const auto s11 = std::hypot(us11, us21);
		const auto s22abs = std::hypot(us12, us22);
		const auto s22 = std::copysign(s22abs, detUs);

		// Compute the cosine and sine of the rotation matrix U.
		const auto [cu, su] = s11 > s22abs ? std::tuple{ us11 / s11, us21 / s11 } : std::tuple{ us22 / s22, -us12 / s22 };

		// Multiply V and Q to undo the RQ decomposition.
		const auto cvq = cv * rq.cq - sv * rq.sq;
		const auto svq = sv * rq.cq + cv * rq.sq;

		return { cu, su, s11, s22, cvq, svq };
	}


	template <class T>
	std::tuple<T, T> DiagonalizeSymmetric2x2(const T& a11, const T& aoff, const T& a22) {
		constexpr auto sqrt2 = T(1.4142135623730950488016887242096980785696718753769480731766797379);
		constexpr auto rsqrt2 = T(0.7071067811865475244008443621048490392848359376884740365883398689);

		const auto z = a11 - a22;
		const auto p1 = std::hypot(z, T(2) * aoff);
		const auto pab = std::sqrt((std::abs(z) + p1) / p1);
		auto cv = p1 != T(0) ? -pab * rsqrt2 : T(1);
		auto sv = p1 != T(0) ? -std::copysign(sqrt2, z) * aoff / (p1 * pab) : T(0);

		// Pick the solution with the smallest rotation angle.
		// I'm not sure if this works.
		// if (sv < T(0)) {
		//	 std::swap(cv, sv);
		// }
		// if (cv > sv) {
		// 	 std::swap(cv, sv);
		//   cv = -cv;
		// }

		return { cv, sv };
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, bool Packed>
	std::tuple<T, T, T> TransposeMultiplyPartial(const Matrix<T, Rows, Columns, Order, eMatrixLayout::COLUMN_MAJOR, Packed>& m, int p, int q) {
		const auto ata11 = Dot(m.stripes[p], m.stripes[p]);
		const auto ataoff = Dot(m.stripes[p], m.stripes[q]);
		const auto ata22 = Dot(m.stripes[q], m.stripes[q]);
		return { ata11, ataoff, ata22 };
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, bool Packed>
	void GivensRotate(Matrix<T, Rows, Columns, Order, eMatrixLayout::COLUMN_MAJOR, Packed>& m, int p, int q, T cv, T sv) {
		auto stripeP = m.stripes[p] * cv + m.stripes[q] * sv;
		auto stripeQ = m.stripes[q] * cv - m.stripes[p] * sv;
		std::tie(m.stripes[p], m.stripes[q]) = std::tuple(std::move(stripeP), std::move(stripeQ));
	}


	// I'm not sure how to solve the separation of U and S from X_{inf} when there are zero singular values in S.
	// This means division by zero when calculating U's rows, plus U's rows may not be properly orthogonal
	// due to poor precision when a singular value is really small.
	// Until this is fixed, there is the two-sided variant.
	template <class T, int Dim, eMatrixOrder Order, bool Packed>
	auto DecomposeSVDJacobiOneSided(const Matrix<T, Dim, Dim, Order, eMatrixLayout::COLUMN_MAJOR, Packed>& A) {
		using MyDecomposition = DecompositionSVD_V2<T, Dim, Dim, Order, eMatrixLayout::COLUMN_MAJOR, Packed>;
		constexpr auto tolerance = T(4) * std::numeric_limits<T>::epsilon();

		const auto scale = Max(Abs(A));
		auto X = A / scale;
		std::decay_t<decltype(A)> V = Identity();

		T maxErrorPrev = std::numeric_limits<T>::max();
		T maxError = std::nextafter(maxErrorPrev, T(0));
		while (maxError < maxErrorPrev) {
			maxErrorPrev = maxError;
			maxError = T(0);
			for (int p = 0; p < Dim; ++p) {
				for (int q = p + 1; q < Dim; ++q) {
					const auto [ata11, ataoff, ata22] = TransposeMultiplyPartial(X, p, q);
					const auto pqError = std::abs(ataoff);
					if (pqError > tolerance) {
						maxError = std::max(maxError, pqError);
						const auto [cv, sv] = DiagonalizeSymmetric2x2(ata11, ataoff, ata22);
						GivensRotate(X, p, q, cv, sv);
						GivensRotate(V, p, q, cv, -sv);
					}
				}
			}
		}

		auto XT = Transpose(X);
		Vector<T, Dim, Packed> S;
		for (int i = 0; i < Dim; ++i) {
			auto& row = XT.stripes[i];
			const auto s = Length(row);
			S(i) = s;
			if (s != T(0)) {
				row /= s;
			}
		}
		for (int i = 0; i < Dim; ++i) {
			if (S(i) == T(0)) {
				const auto& row = XT.stripes[i];
				const auto diagonal = std::sqrt(T(1) - LengthSquared(row));
				XT(i, i) = diagonal;
			}
		}
		const auto U = Transpose(XT);

		return MyDecomposition{ U, S * scale, V };
	}



	template <class T, int Dim, eMatrixOrder Order, bool Packed>
	auto DecomposeSVDJacobiTwoSided(const Matrix<T, Dim, Dim, Order, eMatrixLayout::COLUMN_MAJOR, Packed>& A) {
		using MyDecomposition = DecompositionSVD_V2<T, Dim, Dim, Order, eMatrixLayout::COLUMN_MAJOR, Packed>;
		using Mat = std::decay_t<decltype(A)>;
		constexpr auto tolerance = T(4) * std::numeric_limits<T>::epsilon();

		const auto scale = Max(Abs(A));
		Mat U;
		Mat X = A / scale;
		Mat V;

		T maxErrorPrev = std::numeric_limits<T>::max();
		T maxError = std::nextafter(maxErrorPrev, T(0));
		while (maxError < maxErrorPrev) {
			maxErrorPrev = maxError;
			maxError = T(0);
			for (int p = 0; p < Dim; ++p) {
				for (int q = p + 1; q < Dim; ++q) {
					const auto [xpp, xpq, xqp, xqq] = std::tie(X(p, p), X(p, q), x(q, p), x(q, q));
					const auto error = std::max(std::abs(xpq), std::abs(xqp));
					if (error > tolerance) {
						maxError = std::max(maxError, error);
						const auto svd2x2 = DecomposeSVD2x2({xpp, xpq, xqp, xqq});
						GivensRotate(X, p, q, cv, sv);
						GivensRotate(V, p, q, cv, -sv);
					}
				}
			}
		}
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto DecomposeSVD(Matrix<T, Rows, Columns, Order, Layout, Packed> m, std::true_type) {
		using Mat22 = Matrix<T, 2, 2>;
		using Vec4 = Vector<T, 4>;

		Matrix<T, Rows, Columns, Order, eMatrixLayout::COLUMN_MAJOR, false> B;
		Matrix<T, Rows, Columns, Order, eMatrixLayout::COLUMN_MAJOR, false> U;
		Matrix<T, Columns, Columns, Order, eMatrixLayout::COLUMN_MAJOR, false> V;

		// Precondition with QR if needed
		if (Rows > Columns) {
			auto [Q, R] = DecomposeQR(m);
			B = R;
			U = Q.template Submatrix<Rows, Columns>(0, 0);
			V = Identity();
		}
		else {
			B = m;
			U = Identity();
			V = Identity();
		}

		T tolerance = T(1e-6);

		T N = NormSquared(B);
		T s;

		do {
			s = 0;

			for (int i = 0; i < m.ColumnCount(); ++i) {
				for (int j = i + 1; j < m.ColumnCount(); ++j) {
					s += B(i, j) * B(i, j) + B(j, i) * B(j, i);

					const Mat22 Bsub = {
						B(i, i), B(i, j),
						B(j, i), B(j, j)
					};
					const auto svd2x2 = DecomposeSVD2x2(Bsub);

					// Apply givens rotations given by 2x2 SVD to working matrices
					// B = R(c1,s1)*B*R(c2,-s2)
					const Vec4 givensCoeffs = { svd2x2.cu, -svd2x2.su, svd2x2.su, svd2x2.cu };
					for (int col = 0; col < B.ColumnCount(); ++col) {
						const Vec4 bElems = { B(i, col), B(j, col), B(i, col), B(j, col) };
						const auto product = bElems * givensCoeffs;
						B(i, col) = product(0) + product(1);
						B(j, col) = product(2) + product(3);
					}
					auto coli = B.stripes[i];
					B.stripes[i] = svd2x2.cv * coli + svd2x2.sv * B.stripes[j];
					B.stripes[j] = -svd2x2.sv * coli + svd2x2.cv * B.stripes[j];

					// U = U*R(c1,s1);
					coli = U.stripes[i];
					U.stripes[i] = svd2x2.cu * coli + -svd2x2.su * U.stripes[j];
					U.stripes[j] = svd2x2.su * coli + svd2x2.cu * U.stripes[j];

					// V = V*R(c2,s2);
					auto coliv = V.stripes[i];
					V.stripes[i] = svd2x2.cv * coliv + svd2x2.sv * V.stripes[j];
					V.stripes[j] = -svd2x2.sv * coliv + svd2x2.cv * V.stripes[j];
				}
			}
		} while (s > tolerance * N);

		Matrix<T, Columns, Columns, Order, Layout, Packed> Sout;
		Sout = Zero();
		for (int i = 0; i < B.ColumnCount(); ++i) {
			Sout(i, i) = B(i, i);
		}

		return DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>{ U, Sout, Transpose(V) };
	}

	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto DecomposeSVD(Matrix<T, Rows, Columns, Order, Layout, Packed> m, std::false_type) {
		auto [U, S, V] = DecomposeSVD(Transpose(m), std::true_type{});
		return DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>{ Transpose(V), S, Transpose(U) };
	}


} // namespace impl


/// <summary> Calculates the thin SVD of the matrix. </summary>
/// <remarks> For wide matrices, V is wide while U and S square.
///		For tall matrices, U is tall while S and V square. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecomposeSVD(Matrix<T, Rows, Columns, Order, Layout, Packed> m) {
	return impl::DecomposeSVD(m, std::integral_constant<bool, (Rows >= Columns)>());
}



} // namespace mathter