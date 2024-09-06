// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "DecomposeQR.hpp"

#include <Mathter/Matrix/Matrix.hpp>
#include <Mathter/Transforms/IdentityBuilder.hpp>

#include <algorithm>
#include <array>


namespace mathter {


/// <summary> The singular value decomposition of a matrix. </summary>
/// <remarks> This is a thin/compact SVD. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct DecompositionSVD {
	using Real = remove_complex_t<T>;
	static constexpr auto PDim = std::min(Rows, Columns);

	Matrix<T, Rows, PDim, Order, Layout, Packed> U;
	Vector<Real, PDim, Packed> S;
	Matrix<T, PDim, Columns, Order, Layout, Packed> V;

	/// <summary> Solve multiple linear systems of equations at the same time. </summary>
	/// <remarks> For overdetermined systems, it returns the least squares solution. </remarks>
	template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
	auto Solve(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const;

	/// <summary> Solve a linear systems of equations. </summary>
	/// <remarks> For overdetermined systems, it returns the least squares solution. </remarks>
	template <class T2, bool Packed2>
	auto Solve(const Vector<T2, std::max(Rows, Columns), Packed2>& b) const;

	/// <summary> Compute the inverse or the pseudoinverse of the matrix. </summary>
	/// <remarks> For non-square matrices, the pseudoinverse is computed instead. </remarks>
	auto Inverse() const -> Matrix<T, Columns, Rows, Order, Layout, Packed>;
};


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, class = std::enable_if_t<(Rows >= Columns)>>
DecompositionSVD(const Matrix<T, Rows, Columns, Order, Layout, Packed>&,
				 const Vector<remove_complex_t<T>, Columns, Packed>&,
				 const Matrix<T, Columns, Columns, Order, Layout, Packed>&) -> DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>;


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, class = std::enable_if_t<(Rows < Columns)>>
DecompositionSVD(const Matrix<T, Rows, Rows, Order, Layout, Packed>&,
				 const Vector<remove_complex_t<T>, Rows, Packed>&,
				 const Matrix<T, Rows, Columns, Order, Layout, Packed>&) -> DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>;


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>::Solve(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const {
	if constexpr (Order == eMatrixOrder::PRECEDE_VECTOR) {
		static_assert(Rows == Rows2);
		auto X = ConjTranspose(U) * b;
		for (size_t row = 0; row < Columns; ++row) {
			X.Row(row, X.Row(row) / S(row));
		}
		return ConjTranspose(V) * X;
	}
	else {
		static_assert(Columns == Columns2);
		auto X = b * ConjTranspose(V);
		for (size_t col = 0; col < Rows; ++col) {
			X.Column(col, X.Column(col) / S(col));
		}
		return X * ConjTranspose(U);
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, bool Packed2>
auto DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>::Solve(const Vector<T2, std::max(Rows, Columns), Packed2>& b) const {
	if constexpr (Order == eMatrixOrder::PRECEDE_VECTOR) {
		return Vector(Solve(Matrix<T2, Rows, 1, Order, Layout, Packed>(b)));
	}
	else {
		return Vector(Solve(Matrix<T2, 1, Columns, Order, Layout, Packed>(b)));
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>::Inverse() const -> Matrix<T, Columns, Rows, Order, Layout, Packed> {
	auto Xi = ConjTranspose(V);
	for (size_t i = 0; i < S.Dimension(); ++i) {
		Xi.Column(i, Xi.Column(i) / S(i));
	}
	const auto Ui = ConjTranspose(U);
	return Xi * Ui;
}


namespace impl {

	/// <summary> Find a suitable divisor that scales the input numbers to be roughly 1.0 in magnitude. </summary>
	template <class T1, class... Ts>
	auto ScaleElements(T1 first, Ts... rest) {
		if constexpr (is_complex_v<T1> || (... || is_complex_v<Ts>)) {
			return std::max(ScaleElements(std::real(first), std::real(rest)...),
							ScaleElements(std::imag(first), std::imag(rest)...));
		}
		else {
			if constexpr (sizeof...(rest) == 0) {
				return std::abs(first);
			}
			else {
				return std::max(std::abs(first), ScaleElements(rest...));
			}
		}
	}

	/// <summary> The RQ decomposition of a 2x2 matrix. </summary>
	/// <remarks>
	///	Decomposes A such that A = RQ, where R is upper triangular and Q is a rotation matrix.
	///	The layout is R = {{r11, r12}, {0, r22}} and Q = {{c, -s*}, {s, c*}}.
	/// </remarks>
	template <class T>
	struct DecompositionRQ2x2 {
		T r11;
		T r12;
		T r22;
		T cq;
		T sq;
	};

	/// <summary> The SVD of a 2x2 matrix. </summary>
	/// <remarks>
	/// Decomposes A such that A = USV, where U is a rotation matrix, S is diagonal, and V is also a rotation matrix.
	///	The layout is U = {{cu, -su*}, {su, cu*}}, S = {{s11, 0}, {0, s22}}, and V = {{cv, -sv*}, {sv, cv*}}.
	///	The determinant must be used to multiply either the second column of U or s22. For the real case,
	///	det is essentially a plus or minus, for the complex case, it's an arbitrary complex number of unit length.
	/// </remarks>
	template <class T>
	struct DecompositionSVD2x2 {
		T cu;
		T su;
		remove_complex_t<T> s11;
		remove_complex_t<T> s22;
		T cv;
		T sv;
		T det; // Multiply either s22 or the 2nd column of U by det, depending on whether you want proper rotations or positive singular values.
	};


	/// <summary> Modifies a diagonalizing rotation such that its angle is minimal. </summary>
	/// <remarks> There are multiple rotation matrices V that diagonalize a Hermitian matrix V.
	///		Given one such rotation matrix, this method returns the one that has the smallest
	///		angle of rotation. Minimal angle rotations are required to ensure the convergence
	///		of the both 1- and 2-sided Jacobi SVD algorithms. </remarks>
	template <class T>
	std::tuple<T, T> MinimizeDiagonalizingRotation(T cv, T sv) {
		// Rotate V by Q = [0, -1; 1, 0] as in V' = V*Q.
		// This merely flips the eigenvalues in V^T * A^T * A * V, but leaves the off-diagonals zero.
		auto cvMag = std::abs(cv);
		auto svMag = std::abs(sv);
		if (cvMag < svMag) {
			std::tie(cv, sv) = std::tuple(-conj()(sv), conj()(cv));
			std::swap(cvMag, svMag);
		}

		// Rotate V by Q = [q, 0; 0, conj(q)] as in V' = V*Q.
		// For reals, this flips the sign of cv to be positive.
		// For complexes, this rotates cv in the complex plane such that it's purely a positive real.
		const auto q = conj()(cv / cvMag);
		std::tie(cv, sv) = std::tuple(q * cv, q * sv); // Make cv positive real.

		return { cv, sv };
	}


	/// <summary> Diagonalizes a Hermitian (real or complex) matrix. </summary>
	/// <returns> The cosine & sine (i.e. left column) of the diagonalizing rotation matrix. </returns>
	/// <remarks> The returned unitary matrix V diagonalizes the Hermitian
	///		input matrix H such that V^T * H * V = D is diagonal. In the complex case
	///		^T means conjugate transpose. </remarks>
	template <class T>
	std::tuple<T, T> DiagonalizeHermitian2x2(const remove_complex_t<T>& a11, const T& aoff, const remove_complex_t<T>& a22) {
		using Real = remove_complex_t<T>;

		auto [cv, sv] = std::tuple{ T(1), T(0) };
		if (aoff == T(0)) {
			return { cv, sv };
		}

		const auto z = a22 - a11;
		const auto d = std::hypot(z, Real(2) * std::real(aoff), Real(2) * std::imag(aoff));
		const auto Lma11 = Real(0.5) * (z + std::copysign(d, z));

		const auto aoffMagApprox = ScaleElements(aoff);
		const auto Lma11Mag = std::abs(Lma11);

		std::tie(cv, sv) = aoffMagApprox > Lma11Mag ?
							   std::tuple{ T(1), Lma11 / aoff } :
							   std::tuple{ aoff / Lma11, T(1) };


		const auto scale = std::sqrt(std::norm(cv) + std::norm(sv));

		std::tie(cv, sv) = std::tuple(cv / scale, sv / scale);

		return { cv, sv };
	}


	/// <summary> Diagonalizes an upper triangular matrix. </summary>
	/// <returns> The cosine & sine (i.e. left column) of the diagonalizing rotation matrix. </returns>
	/// <remarks> The returned unitary matrix V diagonalizes the upper triangular
	///		input matrix R such that V^T * R^T * R * V = D is diagonal. In the complex case
	///		^T means conjugate transpose. </remarks>
	template <class T>
	std::tuple<T, T> DiagonalizeTriangular2x2(const T& r11, const T& r12, const T& r22) {
		using Real = remove_complex_t<T>;

		const auto zr = (std::real(r22) + std::real(r11)) * (std::real(r22) - std::real(r11));
		const auto zi = (std::imag(r22) + std::imag(r11)) * (std::imag(r22) - std::imag(r11));
		const auto [zh, zl] = Fast2Sum(zr, zi);
		const auto z = zh + std::norm(r12) + zl;
		const auto d = std::sqrt(Real(4) * std::norm(r11) * std::norm(r12) + z * z);
		const auto u = (z + std::copysign(d, z)) / (Real(2) * r12 * conj()(r11));
		if (std::isfinite(std::real(u)) && std::isfinite(std::imag(u))) {
			const auto scale = std::hypot(Real(1), std::real(u), std::imag(u));
			const auto [cv1, sv1] = std::tuple{ T(Real(1) / scale), u / scale };
			return { cv1, sv1 };
		}
		else {
			return { T(1), T(0) };
		}
	}


	/// <summary> Computes the RQ of a 2x2 matrix. </summary>
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
		const auto scaleA21 = ScaleElements(a21);
		const auto scaleA22 = ScaleElements(a22);
		const auto scaleDiv = ScaleElements(scaleA21, scaleA22); // a21 is never zero. See if statement above.
		const auto a21s = a21 / scaleDiv;
		const auto a22s = a22 / scaleDiv;

		// Compute the cosine and sine of Q and the values of R.
		const auto divisor = std::sqrt(std::norm(a21s) + std::norm(a22s));

		const auto cq = conj()(a22s / divisor);
		const auto sq = a21s / divisor;

		const auto r11 = a11 * conj()(cq) - a12 * sq;
		const auto r12 = a11 * conj()(sq) + a12 * cq;
		const auto r22 = scaleDiv * (scaleA21 > scaleA22 ? a21s / sq : a22s / conj{}(cq));

		return { r11, r12, r22, cq, sq };
	}


	/// <summary> Computes the singular value decomposition of a 2x2 matrix. </summary>
	///	<remarks> This uses an explicit algorithm and is meant to be used as the kernel
	///		of the 2-sided Jacobi SVD algorithm. </remarks>
	template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	DecompositionSVD2x2<T> DecomposeSVD2x2(const Matrix<T, 2, 2, Order, Layout, Packed>& A) {
		using Real = remove_complex_t<T>;

		// Scale A such that its elements are near 1.0.
		// This avoids overflow for the diagonalization and computation of U as well.
		// The RQ decompositions does its own scaling slightly differently, may wanna get rid of that.
		const auto scale = mathter::ScaleElements(A);
		if (scale == Real(0)) {
			return { T(1), T(0), Real(0), Real(0), T(1), T(0), T(1) };
		}
		const auto As = A / scale;

		// RQ preconditioning. This should help preserve floating point accuracy
		// by removing terms from the diagonalization step (see below).
		const auto rq = DecomposeRQ2x2(As);

		// Compute V1: A = U * S * V1^T * Q
		const auto [cv1, sv1] = DiagonalizeTriangular2x2(rq.r11, rq.r12, rq.r22);

		// Compute V = Q^T * V1 to "undo" the RQ decomposition.
		auto cv = conj()(rq.cq) * cv1 + conj()(rq.sq) * sv1;
		auto sv = -rq.sq * cv1 + rq.cq * sv1;
		std::tie(cv, sv) = MinimizeDiagonalizingRotation(cv, sv);

		// Compute U * S = A * V.
		const Matrix<T, 2, 2, Order, Layout, Packed> V = { cv, -conj()(sv), sv, conj()(cv) };
		const auto US = As * V;

		// Compute the singular values, i.e. the diagonal of matrix S.
		const auto s11 = Length(US.Column(0));
		const auto s22 = Length(US.Column(1));

		// Compute the cosine and sine of the rotation matrix U.
		const auto det = Determinant(US);
		const auto detScaled = det / ScaleElements(det, std::numeric_limits<Real>::min());
		const auto sign = detScaled != T(0) ? detScaled / std::abs(detScaled) : T(1);
		const auto chkSign = std::abs(sign);
		const auto s22signed = s22 / sign;
		const auto [cu, su] = s11 > s22 ?
								  std::tuple{ US(0, 0) / s11, US(1, 0) / s11 } :
								  std::tuple{ conj()(US(1, 1)) / s22signed, -conj()(US(0, 1)) / s22signed };

		return { cu, su, scale * s11, scale * s22, conj()(cv), -sv, sign };
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	void GivensRotateRight(Matrix<T, Rows, Columns, Order, Layout, Packed>& m, int p, int q, T cv, T sv, T det = T(1)) {
		const auto colP = m.Column(p) * cv + m.Column(q) * sv;
		const auto colQ = det * (m.Column(q) * conj{}(cv)-m.Column(p) * conj{}(sv));
		m.Column(p, colP);
		m.Column(q, colQ);
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	void GivensRotateLeft(Matrix<T, Rows, Columns, Order, Layout, Packed>& m, int p, int q, T cv, T sv, T det = T(1)) {
		const auto rowP = m.Row(p) * cv - det * m.Row(q) * conj{}(sv);
		const auto rowQ = det * m.Row(q) * conj{}(cv) + m.Row(p) * sv;
		m.Row(p, rowP);
		m.Row(q, rowQ);
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto DecomposeSVDJacobiTwoSided(const Matrix<T, Rows, Columns, Order, Layout, Packed>& A) {
		static_assert(Rows >= Columns);

		using Real = remove_complex_t<T>;
		using MatExt = Matrix<T, Rows, Rows, Order, Layout, Packed>;
		constexpr auto tolerance = Real(1e-1f) * std::numeric_limits<Real>::epsilon();

		const auto scaler = ScaleElements(A);
		MatExt U = Identity();
		MatExt X;
		MatExt V = Identity();
		X.Insert(0, 0, A / scaler);
		for (size_t col = Columns; col < Rows; ++col) {
			for (size_t row = 0; row < Rows; ++row) {
				X(row, col) = T(0);
			}
		}

		auto maxErrorPrev = std::numeric_limits<Real>::max();
		auto maxError = std::nextafter(maxErrorPrev, Real(0));
		while (maxError < maxErrorPrev) {
			maxErrorPrev = maxError;
			maxError = Real(0);
			for (int p = 0; p < Rows; ++p) {
				for (int q = p + 1; q < Rows; ++q) {
					const auto [xpp, xpq, xqp, xqq] = std::tuple(X(p, p), X(p, q), X(q, p), X(q, q));
					const auto error = ScaleElements(xpq, xqp);
					if (error > tolerance) {
						maxError = std::max(maxError, error);
						const auto svd2x2 = DecomposeSVD2x2(Matrix<T, 2, 2>{ xpp, xpq, xqp, xqq });

						GivensRotateRight(U, p, q, svd2x2.cu, svd2x2.su, svd2x2.det);
						GivensRotateLeft(X, p, q, conj()(svd2x2.cu), -conj()(svd2x2.det) * svd2x2.su, conj()(svd2x2.det));
						GivensRotateRight(X, p, q, conj()(svd2x2.cv), -svd2x2.sv);
						GivensRotateLeft(V, p, q, svd2x2.cv, svd2x2.sv);
					}
				}
			}
		}

		Vector<Real, Columns, Packed> S;
		for (size_t i = 0; i < Columns; ++i) {
			S(i) = std::real(X(i, i));
		}

		return DecompositionSVD{ U.template Extract<Rows, Columns>(0, 0), S * scaler, V.template Extract<Columns, Columns>(0, 0) };
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	std::tuple<remove_complex_t<T>, T, remove_complex_t<T>> TransposeMultiplyPartial(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m, int p, int q) {
		const auto ata11 = Sum(Conj(m.Column(p)) * m.Column(p));
		const auto ataoff = Sum(Conj(m.Column(p)) * m.Column(q));
		const auto ata22 = Sum(Conj(m.Column(q)) * m.Column(q));
		return { std::real(ata11), ataoff, std::real(ata22) };
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto DecomposeSVDJacobiOneSided(const Matrix<T, Rows, Columns, Order, Layout, Packed>& A) {
		static_assert(Rows >= Columns);

		using Real = remove_complex_t<T>;

		using MatSquare = Matrix<T, Columns, Columns, Order, Layout, Packed>;
		using MatTall = Matrix<T, Rows, Columns, Order, Layout, Packed>;

		const auto scaler = ScaleElements(A);
		MatTall X = A / scaler;
		MatSquare V = Identity();

		auto maxErrorPrev = std::numeric_limits<Real>::max();
		auto maxError = std::nextafter(maxErrorPrev, Real(0));
		while (maxError < maxErrorPrev) {
			maxErrorPrev = maxError;
			maxError = Real(0);
			for (int p = 0; p < Columns; ++p) {
				for (int q = p + 1; q < Columns; ++q) {
					const auto [ata11, ataoff, ata22] = TransposeMultiplyPartial(X, p, q);
					const auto chk = ConjTranspose(X) * X;
					const auto error = std::abs(ataoff);
					if (error != T(0)) {
						maxError = std::max(maxError, error);
						const auto [cv0, sv0] = DiagonalizeHermitian2x2(ata11, ataoff, ata22);
						const auto [cv, sv] = MinimizeDiagonalizingRotation(cv0, sv0);

						GivensRotateRight(X, p, q, cv, sv);
						GivensRotateLeft(V, p, q, cv, -sv);
					}
				}
			}
		}

		// Sort the singular values in decreasing order.
		std::array<std::pair<Real, size_t>, Columns> sortedS;
		for (size_t i = 0; i < Columns; ++i) {
			sortedS[i] = { Length(X.Column(i)), i };
		}
		std::sort(sortedS.begin(), sortedS.end(), [](const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });

		Vector<Real, Columns, Packed> S;
		for (size_t i = 0; i < Columns; ++i) {
			S(i) = sortedS[i].first;
		}


		/// Sort the rows/columns of X and V to match that of the sorted singular values.
		MatTall sortedX;
		MatSquare sortedV;
		for (size_t i = 0; i < Columns; ++i) {
			const auto from = sortedS[i].second;
			sortedX.Column(i, X.Column(from));
			sortedV.Row(i, V.Row(from));
		}
		X = sortedX;
		V = sortedV;

		// Normalize the columns of X using the singular values.
		const auto normalizationThreshold = S(0) * std::sqrt(std::numeric_limits<Real>::epsilon());
		size_t numOverThreshold = 0;
		for (size_t i = 0; i < Columns; ++i) {
			const auto scale = S(i);
			if (scale != Real(0)) {
				X.Column(i, X.Column(i) / scale);
			}
			numOverThreshold += static_cast<size_t>(scale >= normalizationThreshold);
		}

		// Do a QR decomposition if the singular values have a wide range, as the
		// columns associated with the smallest (or zero) singular values may not be orthogonal
		// to the rest due to numerical errors. The QR decomposition will orthogonalize them.
		// NOTE: worth examining doing a polar decomposition A = UP instead.
		if (numOverThreshold != Columns) {
			const auto [Q, R] = DecomposeQR(X);
			for (size_t i = 0; i < Columns; ++i) {
				S(i) = std::copysign(S(i), std::real(R(i, i)));
			}
			return DecompositionSVD{ Q, S * scaler, V };
		}
		else {
			return DecompositionSVD{ X, S * scaler, V };
		}
	}

} // namespace impl


enum class eSVDAlgorithm {
	JACOBI_ONE_SIDED,
	JACOBI_TWO_SIDED,
};


inline constexpr auto SVDAlgorithmOneSided = std::integral_constant<eSVDAlgorithm, eSVDAlgorithm::JACOBI_ONE_SIDED>{};
inline constexpr auto SVDAlgorithmTwoSided = std::integral_constant<eSVDAlgorithm, eSVDAlgorithm::JACOBI_TWO_SIDED>{};


/// <summary> Calculates the thin SVD of the matrix. </summary>
/// <remarks> For wide matrices, V is wide while U and S square.
///		For tall matrices, U is tall while S and V square. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, eSVDAlgorithm Algorithm = eSVDAlgorithm::JACOBI_TWO_SIDED>
auto DecomposeSVD(Matrix<T, Rows, Columns, Order, Layout, Packed> m,
				  std::integral_constant<eSVDAlgorithm, Algorithm> algorithm = {}) {
	if constexpr (Rows >= Columns) {
		if constexpr (Algorithm == eSVDAlgorithm::JACOBI_ONE_SIDED) {
			return impl::DecomposeSVDJacobiOneSided(m);
		}
		else {
			return impl::DecomposeSVDJacobiTwoSided(m);
		}
	}
	else {
		const auto [VT, S, UT] = DecomposeSVD(FlipLayoutAndOrder(m), algorithm);
		return DecompositionSVD{ FlipLayoutAndOrder(UT), S, FlipLayoutAndOrder(VT) };
	}
}

} // namespace mathter