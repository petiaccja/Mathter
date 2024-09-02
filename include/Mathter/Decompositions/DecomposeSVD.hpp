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


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct DecompositionSVD {
	using Real = remove_complex_t<T>;
	static constexpr int PDim = std::min(Rows, Columns);

	Matrix<T, Rows, PDim, Order, Layout, Packed> U;
	Vector<Real, PDim, Packed> S;
	Matrix<T, PDim, Columns, Order, Layout, Packed> V;
};


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, class = std::enable_if_t<(Rows >= Columns)>>
DecompositionSVD(const Matrix<T, Rows, Columns, Order, Layout, Packed>&,
				 const Vector<remove_complex_t<T>, Columns, Packed>&,
				 const Matrix<T, Columns, Columns, Order, Layout, Packed>&) -> DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>;


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, class = std::enable_if_t<(Rows < Columns)>>
DecompositionSVD(const Matrix<T, Rows, Rows, Order, Layout, Packed>&,
				 const Vector<remove_complex_t<T>, Rows, Packed>&,
				 const Matrix<T, Rows, Columns, Order, Layout, Packed>&) -> DecompositionSVD<T, Rows, Columns, Order, Layout, Packed>;


namespace impl {

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


	/// <summary> Find the Givens rotation with the smallest angle. </summary>
	/// <remarks> The 2x2 SVD and the symmatric diagonalizer both have 4 solutions, offset by
	///		90 degrees. The one of the 4 that has the smallest rotation angle (i.e. largest cosine)
	///		may be the fastest to converge on Jacobi iterations because it applies the least modification
	///		to the rows/columns updated by the Givens rotation while still zeroing out the
	///		required elements. </remarks>
	template <class T>
	std::tuple<T, T> GetMinimalRotation(T cv, T sv) {
		using Real = remove_complex_t<T>;
		if (std::abs(cv) < std::abs(sv)) {
			std::tie(cv, sv) = std::tuple(-conj{}(sv), conj{}(cv)); // Rotate by 90 degrees to make cv bigger in absolute value.
		}
		if (std::real(cv) < static_cast<Real>(0)) {
			std::tie(cv, sv) = std::tuple(-cv, -sv); // Rotate by 180 degrees to make cv positive.
		}
		return { cv, sv };
	}


	/// <summary> Recompute the sine so that cv^2 + sv^2 = 1 holds precisely. </summary>
	template <class T>
	std::tuple<T, T> NormalizeRotation(T cv, T sv) {
		using Real = remove_complex_t<T>;
		// Assuming that cv^2 + sv^2 = 1 - epsilon, where epsilon << 1.
		// Instead of doing the full 1 / sqrt(cv^2 + sv^2), uses a second-order
		// series expansion of 1 / sqrt(1 - epsilon), since epsilon should be super small.
		const auto cvMag = std::abs(cv);
		const auto svMag = std::abs(sv);
		const auto epsilon = std::fma(-cvMag, cvMag, Real(1)) - svMag * svMag;
		const auto invNorm = Real(0.375) * epsilon * epsilon + Real(0.5) * epsilon + T(1); // Accumulate smallest first!
		return { cv * invNorm, sv * invNorm };
	}


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


		const auto scale = std::sqrt((std::real(cv) * std::real(cv)
									  + std::imag(cv) * std::imag(cv))
									 + (std::real(sv) * std::real(sv)
										+ std::imag(sv) * std::imag(sv)));

		std::tie(cv, sv) = std::tuple(cv / scale, sv / scale);
		std::tie(cv, sv) = GetMinimalRotation(cv, sv);

		return { cv, sv };
	}


	template <class T, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	DecompositionRQ2x2<T> DecomposeRQ2x2(const Matrix<T, 2, 2, Order, Layout, Packed>& A) {
		using Real = remove_complex_t<T>;

		const auto a11 = A(0, 0);
		const auto a12 = A(0, 1);
		const auto a21 = A(1, 0);
		const auto a22 = A(1, 1);

		// If a21 is not precisely zero, the scaling will take care of it and the rest works fine.
		if (a21 == T(0)) {
			return { a11, a12, a22, T(1), T(0) };
		}

		// Rescale matrix elements to avoid underflow and overflow.
		const auto scaleNum = ScaleElements(a11, a12, T(std::numeric_limits<Real>::min())); // Avoid NaN.
		const auto scaleDiv = ScaleElements(a21, a22); // a21 is never zero. See if statement above.
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
		const auto scaler = ScaleElements(a11, a12, a22);
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

		auto cv = p1 != 0 ? std::clamp(-pab * rsqrt2, -T(1), T(1)) : T(1);
		auto sv = p1 != 0 ? std::clamp(zsign * sqrt2 * g / (p1 * pab), -T(1), T(1)) : T(0);
		std::tie(cv, sv) = GetMinimalRotation(cv, sv);
		std::tie(cv, sv) = NormalizeRotation(cv, sv);

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


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	void GivensRotateRight(Matrix<T, Rows, Columns, Order, Layout, Packed>& m, int p, int q, T cv, T sv) {
		const auto colP = m.Column(p) * cv + m.Column(q) * sv;
		const auto colQ = m.Column(q) * conj{}(cv)-m.Column(p) * conj{}(sv);
		m.Column(p, colP);
		m.Column(q, colQ);
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	void GivensRotateLeft(Matrix<T, Rows, Columns, Order, Layout, Packed>& m, int p, int q, T cv, T sv) {
		const auto rowP = m.Row(p) * cv - m.Row(q) * conj{}(sv);
		const auto rowQ = m.Row(q) * conj{}(cv) + m.Row(p) * sv;
		m.Row(p, rowP);
		m.Row(q, rowQ);
	}


	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto DecomposeSVDJacobiTwoSided(const Matrix<T, Rows, Columns, Order, Layout, Packed>& A) {
		using Real = remove_complex_t<T>;
		using MatSquare = Matrix<T, Rows, Columns, Order, Layout, Packed>;
		using MatTall = Matrix<T, Rows, Rows, Order, Layout, Packed>;
		constexpr auto tolerance = Real(1) * std::numeric_limits<Real>::epsilon();

		const auto scaler = ScaleElements(A);
		MatTall U = Identity();
		MatSquare X = A / scaler;
		MatSquare V = Identity();

		auto maxErrorPrev = std::numeric_limits<Real>::max();
		auto maxError = std::nextafter(maxErrorPrev, Real(0));
		while (maxError < maxErrorPrev) {
			maxErrorPrev = maxError;
			maxError = Real(0);
			for (int p = 0; p < Columns; ++p) {
				for (int q = p + 1; q < Columns; ++q) {
					const auto [xpp, xpq, xqp, xqq] = std::tie(X(p, p), X(p, q), X(q, p), X(q, q));
					const auto error = std::max(std::abs(xpq), std::abs(xqp));
					if (error > tolerance) {
						maxError = std::max(maxError, error);
						const auto svd2x2 = DecomposeSVD2x2(Matrix<T, 2, 2>{ xpp, xpq, xqp, xqq });

						GivensRotateRight(U, p, q, svd2x2.cu, svd2x2.su);
						GivensRotateLeft(X, p, q, svd2x2.cu, -svd2x2.su);
						GivensRotateRight(X, p, q, svd2x2.cv, -svd2x2.sv);
						GivensRotateLeft(V, p, q, svd2x2.cv, svd2x2.sv);
					}
				}
			}
		}

		Vector<Real, Columns, Packed> S;
		for (size_t i = 0; i < Columns; ++i) {
			S(i) = std::real(X(i, i));
		}

		return DecompositionSVD{ U, S * scaler, V };
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
		using Real = remove_complex_t<T>;

		using MatSquare = Matrix<T, Rows, Columns, Order, Layout, Packed>;
		using MatTall = Matrix<T, Rows, Rows, Order, Layout, Packed>;

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
						const auto [cv, sv] = DiagonalizeHermitian2x2(ata11, ataoff, ata22);

						GivensRotateRight(X, p, q, cv, sv);

						const auto [_ignore1, postError, _ignore2] = TransposeMultiplyPartial(X, p, q);

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