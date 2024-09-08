// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/MathUtil.hpp"
#include "../Matrix/Algorithm.hpp"
#include "../Matrix/Arithmetic.hpp"
#include "../Matrix/Cast.hpp"
#include "../Matrix/Math.hpp"
#include "../Matrix/Matrix.hpp"
#include "../Transforms/IdentityBuilder.hpp"
#include "../Transforms/RandomBuilder.hpp"

#include <random>


namespace mathter {


/// <summary> The thin QR decomposition of a matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct DecompositionQR {
	static_assert(Rows >= Columns);

	using MatQ = Matrix<T, Rows, Columns, Order, Layout, Packed>;
	using MatR = Matrix<T, Columns, Columns, Order, Layout, Packed>;
	using Real = remove_complex_t<T>;

	MatQ Q;
	MatR R;

	/// <summary> Solve multiple linear systems of equations at the same time. </summary>
	/// <remarks> For overdetermined systems, it returns the least squares solution.
	///		Only applies to matrices that pre-multiply vectors, because the QR decomposition
	///		is only defined for tall matrices, and tall matrices represent least squares
	///		problems only for pre-multiplication. Use the LQ decomposition for
	///		post-multiplication and/or wide matrices. </remarks>
	template <class T2, int Columns2, eMatrixLayout Layout2, bool Packed2>
	auto Solve(const Matrix<T2, Rows, Columns2, Order, Layout2, Packed2>& b) const;

	/// <summary> Solve a linear systems of equations. </summary>
	/// <remarks> For overdetermined systems, it returns the least squares solution.
	///		Only applies to matrices that pre-multiply vectors, because the QR decomposition
	///		is only defined for tall matrices, and tall matrices represent least squares
	///		problems only for pre-multiplication. Use the LQ decomposition for
	///		post-multiplication and/or wide matrices. </remarks>
	template <class T2, bool Packed2>
	auto Solve(const Vector<T2, Rows, Packed2>& b) const;

	/// <summary> Compute the inverse or the pseudoinverse of the matrix. </summary>
	/// <remarks> For non-square matrices, the pseudoinverse is computed instead. </remarks>
	auto Inverse() const -> Matrix<T, Columns, Rows, Order, Layout, Packed>;

private:
	/// <summary> Solve multiple linear systems of equations at the same time. </summary>
	/// <remarks> This is a general routine that ignores matrix ordering.
	///		For internal use only. </remarks>
	template <class T2, int Columns2, eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2>
	auto SolveUnordered(const Matrix<T2, Rows, Columns2, Order2, Layout2, Packed2>& b) const;

#define MATHTER_QR_SOLVE_ORDER_ERROR "Solving equations with the QR decomposition is only applicable" \
									 " to matrices that precede vectors in multiplication."           \
									 " Use the LQ or QRorLQ decompositions to solve this system."
};


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
DecompositionQR(const Matrix<T, Rows, Columns, Order, Layout, Packed>&,
				const Matrix<T, Columns, Columns, Order, Layout, Packed>&) -> DecompositionQR<T, Rows, Columns, Order, Layout, Packed>;


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Columns2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionQR<T, Rows, Columns, Order, Layout, Packed>::Solve(const Matrix<T2, Rows, Columns2, Order, Layout2, Packed2>& b) const {
	static_assert(Order == eMatrixOrder::PRECEDE_VECTOR, MATHTER_QR_SOLVE_ORDER_ERROR);
	return SolveUnordered(b);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, bool Packed2>
auto DecompositionQR<T, Rows, Columns, Order, Layout, Packed>::Solve(const Vector<T2, Rows, Packed2>& b) const {
	static_assert(Order == eMatrixOrder::PRECEDE_VECTOR, MATHTER_QR_SOLVE_ORDER_ERROR);
	return Vector(SolveUnordered(Matrix<T2, Rows, 1, Order, Layout, Packed>(b)));
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecompositionQR<T, Rows, Columns, Order, Layout, Packed>::Inverse() const -> Matrix<T, Columns, Rows, Order, Layout, Packed> {
	const auto ordR = SetOrder<eMatrixOrder::PRECEDE_VECTOR>(R, std::false_type{});
	const auto ordI = SetOrder<eMatrixOrder::PRECEDE_VECTOR>(MatR(Identity()), std::false_type{});
	const auto Rinv = SetOrder<Order>(SolveUpperTriangular(ordR, ordI), std::false_type{});
	const auto Qinv = ConjTranspose(Q);
	return Rinv * Qinv;
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Columns2, eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionQR<T, Rows, Columns, Order, Layout, Packed>::SolveUnordered(const Matrix<T2, Rows, Columns2, Order2, Layout2, Packed2>& b) const {
	return SolveUpperTriangular(R, ConjTranspose(Q) * b);
}


/// <summary> The thin LQ decomposition of a matrix. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct DecompositionLQ {
	static_assert(Rows <= Columns);

	using MatL = Matrix<T, Rows, Rows, Order, Layout, Packed>;
	using MatQ = Matrix<T, Rows, Columns, Order, Layout, Packed>;
	using Real = remove_complex_t<T>;

	MatL L;
	MatQ Q;

	/// <summary> Solve multiple linear systems of equations at the same time. </summary>
	/// <remarks> For overdetermined systems, it returns the least squares solution.
	///		Only applies to matrices that post-multiply vectors, because the LQ decomposition
	///		is only defined for wide matrices, and wide matrices represent least squares
	///		problems only for post-multiplication. Use the QR decomposition for
	///		pre-multiplication and/or thin matrices. </remarks>
	template <class T2, int Rows2, eMatrixLayout Layout2, bool Packed2>
	auto Solve(const Matrix<T2, Rows2, Rows, Order, Layout2, Packed2>& b) const;

	/// <summary> Solve a linear systems of equations. </summary>
	/// <remarks> For overdetermined systems, it returns the least squares solution.
	///		Only applies to matrices that post-multiply vectors, because the LQ decomposition
	///		is only defined for wide matrices, and wide matrices represent least squares
	///		problems only for post-multiplication. Use the QR decomposition for
	///		pre-multiplication and/or thin matrices. </remarks>
	template <class T2, bool Packed2>
	auto Solve(const Vector<T2, Columns, Packed2>& b) const;

	/// <summary> Compute the inverse or the pseudoinverse of the matrix. </summary>
	/// <remarks> For non-square matrices, the pseudoinverse is computed instead. </remarks>
	auto Inverse() const -> Matrix<T, Columns, Rows, Order, Layout, Packed>;

#define MATHTER_LQ_SOLVE_ORDER_ERROR "Solving equations with the RQ decomposition is only applicable" \
									 " to matrices that follow vectors in multiplication."            \
									 " Use the QR or QRorLQ decompositions to solve this system."
};


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
DecompositionLQ(const Matrix<T, Rows, Rows, Order, Layout, Packed>&,
				const Matrix<T, Rows, Columns, Order, Layout, Packed>&) -> DecompositionLQ<T, Rows, Columns, Order, Layout, Packed>;


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Rows2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionLQ<T, Rows, Columns, Order, Layout, Packed>::Solve(const Matrix<T2, Rows2, Rows, Order, Layout2, Packed2>& b) const {
	static_assert(Order == eMatrixOrder::FOLLOW_VECTOR, MATHTER_LQ_SOLVE_ORDER_ERROR);
	return FlipLayoutAndOrder(DecompositionQR{ FlipLayoutAndOrder(Q), FlipLayoutAndOrder(L) }
								  .Solve(FlipLayoutAndOrder(b)));
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, bool Packed2>
auto DecompositionLQ<T, Rows, Columns, Order, Layout, Packed>::Solve(const Vector<T2, Columns, Packed2>& b) const {
	static_assert(Order == eMatrixOrder::FOLLOW_VECTOR, MATHTER_LQ_SOLVE_ORDER_ERROR);
	return DecompositionQR{ FlipLayoutAndOrder(Q), FlipLayoutAndOrder(L) }.Solve(b);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecompositionLQ<T, Rows, Columns, Order, Layout, Packed>::Inverse() const -> Matrix<T, Columns, Rows, Order, Layout, Packed> {
	return FlipLayoutAndOrder(DecompositionQR{ FlipLayoutAndOrder(Q), FlipLayoutAndOrder(L) }
								  .Inverse());
}


namespace impl {

	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto HouseholderReflection(const Matrix<T, Rows, Columns, Order, Layout, Packed>& r, size_t k) {
		using Real = remove_complex_t<T>;

		auto x = r.Column(k);
		std::fill_n(x.begin(), k, static_cast<T>(0));
		const auto normX = LengthPrecise(x);
		const auto pivot = x[k];
		const auto mag = std::abs(pivot);
		const auto sign = mag != Real(0) ? -(pivot / mag) : T(1);
		const auto alpha = sign * normX;

		auto u = x;
		u[k] -= alpha;

		const auto normU = LengthPrecise(u);

		if (normU != static_cast<Real>(0)) {
			return u / normU;
		}
		else {
			u[k] = static_cast<T>(1);
			return u;
		}
	}

} // namespace impl


/// <summary> Calculates the QR decomposition of the matrix. </summary>
/// <remarks> The QR decomposition is applicable to square or tall matrices (i.e. Rows >= Columns).
///		For wide matrices, use the LQ or QRorLQ decompositions. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecomposeQR(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	static_assert(Rows >= Columns);

	using ColumnMat = Matrix<T, Rows, 1, Order, Layout, Packed>;
	using RowMat = Matrix<T, 1, Rows, Order, Layout, Packed>;
	using MatQ = Matrix<T, Rows, Rows, Order, Layout, Packed>;
	using MatR = Matrix<T, Rows, Columns, Order, Layout, Packed>;
	using Real = remove_complex_t<T>;

	const auto scaler = std::max(std::numeric_limits<Real>::min(), ScaleElements(m));
	MatQ QT = Identity();
	MatR R = m / scaler;

	// Using Householder reflections.
	// Algorithm from Wikipedia's page: https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections
	// and Higham's book "Accuracy and Stability of Numerical Algorithms"
	for (size_t k = 0; k < m.ColumnCount(); ++k) {
		const auto v = impl::HouseholderReflection(R, k);
		const auto vm = ColumnMat(v);
		const auto vmT = RowMat(Conj(v));
		R = R - Real(2) * vm * (vmT * R);
		QT = QT - Real(2) * vm * (vmT * QT);
	}

	// Completely zero out lower triangle. These elements are near zero after householder transforms.
	for (size_t columnIdx = 0; columnIdx < Columns; ++columnIdx) {
		for (size_t rowIdx = columnIdx + 1; rowIdx < Columns; ++rowIdx) {
			R(rowIdx, columnIdx) = static_cast<T>(0);
		}
	}

	return DecompositionQR{
		ConjTranspose(QT.template Extract<Columns, Rows>(0, 0)),
		scaler * R.template Extract<Columns, Columns>(0, 0)
	};
}


/// <summary> Calculates the LQ decomposition of the matrix. </summary>
/// <remarks> The LQ decomposition is applicable to square or wide matrices (i.e. Rows <= Columns).
///		For tall matrices, use the QR or QRorLQ decompositions. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecomposeLQ(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	static_assert(Rows <= Columns);

	const auto [Q, R] = DecomposeQR(FlipLayoutAndOrder(m));
	return DecompositionLQ{ FlipLayoutAndOrder(R), FlipLayoutAndOrder(Q) };
}


/// <summary> Calculates the QR or LQ decomposition of the matrix. </summary>
/// <remarks> Chooses the QR or LQ decomposition for the provided matrix, depending
///		on which one is (more) suitable. Selection depends on the matrix's multiplication order
///		and shape. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecomposeQRorLQ(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	if constexpr (Rows == Columns) {
		if constexpr (Order == eMatrixOrder::PRECEDE_VECTOR) {
			return DecomposeQR(m);
		}
		else {
			return DecomposeLQ(m);
		}
	}
	else if constexpr (Rows > Columns) {
		return DecomposeQR(m);
	}
	else {
		return DecomposeLQ(m);
	}
}

} // namespace mathter