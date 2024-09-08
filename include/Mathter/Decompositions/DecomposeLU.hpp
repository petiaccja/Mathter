// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/Types.hpp"
#include "../Matrix/Algorithm.hpp"
#include "../Matrix/Matrix.hpp"
#include "../Transforms/IdentityBuilder.hpp"
#include "../Transforms/ZeroBuilder.hpp"

#include <numeric>


namespace mathter {


/// <summary> The LU decomposition of a matrix. </summary>
template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct DecompositionLU {
	using Mat = Matrix<T, Dim, Dim, Order, Layout, Packed>;

	Mat L;
	Mat U;

	/// <summary> Solve multiple linear systems of equations at the same time. </summary>
	template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
	auto Solve(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const;

	/// <summary> Solve a linear systems of equations. </summary>
	template <class T2, bool Packed2>
	auto Solve(const Vector<T2, Dim, Packed2>& b) const;

	/// <summary> Compute the inverse or the pseudoinverse of the matrix. </summary>
	auto Inverse() const -> Matrix<T, Dim, Dim, Order, Layout, Packed>;
};


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
DecompositionLU(const Matrix<T, Dim, Dim, Order, Layout, Packed>&,
				const Matrix<T, Dim, Dim, Order, Layout, Packed>&) -> DecompositionLU<T, Dim, Order, Layout, Packed>;


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionLU<T, Dim, Order, Layout, Packed>::Solve(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const {
	if constexpr (Order == eMatrixOrder::PRECEDE_VECTOR) {
		static_assert(Rows2 == Dim, "Incorrect shape for system of equations right-hand side.");
		const auto Ux = SolveLowerTriangular(L, b);
		const auto x = SolveUpperTriangular(U, Ux);
		return x;
	}
	else {
		static_assert(Columns2 == Dim, "Incorrect shape for system of equations right-hand side.");
		const auto xL = SolveUpperTriangular(U, b);
		const auto x = SolveLowerTriangular(L, xL);
		return x;
	}
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, bool Packed2>
auto DecompositionLU<T, Dim, Order, Layout, Packed>::Solve(const Vector<T2, Dim, Packed2>& b) const {
	using Vec = Vector<T2, Dim, Packed2>;
	constexpr auto vecMatRows = Order == eMatrixOrder::PRECEDE_VECTOR ? Dim : 1;
	constexpr auto vecMatCols = Order == eMatrixOrder::PRECEDE_VECTOR ? 1 : Dim;
	using VecMat = Matrix<T2, vecMatRows, vecMatCols, Order, Layout, Packed2>;
	return Vec(Solve(VecMat(b)));
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto mathter::DecompositionLU<T, Dim, Order, Layout, Packed>::Inverse() const -> Matrix<T, Dim, Dim, Order, Layout, Packed> {
	return Solve(Mat(Identity()));
}


/// <summary> The LUP decomposition of a matrix. </summary>
template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct DecompositionLUP {
	using Mat = Matrix<T, Dim, Dim, Order, Layout, Packed>;
	using Perm = Vector<uint32_t, Dim, Packed>;

	Mat L;
	Mat U;
	Perm P;

	/// <summary> Solve multiple linear systems of equations at the same time. </summary>
	template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
	auto Solve(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const;

	/// <summary> Solve a linear systems of equations. </summary>
	template <class T2, bool Packed2>
	auto Solve(const Vector<T2, Dim, Packed2>& b) const;

	/// <summary> Compute the inverse or the pseudoinverse of the matrix. </summary>
	auto Inverse() const -> Matrix<T, Dim, Dim, Order, Layout, Packed>;

	/// <summary> Expand the permutation vector to an orthonormal matrix. </summary>
	Mat ExpandPermutation() const;

	/// <summary> Expand the permutation vector to an orthonormal matrix. </summary>
	static Mat ExpandPermutation(Perm P);

	/// <summary> Applies the inverse permutation on the rows/columns of a matrix. </summary>
	/// <remarks> This can be used to permute x after solving L*U*P^-1*x = b. </remarks>
	template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
	auto InversePermute(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const;

	/// <summary> Applies the permutation on the rows/columns of a matrix. </summary>
	template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
	auto Permute(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const;
};


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, class Int>
DecompositionLUP(const Matrix<T, Dim, Dim, Order, Layout, Packed>&,
				 const Matrix<T, Dim, Dim, Order, Layout, Packed>&,
				 const Vector<Int, Dim, Packed>&) -> DecompositionLUP<T, Dim, Order, Layout, Packed>;


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionLUP<T, Dim, Order, Layout, Packed>::Solve(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const {
	if constexpr (Order == eMatrixOrder::PRECEDE_VECTOR) {
		return DecompositionLU{ L, U }.Solve(Permute(b));
	}
	else {
		return InversePermute(DecompositionLU{ L, U }.Solve(b));
	}
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, bool Packed2>
auto mathter::DecompositionLUP<T, Dim, Order, Layout, Packed>::Solve(const Vector<T2, Dim, Packed2>& b) const {
	using Vec = Vector<T2, Dim, Packed2>;
	constexpr auto vecMatRows = Order == eMatrixOrder::PRECEDE_VECTOR ? Dim : 1;
	constexpr auto vecMatCols = Order == eMatrixOrder::PRECEDE_VECTOR ? 1 : Dim;
	using VecMat = Matrix<T2, vecMatRows, vecMatCols, Order, Layout, Packed2>;
	return Vec(Solve(VecMat(b)));
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto mathter::DecompositionLUP<T, Dim, Order, Layout, Packed>::Inverse() const -> Matrix<T, Dim, Dim, Order, Layout, Packed> {
	return Solve(Mat(Identity()));
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecompositionLUP<T, Dim, Order, Layout, Packed>::ExpandPermutation() const -> Mat {
	return ExpandPermutation(P);
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecompositionLUP<T, Dim, Order, Layout, Packed>::ExpandPermutation(Perm P) -> Mat {
	Mat PM = Zero();
	for (size_t i = 0; i < Dim; ++i) {
		PM(i, P[i]) = static_cast<T>(1);
	}
	return PM;
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionLUP<T, Dim, Order, Layout, Packed>::InversePermute(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const {
	std::decay_t<decltype(b)> PinvB;
	for (size_t i = 0; i < Dim; ++i) {
		if constexpr (Order == eMatrixOrder::PRECEDE_VECTOR) {
			static_assert(Rows2 == Dim, "Incorrect shape for matrix to permute.");
			PinvB.Row(P[i], b.Row(i));
		}
		else {
			static_assert(Columns2 == Dim, "Incorrect shape for matrix to permute.");
			PinvB.Column(P[i], b.Column(i));
		}
	}
	return PinvB;
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, int Rows2, int Columns2, eMatrixLayout Layout2, bool Packed2>
auto DecompositionLUP<T, Dim, Order, Layout, Packed>::Permute(const Matrix<T2, Rows2, Columns2, Order, Layout2, Packed2>& b) const {
	std::decay_t<decltype(b)> PinvB;
	for (size_t i = 0; i < Dim; ++i) {
		if constexpr (Order == eMatrixOrder::PRECEDE_VECTOR) {
			static_assert(Rows2 == Dim, "Incorrect shape for matrix to permute.");
			PinvB.Row(i, b.Row(P[i]));
		}
		else {
			static_assert(Columns2 == Dim, "Incorrect shape for matrix to permute.");
			PinvB.Column(i, b.Column(P[i]));
		}
	}
	return PinvB;
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecomposeLU(const Matrix<T, Dim, Dim, Order, Layout, Packed>& m) {
	using Mat = std::decay_t<decltype(m)>;

	Mat L = Identity();
	Mat U = m;

	for (size_t zeroedColIdx = 0; zeroedColIdx < Dim - 1; ++zeroedColIdx) {
		const auto pivotRow = U.Row(zeroedColIdx);
		const auto pivot = pivotRow(zeroedColIdx);

		for (size_t zeroedRowIdx = zeroedColIdx + 1; zeroedRowIdx < Dim; ++zeroedRowIdx) {
			const auto target = U(zeroedRowIdx, zeroedColIdx);
			const auto scale = target / pivot;
			U.Row(zeroedRowIdx, U.Row(zeroedRowIdx) - pivotRow * scale);
			U(zeroedRowIdx, zeroedColIdx) = static_cast<T>(0); // Just to be sure that it is exactly zero.
			L(zeroedRowIdx, zeroedColIdx) = scale;
		}
	}

	return DecompositionLU{ L, U };
}


template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecomposeLUP(const Matrix<T, Dim, Dim, Order, Layout, Packed>& m) {
	using MatResult = Matrix<T, Dim, Dim, Order, Layout, Packed>;
	using MatWorking = Matrix<T, Dim, Dim, Order, eMatrixLayout::ROW_MAJOR, Packed>;
	using Perm = Vector<uint32_t, Dim, Packed>;
	using Real = remove_complex_t<T>;

	MatResult L = Identity();
	MatWorking U = m;
	Perm P;
	std::iota(P.begin(), P.end(), uint32_t(0));

	for (size_t zeroedColIdx = 0; zeroedColIdx < Dim - 1; ++zeroedColIdx) {
		const auto zeroedColAbs = Abs(U.Column(zeroedColIdx));
		const auto pivotRowIt = std::max_element(zeroedColAbs.begin() + zeroedColIdx, zeroedColAbs.end());
		const auto pivotRowIdx = pivotRowIt - zeroedColAbs.begin();
		const auto pivotRow = U.Row(pivotRowIdx);
		const auto pivot = pivotRow(zeroedColIdx);

		if (pivot != static_cast<Real>(0)) {
			// Swap pivot row and the current row.
			U.Row(pivotRowIdx, U.Row(zeroedColIdx));
			U.Row(zeroedColIdx, pivotRow);
			std::swap(P[pivotRowIdx], P[zeroedColIdx]);
			for (size_t i = 0; i < zeroedColIdx; ++i) {
				std::swap(L(pivotRowIdx, i), L(zeroedColIdx, i));
			}

			// Zero out column as usual business.
			for (size_t zeroedRowIdx = zeroedColIdx + 1; zeroedRowIdx < Dim; ++zeroedRowIdx) {
				const auto target = U(zeroedRowIdx, zeroedColIdx);
				const auto scale = target / pivot;
				U.Row(zeroedRowIdx, U.Row(zeroedRowIdx) - pivotRow * scale);
				U(zeroedRowIdx, zeroedColIdx) = static_cast<T>(0); // Just to be sure that it is exactly zero.
				L(zeroedRowIdx, zeroedColIdx) = scale;
			}
		}
	}

	return DecompositionLUP{ L, MatResult(U), P };
}

} // namespace mathter