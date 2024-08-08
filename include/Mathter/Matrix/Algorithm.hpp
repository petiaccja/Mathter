// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Matrix.hpp"


namespace mathter {

template <class T1, int Dim, eMatrixOrder Order1, eMatrixLayout Layout1, bool Packed1,
		  class T2, int Columns2, eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2>
auto SolveUpperTriangular(const Matrix<T1, Dim, Dim, Order1, Layout1, Packed1>& A,
						  const Matrix<T2, Dim, Columns2, Order2, Layout2, Packed2>& b) {
	Matrix<T1, Dim, Dim, Order1, eMatrixLayout::ROW_MAJOR, Packed1> AWorking(A);
	Matrix<T2, Dim, Columns2, Order2, eMatrixLayout::ROW_MAJOR, Packed2> bWorking(b);

	for (ptrdiff_t rowIdx = Dim - 1; rowIdx >= 0; --rowIdx) {
		const auto& leading = AWorking(rowIdx, rowIdx);
		const auto rowANorm = AWorking.Row(rowIdx) / leading;
		const auto rowBNorm = bWorking.Row(rowIdx) / leading;

		AWorking.Row(rowIdx, rowANorm);
		bWorking.Row(rowIdx, rowBNorm);

		for (ptrdiff_t targetIdx = rowIdx - 1; targetIdx >= 0; --targetIdx) {
			const auto targetA = AWorking.Row(targetIdx);
			const auto targetB = bWorking.Row(targetIdx);
			const auto multiplier = targetA[rowIdx];
			AWorking.Row(targetIdx, targetA - multiplier * rowANorm);
			bWorking.Row(targetIdx, targetB - multiplier * rowBNorm);
		}
	}

	return std::decay_t<decltype(b)>(bWorking);
}

} // namespace mathter
