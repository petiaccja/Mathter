// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Matrix.hpp"


namespace mathter {

template <int Rows, int Columns, class T1, class T2, eMatrixOrder Order, eMatrixLayout Layout1, eMatrixLayout Layout2, bool Packed1, bool Packed2>
bool operator==(const Matrix<T1, Rows, Columns, Order, Layout1, Packed1>& lhs, const Matrix<T2, Columns, Rows, opposite_order_v<Order>, Layout2, Packed2>& rhs) {
	bool equal = true;
	for (int i = 0; i < Rows; ++i) {
		for (int j = 0; j < Columns; ++j) {
			equal = equal && lhs(i, j) == rhs(j, i);
		}
	}
	return equal;
}


template <int Rows, int Columns, class T1, class T2, eMatrixOrder Order, eMatrixLayout Layout1, eMatrixLayout Layout2, bool Packed1, bool Packed2>
bool operator==(const Matrix<T1, Rows, Columns, Order, Layout1, Packed1>& lhs, const Matrix<T2, Rows, Columns, Order, Layout2, Packed2>& rhs) {
	bool equal = true;
	for (int i = 0; i < Rows; ++i) {
		for (int j = 0; j < Columns; ++j) {
			equal = equal && lhs(i, j) == rhs(i, j);
		}
	}
	return equal;
}


template <int Rows, int Columns, class T1, class T2, eMatrixOrder Order, eMatrixLayout Layout1, eMatrixLayout Layout2, bool Packed1, bool Packed2>
bool operator!=(const Matrix<T1, Rows, Columns, Order, Layout1, Packed1>& lhs, const Matrix<T2, Columns, Rows, opposite_order_v<Order>, Layout2, Packed2>& rhs) {
	return !(lhs == rhs);
}


template <int Rows, int Columns, class T1, class T2, eMatrixOrder Order, eMatrixLayout Layout1, eMatrixLayout Layout2, bool Packed1, bool Packed2>
bool operator!=(const Matrix<T1, Rows, Columns, Order, Layout1, Packed1>& lhs, const Matrix<T2, Rows, Columns, Order, Layout2, Packed2>& rhs) {
	return !(lhs == rhs);
}

} // namespace mathter