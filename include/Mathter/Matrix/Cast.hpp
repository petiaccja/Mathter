// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Matrix.hpp"


namespace mathter {


template <class Mat>
struct flip_layout_and_order;


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct flip_layout_and_order<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	using type = Matrix<T, Columns, Rows, opposite_order_v<Order>, opposite_layout_v<Layout>, Packed>;
};


template <class Mat>
using flip_layout_and_order_t = typename flip_layout_and_order<Mat>::type;


/// <summary> Flip both the order and the layout of the matrix. </summary>
/// <remarks> This essentially transposes the matrix, but leaves the stripes the same.
///		Should be optimized out by compilers. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto FlipLayoutAndOrder(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) {
	return flip_layout_and_order_t<std::decay_t<decltype(m)>>(m);
}


/// <summary> Flip  the order of the matrix. </summary>
/// <typeparam name="PreserveTransform"> If true, transposes the matrix to preserve transform,
///		if false, preserves stripes and is essentially a noop. </typeparam>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, bool PreserveTransform>
auto FlipOrder(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m, std::integral_constant<bool, PreserveTransform>) {
	using Mat = Matrix<T, Columns, Rows, opposite_order_v<Order>, Layout, Packed>;
	if constexpr (PreserveTransform) {
		return Mat(m);
	}
	else {
		Mat out;
		out.stripes = m.stripes;
		return out;
	}
}


/// <summary> Flip  the order of the matrix. </summary>
/// <typeparam name="PreserveTransform"> If true, transposes the matrix to preserve transform,
///		if false, preserves stripes and is essentially a noop. </typeparam>
template <eMatrixOrder DesiredOrder, class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, bool PreserveTransform>
auto SetOrder(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m, std::integral_constant<bool, PreserveTransform>) {
	using Mat = Matrix<T, Columns, Rows, DesiredOrder, Layout, Packed>;
	if constexpr (PreserveTransform) {
		return Mat(m);
	}
	else {
		Mat out;
		out.stripes = m.stripes;
		return out;
	}
}

} // namespace mathter