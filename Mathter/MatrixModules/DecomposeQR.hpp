//==============================================================================
// This software is distributed under The Unlicense.
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma once

#include "../Common/MathUtil.hpp"


namespace mathter {


/// <summary> A utility class that can do common operations with the QR decomposition,
///		i.e. solving equation systems. </summary>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class DecompositionQR {
	using MatrixT = Matrix<T, Rows, Columns, Order, Layout, Packed>;

public:
	DecompositionQR(Matrix<T, Rows, Rows, Order, Layout, Packed> Q,
					Matrix<T, Rows, Columns, Order, Layout, Packed> R) : Q(Q), R(R) {}

	Matrix<T, Rows, Rows, Order, Layout, Packed> Q;
	Matrix<T, Rows, Columns, Order, Layout, Packed> R;
};


/// <summary> Calculates the QR decomposition of the matrix using Householder transforms. </summary>
/// <remarks> The matrix must have Rows &lt;= Columns. It's a full QR decomposition, not a thin one. </remarks>
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto DecomposeQR(Matrix<T, Rows, Columns, Order, Layout, Packed> m) {
	static_assert(Rows >= Columns);

	Matrix<T, Rows, Rows, Order, Layout, Packed> Q;
	Matrix<T, Rows, Columns, Order, Layout, Packed> R;

	R = m;
	Q.SetIdentity();

	Matrix<T, Rows, Rows, Order, Layout, Packed> Qi;
	Vector<T, Rows, Packed> u;
	Matrix<T, Rows, 1, Order, eMatrixLayout::ROW_MAJOR, Packed> v;

	for (int col = 0; col < m.ColumnCount(); ++col) {
		u = R.Column(col);
		for (int i = 0; i < col; ++i) {
			u(i) = T(0);
		}
		T alpha = sign(R(col, col)) * u.LengthPrecise();
		u(col) -= alpha;
		T norm = u.LengthPrecise();
		if (norm == 0) {
			continue;
		}
		u /= norm;
		v = u;
		Qi = (T(-2) * v) * v.Transposed();
		for (int i = 0; i < Q.ColumnCount(); ++i) {
			Qi(i, i) += T(1);
		}
		R = Qi * R;
		Q = Qi * Q;
	}
	Q = Q.Transposed();

	return DecompositionQR<T, Rows, Columns, Order, Layout, Packed>{ Q, R };
}


} // namespace mathter