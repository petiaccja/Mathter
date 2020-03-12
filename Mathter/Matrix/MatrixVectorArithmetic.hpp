#pragma once

#include "MatrixImpl.hpp"

namespace mathter {


//------------------------------------------------------------------------------
// Matrix-vector arithmetic
//------------------------------------------------------------------------------

// v*M
template <class Vt, class Mt, int Vd, int Mcol, eMatrixOrder Morder, bool Packed, class Rt>
inline Vector<Rt, Mcol, Packed> operator*(const Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Vd, Mcol, Morder, eMatrixLayout::ROW_MAJOR, Packed>& mat) {
	Vector<Rt, Mcol, Packed> result;

	result = vec(0) * mat.stripes[0];
	for (int i = 1; i < Vd; ++i) {
		result += vec(i) * mat.stripes[i];
	}
	return result;
}

template <class Vt, class Mt, int Vd, int Mcol, eMatrixOrder Morder, bool Packed, class Rt>
inline Vector<Rt, Mcol, Packed> operator*(const Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Vd, Mcol, Morder, eMatrixLayout::COLUMN_MAJOR, Packed>& mat) {
	Vector<Rt, Mcol, Packed> result;

	for (int i = 0; i < Vd; ++i) {
		result(i) = Dot(vec, mat.stripes[i]);
	}
	return result;
}

template <class Vt, class Mt, int Vd, int Mcol, eMatrixOrder Morder, eMatrixLayout Mlayout, bool Packed>
inline Vector<traits::MatMulElemT<Vt, Mt>, Mcol, Packed> operator*(const Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Vd, Mcol, Morder, Mlayout, Packed>& mat) {
	using Rt = traits::MatMulElemT<Vt, Mt>;
	return operator*<Vt, Mt, Vd, Mcol, Morder, Packed, Rt>(vec, mat);
}


// (v|1)*M
template <class Vt, class Mt, int Vd, eMatrixLayout Mlayout, eMatrixOrder Morder, bool Packed, class Rt = traits::MatMulElemT<Vt, Mt>>
Vector<Rt, Vd, Packed> operator*(const Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Vd + 1, Vd, Morder, Mlayout, Packed>& mat) {
	return (vec | Vt(1)) * mat;
}

template <class Vt, class Mt, int Vd, eMatrixLayout Mlayout, eMatrixOrder Morder, bool Packed, class Rt = traits::MatMulElemT<Vt, Mt>>
Vector<Rt, Vd, Packed> operator*(const Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Vd + 1, Vd + 1, Morder, Mlayout, Packed>& mat) {
	auto res = (vec | Vt(1)) * mat;
	res /= res(res.Dimension() - 1);
	return (Vector<Rt, Vd, Packed>)res;
}

// M*v
template <class Vt, class Mt, int Vd, int Mrow, eMatrixOrder Morder, bool Packed, class Rt>
inline Vector<Rt, Mrow, Packed> operator*(const Matrix<Mt, Mrow, Vd, Morder, eMatrixLayout::ROW_MAJOR, Packed>& mat, const Vector<Vt, Vd, Packed>& vec) {
	Vector<Rt, Mrow, Packed> result;
	
	for (int i = 0; i < Mrow; ++i) {
		result(i) = Dot(vec, mat.stripes[i]);
	}
	return result;
}

template <class Vt, class Mt, int Vd, int Mrow, eMatrixOrder Morder, bool Packed, class Rt>
inline Vector<Rt, Mrow, Packed> operator*(const Matrix<Mt, Mrow, Vd, Morder, eMatrixLayout::COLUMN_MAJOR, Packed>& mat, const Vector<Vt, Vd, Packed>& vec) {
	Vector<Rt, Mrow, Packed> result;

	result = vec(0) * mat.stripes[0];
	for (int i = 1; i < Vd; ++i) {
		result += vec(i) * mat.stripes[i];
	}
	return result;
}

template <class Vt, class Mt, int Vd, int Mrow, eMatrixOrder Morder, eMatrixLayout Mlayout, bool Packed>
Vector<traits::MatMulElemT<Vt, Mt>, Mrow, Packed> operator*(const Matrix<Mt, Mrow, Vd, Morder, Mlayout, Packed>& mat, const Vector<Vt, Vd, Packed>& vec) {
	using Rt = traits::MatMulElemT<Vt, Mt>;
	return operator*<Vt, Mt, Vd, Mrow, Morder, Packed, Rt>(mat, vec);
}



// M*(v|1)
template <class Vt, class Mt, int Vd, eMatrixLayout Mlayout, eMatrixOrder Morder, bool Packed, class Rt = traits::MatMulElemT<Vt, Mt>>
Vector<Rt, Vd, Packed> operator*(const Matrix<Mt, Vd, Vd + 1, Morder, Mlayout, Packed>& mat, const Vector<Vt, Vd, Packed>& vec) {
	return mat * (vec | Vt(1));
}

template <class Vt, class Mt, int Vd, eMatrixLayout Mlayout, eMatrixOrder Morder, bool Packed, class Rt = traits::MatMulElemT<Vt, Mt>>
Vector<Rt, Vd, Packed> operator*(const Matrix<Mt, Vd + 1, Vd + 1, Morder, Mlayout, Packed>& mat, const Vector<Vt, Vd, Packed>& vec) {
	auto res = (vec | Vt(1)) * mat;
	res /= res(res.Dimension() - 1);
	return (Vector<Rt, Vd, Packed>)res;
}

// v*=M
template <class Vt, class Mt, int Vd, eMatrixOrder Morder, eMatrixLayout Layout, bool Packed>
Vector<Vt, Vd, Packed>& operator*=(Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Vd, Vd, Morder, Layout, Packed>& mat) {
	vec = operator*<Vt, Mt, Vd, Vd, Morder, Packed, Vt>(vec, mat);
	return vec;
}


} // namespace mathter