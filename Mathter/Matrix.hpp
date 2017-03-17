#pragma once

#include <type_traits>
#include <iostream> // debug only

#include "Vector.hpp"

namespace mathter {


enum class eMatrixLayout {
	ROW_MAJOR,
	COLUMN_MAJOR,
};


//------------------------------------------------------------------------------
// Matrix base class only allocating the memory
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR>
class MatrixData {
public:
	constexpr int ColumnCount() const {
		return Columns;
	}
	constexpr int RowCount() const {
		return Rows;
	}
	constexpr int Width() const {
		return Columns;
	}
	constexpr int Height() const {
		return Rows;
	}
protected:
	// Rows equal height, Columns equal width, row-major has column-sized stripes
	static constexpr int StripeDim = Columns;
	static constexpr int StripeCount = Rows;

	Vector<T, StripeDim> stripes[StripeCount];

	// Get element
	inline T& GetElement(int col, int row) {
		return stripes[row][col];
	}

	inline T GetElement(int col, int row) const {
		return stripes[row][col];
	}
};


//------------------------------------------------------------------------------
// Matrix operations class extending MatrixData by special operations
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order, bool Square = (Columns == Rows)>
class MatrixOps;

template <class T, int Columns, int Rows, eMatrixOrder Order>
class MatrixOps<T, Columns, Rows, Order, false> : public MatrixData<T, Columns, Rows, Order> {
protected:

};

template <class T, int Dim, eMatrixOrder Order>
class MatrixOps<T, Dim, Dim, Order, true> : public MatrixOps<T, Dim, Dim, Order, false> {
	using MyMatrixT = Matrix<T, Dim, Dim, Order>;
public:
	template <class T2, eMatrixOrder Order2>
	MyMatrixT& operator*=(const Matrix<T2, Dim, Dim, Order2>& rhs) {
		auto& thisMat = static_cast<MyMatrixT&>(*this);
		thisMat = operator*<T, T2, Dim, Dim, Dim, Order, Order2, T>(thisMat, rhs);
		return thisMat;
	}

	T Trace() const;
	T Determinant() const;
	MyMatrixT& Transpose();
	MyMatrixT& Invert();
	MyMatrixT Inverted() const;

	static MyMatrixT Identity();
	MyMatrixT& SetIdentity();
};


//------------------------------------------------------------------------------
// Global Matrix function prototypes
//------------------------------------------------------------------------------

template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V = decltype(T() + U())>
Matrix<U, Columns, Rows, Order1> operator+(const Matrix<T, Columns, Rows, Order1>&, const Matrix<U, Columns, Rows, Order2>&);

template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V = decltype(T() + U())>
Matrix<U, Columns, Rows, Order1> operator-(const Matrix<T, Columns, Rows, Order1>&, const Matrix<U, Columns, Rows, Order2>&);


//------------------------------------------------------------------------------
// Matrix class providing the common interface for all matrices
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order>
class Matrix : public MatrixOps<T, Columns, Rows, Order> {
protected:
	using MatrixData<T, Columns, Rows, Order>::GetElement;
	using MatrixData::stripes;
public:
	//--------------------------------------------
	// Constructors
	//--------------------------------------------

	Matrix() = default;

	template <class H, class... Args, class = std::enable_if<impl::All<impl::IsScalar, Args...>::value, int>::type>
	Matrix(H h, Args... args) {
		static_assert(1 + sizeof...(Args) == Columns*Rows, "All elements of matrix have to be initialized.");
		Assign<0, 0>(h, args...);
	}


	//--------------------------------------------
	// Accessors
	//--------------------------------------------

	// General matrix indexing
	T& operator()(int col, int row) {
		return GetElement(col, row);
	}
	T operator()(int col, int row) const {
		return GetElement(col, row);
	}

	// Column and row vector simple indexing
	template <class = typename std::enable_if<(Columns == 1 && Rows > 1) || (Columns > 1 && Rows == 1)>::type>
	T& operator()(int idx) {
		return GetElement(Columns == 1 ? 0 : idx, Rows == 1 ? 0 : idx);
	}
	template <class = typename std::enable_if<(Columns == 1 && Rows > 1) || (Columns > 1 && Rows == 1)>::type>
	T operator()(int idx) const {
		return GetElement(Columns == 1 ? 0 : idx, Rows == 1 ? 0 : idx);
	}


	//--------------------------------------------
	// Arithmetic
	//--------------------------------------------

	// Non-modifying external operators
	template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V>
	friend Matrix<U, Columns, Rows, Order1> operator+(const Matrix<T, Columns, Rows, Order1>& lhs, const Matrix<U, Columns, Rows, Order2>& rhs);

	template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V>
	friend Matrix<U, Columns, Rows, Order1> operator-(const Matrix<T, Columns, Rows, Order1>& lhs, const Matrix<U, Columns, Rows, Order2>& rhs);

	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
		->Matrix<V, Columns2, Rows1, Order1>;

	// Internal modifying operators
	template <class U, eMatrixOrder Order2>
	Matrix& operator+=(const Matrix<U, Columns, Rows, Order2>& rhs) {
		*this = operator+<T, U, Columns, Rows, Order, Order2, T>(*this, rhs);
	}

	template <class U, eMatrixOrder Order2>
	Matrix& operator-=(const Matrix<U, Columns, Rows, Order2>& rhs) {
		*this = operator-<T, U, Columns, Rows, Order, Order2, T>(*this, rhs);
	}


	//--------------------------------------------
	// Matrix functions
	//--------------------------------------------
	auto Transposed() const -> Matrix<T, Rows, Columns, Order> {
		Matrix<T, Rows, Columns, Order> result;
		for (int y = 0; y < Rows; ++y) {
			for (int x = 0; x < Rows; ++x) {
				result(y, x) = (*this)(x, y);
			}
		}
		return result;
	}

	//--------------------------------------------
	// Matrix-vector arithmetic
	//--------------------------------------------
	template <class Vt, class Mt, int Vd, int Mcol, eMatrixOrder Morder, class Rt>
	friend Vector<Rt, Mcol> operator*(const Vector<Vt, Vd>& vec, const Matrix<Mt, Mcol, Vd, Morder>& mat);

	template <class Vt, class Mt, int Vd, int Mrow, eMatrixOrder Morder, class Rt>
	friend Vector<Rt, Mrow> operator*(const Matrix<Mt, Vd, Mrow, Morder>& mat, const Vector<Vt, Vd>& vec);
	

protected:
	//--------------------------------------------
	// Helpers
	//--------------------------------------------

	template <int x, int y, class Head, class... Args>
	void Assign(Head head, Args... args) {
		(*this)(x, y) = head;
		Assign<((x + 1) % Columns), ((x != Columns - 1) ? y : (y + 1))>(args...);
	}

	template <int, int>
	void Assign() {}
};



//------------------------------------------------------------------------------
// Matrix-Matrix arithmetic
//------------------------------------------------------------------------------

template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V = MatMulElemT<T, U>, class = void>
auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
-> Matrix<V, Columns2, Rows1, Order1>
{
	using ResultT = Matrix<V, Columns2, Rows1, Order1>;

	ResultT result;

	VectorSpec<V, Columns2> scalarMultiplier;
	for (int y = 0; y < Rows1; ++y) {
		scalarMultiplier = rhs.stripes[0];
		scalarMultiplier.mul(lhs(0, y));
		static_cast<VectorSpec<V, Columns2>&>(result.stripes[y]) = scalarMultiplier;
	}
	for (int x = 1; x < Match; ++x) {
		for (int y = 0; y < Rows1; ++y) {
			scalarMultiplier = rhs.stripes[x];
			scalarMultiplier.mul(lhs(x, y));
			static_cast<VectorSpec<V, Columns2>&>(result.stripes[y]).add(scalarMultiplier);
		}
	}

	return result;
}


template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V>
Matrix<U, Columns, Rows, Order1> operator+(const Matrix<T, Columns, Rows, Order1>& lhs, const Matrix<U, Columns, Rows, Order2>& rhs) {
	Matrix<U, Columns, Rows, Order1> result;
	for (int i = 0; i < Rows; ++i) {
		result.stripes[i] = lhs.stripes[i] + rhs.stripes[i];
	}
	return result;
}

template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V>
Matrix<U, Columns, Rows, Order1> operator-(const Matrix<T, Columns, Rows, Order1>& lhs, const Matrix<U, Columns, Rows, Order2>& rhs) {
	Matrix<U, Columns, Rows, Order1> result;
	for (int i = 0; i < Rows; ++i) {
		result.stripes[i] = lhs.stripes[i] - rhs.stripes[i];
	}
	return result;
}


//------------------------------------------------------------------------------
// MatrixOps : Square matrices
//------------------------------------------------------------------------------

template <class T, int Dim, eMatrixOrder Order>
T MatrixOps<T, Dim, Dim, Order, true>::Trace() const {
	T sum = GetElement(0, 0);
	for (int i = 1; i < Dim; ++i) {
		sum += GetElement(i, i);
	}
	return sum;
}

template <class T, int Dim, eMatrixOrder Order>
T MatrixOps<T, Dim, Dim, Order, true>::Determinant() const {
	static_assert(false, "Determinant not implemented yet.");
	return T();
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::Transpose() -> MyMatrixT& {
	static_cast<MyMatrixT&>(*this) = static_cast<MyMatrixT&>(*this).Transposed();
	return static_cast<MyMatrixT&>(*this);
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::Invert() -> MyMatrixT& {
	static_assert(false, "Inverse not implemented yet.");
	return static_cast<MyMatrixT&>(*this);
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::Inverted() const -> MyMatrixT {
	static_assert(false, "Inverse not implemented yet.");
	auto copy = return static_cast<MyMatrixT&>(*this);
	return copy.Invert();
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::Identity() -> MyMatrixT {
	MyMatrixT res;

	for (int i = 0; i < Dim; ++i) {
		static_cast<MatrixOps&>(res).stripes[i].Spread(T(0));
		res(i, i) = T(1);
	}

	return res;
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::SetIdentity() -> MyMatrixT& {
	static_cast<MyMatrixT&>(*this) = Identity();
}



//------------------------------------------------------------------------------
// Matrix-vector arithmetic
//------------------------------------------------------------------------------

// v*M
template <class Vt, class Mt, int Vd, int Mcol, eMatrixOrder Morder, class Rt = MatMulElemT<Vt, Mt>>
Vector<Rt, Mcol> operator*(const Vector<Vt, Vd>& vec, const Matrix<Mt, Mcol, Vd, Morder>& mat) {
	Vector<Rt, Mcol> result;
	result = vec(0) * mat.stripes[0];
	for (int i = 1; i < Vd; ++i) {
		result += vec(i) * mat.stripes[i];
	}
	return result;
}

// M*v
template <class Vt, class Mt, int Vd, int Mrow, eMatrixOrder Morder, class Rt = MatMulElemT<Vt, Mt>>
Vector<Rt, Mrow> operator*(const Matrix<Mt, Vd, Mrow, Morder>& mat, const Vector<Vt, Vd>& vec) {
	Vector<Rt, Mrow> result;
	for (int i = 0; i < Mrow; ++i) {
		result(i) = vec.Dot(vec, mat.stripes[i]);
	}
	return result;
}

// v*=M
template <class Vt, class Mt, int Vd, eMatrixOrder Morder>
Vector<Vt, Vd>& operator*=(Vector<Vt, Vd>& vec, const Matrix<Mt, Vd, Vd, Morder>& mat) {
	vec = operator*<Vt, Mt, Vd, Vd, Morder, Vt>(vec, mat);
	return vec;
}



} // namespace mathter
