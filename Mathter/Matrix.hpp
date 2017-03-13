#pragma once

#include <type_traits>
#include <iostream> // debug only

#include "Vector.hpp"



enum class eMatrixLayout {
	ROW_MAJOR,
	COLUMN_MAJOR,
};



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
public:
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


template <class T, int Columns, int Rows, eMatrixOrder Order, bool Square = (Columns == Rows)>
class MatrixOps;

template <class T, int Columns, int Rows, eMatrixOrder Order>
class MatrixOps<T, Columns, Rows, Order, false> : public MatrixData<T, Columns, Rows, Order> {
protected:

};

template <class T, int Dim, eMatrixOrder Order>
class MatrixOps<T, Dim, Dim, Order, true> : public MatrixOps<T, Dim, Dim, Order, false> {
	using MyMatrixT = Matrix<T, Dim, Dim, Order>;
protected:
	template <class T2, eMatrixOrder Order2>
	MyMatrixT& operator*=(const Matrix<T2, Dim, Dim, Order2>& rhs) {
		auto& thisMat = static_cast<MyMatrixT&>(*this);
		thisMat = operator*<T, T2, Dim, Dim, Dim, Order, Order2, T>(thisMat, rhs);
		return thisMat;
	}

	float Trace() const;
	float Determinant() const;
	MyMatrixT& Transpose();
	MyMatrixT& Invert();
	MyMatrixT Inverted() const;
};


template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V = decltype(T() + U())>
Matrix<U, Columns, Rows, Order1> operator+(const Matrix<T, Columns, Rows, Order1>&, const Matrix<U, Columns, Rows, Order2>&);

template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V = decltype(T() + U())>
Matrix<U, Columns, Rows, Order1> operator-(const Matrix<T, Columns, Rows, Order1>&, const Matrix<U, Columns, Rows, Order2>&);



template <class T>
struct IsMatrix {
	static constexpr bool value = false;
};

template <class T, int Columns, int Rows, eMatrixOrder Order>
struct IsMatrix<Matrix<T, Columns, Rows, Order>> {
	static constexpr bool value = true;
};

template <class T>
struct NotMatrix {
	static constexpr bool value = !IsMatrix<T>::value;
};


template <class T>
struct IsScalar {
	static constexpr bool value = !IsMatrix<T>::value && !IsVector<T>::value;
};



template <class T, int Columns, int Rows, eMatrixOrder Order>
class Matrix : public MatrixOps<T, Columns, Rows, Order> {
protected:
	using MatrixData<T, Columns, Rows, Order>::GetElement;
public:
	//--------------------------------------------
	// Constructors
	//--------------------------------------------

	Matrix() = default;
	
	template <class... Args, typename std::enable_if<All<IsScalar, Args...>::value, int>::type = 0>
	Matrix(Args... args);


	//--------------------------------------------
	// Accessors
	//--------------------------------------------

	T& operator()(int col, int row) {
		return GetElement(col, row);
	}
	T operator()(int col, int row) const {
		return GetElement(col, row);
	}

	//--------------------------------------------
	// Arithmetic
	//--------------------------------------------

	template <class T, class U, eMatrixOrder Order2, class V>
	friend Matrix<U, Columns, Rows, Order> operator+(const Matrix<T, Columns, Rows, Order>&, const Matrix<U, Columns, Rows, Order2>&);

	template <class T, class U, eMatrixOrder Order2, class V>
	friend Matrix<U, Columns, Rows, Order> operator-(const Matrix<T, Columns, Rows, Order>&, const Matrix<U, Columns, Rows, Order2>&);


	//--------------------------------------------
	// Matrix functions
	//--------------------------------------------
	Matrix<T, Rows, Columns, Order> Transposed() const;


private:

};


template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V = MatMulElemT<T,U>, class = void>
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


template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V = decltype(T() + U())>
Matrix<U, Columns, Rows, Order1> operator+(const Matrix<T, Columns, Rows, Order1>& lhs, const Matrix<U, Columns, Rows, Order2>& rhs) {
	Matrix<U, Columns, Rows, Order1> result;
	for (int i = 0; i < Rows; ++i) {
		result.stripes[i] = lhs.stripes[i] + rhs.stripes[i];
	}
	return result;
}

template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, class V = decltype(T() + U())>
Matrix<U, Columns, Rows, Order1> operator-(const Matrix<T, Columns, Rows, Order1>& lhs, const Matrix<U, Columns, Rows, Order2>& rhs) {
	Matrix<U, Columns, Rows, Order1> result;
	for (int i = 0; i < Rows; ++i) {
		result.stripes[i] = lhs.stripes[i] - rhs.stripes[i];
	}
	return result;
}