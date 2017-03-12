#pragma once

#include <type_traits>
#include <iostream> // debug only

#include "Vector.hpp"

enum class eMatrixLayout {
	ROW_MAJOR,
	COLUMN_MAJOR,
};

enum class eMatrixOrder {
	PRECEDE_VECTOR,
	FOLLOW_VECTOR,
};

template <class T, int Columns, int Rows, eMatrixLayout Layout = eMatrixLayout::ROW_MAJOR, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR>
class Matrix;


template <class T, int Columns, int Rows, eMatrixLayout Layout = eMatrixLayout::ROW_MAJOR, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR>
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
	static constexpr int StripeDim = Layout == eMatrixLayout::ROW_MAJOR ? Columns : Rows;
	static constexpr int StripeCount = Layout == eMatrixLayout::ROW_MAJOR ? Rows : Columns;

	Vector<T, StripeDim> stripes[StripeCount];

	// Get element
	inline T& GetElement(int col, int row, std::true_type) {
		return stripes[row][col];
	}
	inline T& GetElement(int col, int row, std::false_type) {
		return stripes[col][row];
	}
	inline T GetElement(int col, int row, std::true_type) const {
		return stripes[row][col];
	}
	inline T GetElement(int col, int row, std::false_type) const {
		return stripes[col][row];
	}
	inline T& GetElement(int col, int row) {
		return GetElement(col, row, std::integral_constant<bool, Layout == eMatrixLayout::ROW_MAJOR>());
	}
	inline T GetElement(int col, int row) const {
		return GetElement(col, row, std::integral_constant<bool, Layout == eMatrixLayout::ROW_MAJOR>());
	}
};


template <class T, int Columns, int Rows, eMatrixLayout Layout = eMatrixLayout::ROW_MAJOR, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR>
class Matrix : public MatrixData<T, Columns, Rows, Layout, Order> {
protected:
	using MatrixData<T, Columns, Rows, Layout, Order>::GetElement;
public:

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

private:

};



template <class T,
		class U, int Match,
		int Rows1,
		int Columns2,
		eMatrixLayout Layout1,
		eMatrixLayout Layout2,
		eMatrixOrder Order1,
		eMatrixOrder Order2,
		typename std::enable_if<Layout1 == eMatrixLayout::ROW_MAJOR && Layout2 == eMatrixLayout::ROW_MAJOR && Columns2 <= 4, int>::type = 0>
auto operator*(const Matrix<T, Match, Rows1, Layout1, Order1>& lhs, const Matrix<U, Columns2, Match, Layout2, Order2>& rhs) {
	using ElemT = decltype(T() * U() + T() + U());
	using ResultT = Matrix<ElemT, Columns2, Rows1, eMatrixLayout::ROW_MAJOR, Order1>;
	using LeftT = Matrix<T, Match, Rows1, eMatrixLayout::ROW_MAJOR, Order1>;
	using RightT = Matrix<T, Columns2, Match, eMatrixLayout::ROW_MAJOR, Order2>;

	ResultT result;
	// With less than 4 sized vector ops, try to use SSE spread externally.
	if (Columns2 <= 4) {
		Vector<ElemT, Columns2> scalarMultiplier;
		for (int y = 0; y < Rows1; ++y) {
			scalarMultiplier.Spread(lhs(0, y));
			result.stripes[y] = scalarMultiplier * rhs.stripes[0];
		}
		for (int x = 1; x < Columns2; ++x) {
			for (int y = 0; y < Rows1; ++y) {
				scalarMultiplier.Spread(lhs(x, y));
				result.stripes[y] += scalarMultiplier * rhs.stripes[x];
			}
		}
	}
	// With more than 4 sized vector ops, SSE spread will be used by vectors internally.
	else {
		ElemT scalarMultiplier;
		for (int y = 0; y < Rows1; ++y) {
			scalarMultiplier = lhs(0, y);
			result.stripes[y] = scalarMultiplier * rhs.stripes[y];
		}
		for (int x = 1; x < Columns2; ++x) {
			for (int y = 0; y < Rows1; ++y) {
				scalarMultiplier = lhs(x, y);
				result.stripes[y] += scalarMultiplier * rhs.stripes[y];
			}
		}
	}
	return result;
}