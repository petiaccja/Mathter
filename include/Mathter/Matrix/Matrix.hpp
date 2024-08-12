// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once


#include "../Common/TypeTraits.hpp"
#include "../Common/Types.hpp"
#include "../Vector/Vector.hpp"

#include <array>
#include <cassert>
#include <type_traits>


namespace mathter {


//------------------------------------------------------------------------------
// Matrix base class only allocating the memory
//------------------------------------------------------------------------------

template <class T, int Rows, int Columns, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout Layout = eMatrixLayout::ROW_MAJOR, bool Packed = false>
struct MatrixStorage {
	static constexpr int stripeDim = Layout == eMatrixLayout::ROW_MAJOR ? Columns : Rows;
	static constexpr int stripeCount = Layout == eMatrixLayout::ROW_MAJOR ? Rows : Columns;

	using Stripe = Vector<T, stripeDim, Packed>;
	std::array<Stripe, stripeCount> stripes;
};

struct StripeArgType {};
constexpr StripeArgType stripeArg;


//------------------------------------------------------------------------------
// Matrix class providing the common interface for all matrices
//------------------------------------------------------------------------------

template <class T, int Rows, int Columns, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout Layout = eMatrixLayout::ROW_MAJOR, bool Packed = false>
class Matrix : public MatrixStorage<T, Rows, Columns, Order, Layout, Packed> {
	static_assert(Columns >= 1 && Rows >= 1, "Dimensions must be positive integers.");

public:
	using MatrixStorage<T, Rows, Columns, Order, Layout, Packed>::stripes;
	using typename MatrixStorage<T, Rows, Columns, Order, Layout, Packed>::Stripe;
	using MatrixStorage<T, Rows, Columns, Order, Layout, Packed>::stripeDim;
	using MatrixStorage<T, Rows, Columns, Order, Layout, Packed>::stripeCount;
	static constexpr bool isVector = std::min(Rows, Columns) == 1;
	static constexpr int vectorDim = isVector ? std::max(Rows, Columns) : 0;

public:
	//--------------------------------------------
	// Constructors
	//--------------------------------------------

	Matrix() = default;

	/// <summary> Copies and converts the element of <paramref name="rhs"/>. </summary>
	template <class T2, eMatrixLayout Layout2, bool Packed2>
	Matrix(const Matrix<T2, Rows, Columns, Order, Layout2, Packed2>& rhs);

	/// <summary> Copies and converts the element of <paramref name="rhs"/>. </summary>
	template <class T2, eMatrixLayout Layout2, bool Packed2>
	Matrix(const Matrix<T2, Columns, Rows, opposite_order_v<Order>, Layout2, Packed2>& rhs);

	/// <summary> Creates a matrix from its elements. </summary>
	template <class... Scalars, class = std::enable_if_t<(sizeof...(Scalars) == size_t(Rows * Columns)) && (... && std::is_convertible_v<std::decay_t<Scalars>, T>), int>>
	Matrix(Scalars&&... elements);

	/// <summary> Constructs a row or column matrix from a vector. </summary>
	/// <remarks> This can only be used for row or column matrices. </remarks>
	template <class TOther, bool PackedOther, class = std::enable_if_t<isVector && sizeof(TOther) != 0>>
	explicit Matrix(const Vector<TOther, vectorDim, PackedOther>& v);

	/// <summary> Used by internal methods. </summary>
	template <class... Stripes>
	Matrix(StripeArgType, Stripes&&... stripes)
		: MatrixStorage<T, Rows, Columns, Order, Layout, Packed>{ std::forward<Stripes>(stripes)... } {}


	//--------------------------------------------
	// Properties
	//--------------------------------------------

	/// <summary> Returns the number of columns of the matrix. </summary>
	constexpr size_t ColumnCount() const;

	/// <summary> Returns the number of rows of the matrix. </summary>
	constexpr size_t RowCount() const;

	/// <summary> Returns the number of columns of the matrix. </summary>
	constexpr size_t Width() const;

	/// <summary> Returns the number of rows of the matrix. </summary>
	constexpr size_t Height() const;

	//--------------------------------------------
	// Accessors
	//--------------------------------------------

	/// <summary> Returns the <paramref name="col"/>-th element of the <paramref name="row"/>-th row. </summary>
	T& operator()(size_t row, size_t col);

	/// <summary> Returns the <paramref name="col"/>-th element of the <paramref name="row"/>-th row. </summary>
	const T& operator()(size_t row, size_t col) const;

	/// <summary> Returns the <paramref name="idx"/>-th element of a row or column matrix. </summary>
	template <class TAlias = T>
	std::enable_if_t<isVector && std::is_same_v<T, TAlias>, TAlias>& operator()(size_t idx);

	/// <summary> Returns the <paramref name="idx"/>-th element of a row or column matrix. </summary>
	template <class TAlias = T>
	const std::enable_if_t<isVector && std::is_same_v<T, TAlias>, TAlias>& operator()(size_t idx) const;

	/// <summary> Return the matrix's column as a vector. </summary>
	auto Column(size_t colIdx) const;

	/// <summary> Return the matrix's row as a vector. </summary>
	auto Row(size_t rowIdx) const;

	/// <summary> Return the matrix's column as a vector. </summary>
	void Column(size_t colIdx, const Vector<T, Rows, Packed>& column);

	/// <summary> Return the matrix's row as a vector. </summary>
	void Row(size_t rowIdx, const Vector<T, Columns, Packed>& row);

	/// <summary> Extract a submatrix. </summary>
	/// <typeparam name="ExtractedRows"> Number of rows of the submatrix. </typeparam>
	/// <typeparam name="ExtractedColumns"> Number of columns of the submatrix. </typeparam>
	/// <param name="whereRow"> First row to extract. </param>
	/// <param name="whereCol"> First columns to extract. </param>
	template <int ExtractedRows, int ExtractedColumns>
	Matrix<T, ExtractedRows, ExtractedColumns, Order, Layout, Packed> Extract(int whereRow, int whereCol) const;

	/// <summary> Insert a submatrix. </summary>
	/// <typeparam name="ExtractedRows"> Number of rows of the submatrix. </typeparam>
	/// <typeparam name="ExtractedColumns"> Number of columns of the submatrix. </typeparam>
	/// <param name="whereRow"> First row to extract. </param>
	/// <param name="whereCol"> First columns to extract. </param>
	template <int InsertedRows, int InsertedColumns>
	void Insert(int whereRow, int whereCol, const Matrix<T, InsertedRows, InsertedColumns, Order, Layout, Packed>& submatrix);

	/// <summary> Convert a row or column matrix to a vector. </summary>
	template <class TOther, bool PackedOther, class = std::enable_if_t<isVector && std::is_convertible_v<T, TOther>>>
	explicit operator Vector<TOther, vectorDim, PackedOther>() const;

protected:
	template <int i, int j, class Head, class... Args>
	void Assign(Head head, Args... args);

	template <int, int>
	void Assign() {
		// Overload to terminate recursion when assign is called with zero arguments.
	}
};


// This deduction guide is needed to complement Matrix's conversion operator to Vectors.
// While Vector does not have a constructor from a Matrix, the deduction guide seems to work
// on all compilers.
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
Vector(const Matrix<T, Rows, Columns, Order, Layout, Packed>& m) -> Vector<T, Matrix<T, Rows, Columns, Order, Layout, Packed>::vectorDim, Packed>;


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class TOther, eMatrixLayout LayoutOther, bool PackedOther>
Matrix<T, Rows, Columns, Order, Layout, Packed>::Matrix(const Matrix<TOther, Rows, Columns, Order, LayoutOther, PackedOther>& rhs) {
	if constexpr (Layout == LayoutOther) {
		for (size_t i = 0; i < stripeCount; ++i) {
			stripes[i] = Stripe(rhs.stripes[i]);
		}
	}
	else {
		for (size_t i = 0; i < RowCount(); ++i) {
			for (size_t j = 0; j < ColumnCount(); ++j) {
				(*this)(i, j) = rhs(i, j);
			}
		}
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class T2, eMatrixLayout LayoutOther, bool Packed2>
Matrix<T, Rows, Columns, Order, Layout, Packed>::Matrix(const Matrix<T2, Columns, Rows, opposite_order_v<Order>, LayoutOther, Packed2>& rhs) {
	if constexpr (Layout != LayoutOther) {
		for (size_t i = 0; i < stripeCount; ++i) {
			stripes[i] = Stripe(rhs.stripes[i]);
		}
	}
	else {
		for (size_t i = 0; i < RowCount(); ++i) {
			for (size_t j = 0; j < ColumnCount(); ++j) {
				(*this)(i, j) = rhs(j, i);
			}
		}
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class... Scalars, class>
Matrix<T, Rows, Columns, Order, Layout, Packed>::Matrix(Scalars&&... elements) {
	Assign<0, 0>(std::forward<Scalars>(elements)...);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class TOther, bool PackedOther, class>
Matrix<T, Rows, Columns, Order, Layout, Packed>::Matrix(const Vector<TOther, vectorDim, PackedOther>& v) {
	if constexpr (stripeDim == vectorDim) {
		stripes[0] = Stripe(v);
	}
	else {
		for (size_t i = 0; i < stripes.size(); ++i) {
			stripes[i][0] = v[i];
		}
	}
}

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
constexpr size_t Matrix<T, Rows, Columns, Order, Layout, Packed>::ColumnCount() const {
	return static_cast<size_t>(Columns);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
constexpr size_t Matrix<T, Rows, Columns, Order, Layout, Packed>::RowCount() const {
	return static_cast<size_t>(Rows);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
constexpr size_t Matrix<T, Rows, Columns, Order, Layout, Packed>::Width() const {
	return static_cast<size_t>(Columns);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
constexpr size_t Matrix<T, Rows, Columns, Order, Layout, Packed>::Height() const {
	return static_cast<size_t>(Rows);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T& Matrix<T, Rows, Columns, Order, Layout, Packed>::operator()(size_t row, size_t col) {
	assert(row < RowCount());
	assert(col < ColumnCount());
	return Layout == eMatrixLayout::ROW_MAJOR ? stripes[row][col] : stripes[col][row];
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
const T& Matrix<T, Rows, Columns, Order, Layout, Packed>::operator()(size_t row, size_t col) const {
	assert(row < RowCount());
	assert(col < ColumnCount());
	return Layout == eMatrixLayout::ROW_MAJOR ? stripes[row][col] : stripes[col][row];
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class TAlias>
auto Matrix<T, Rows, Columns, Order, Layout, Packed>::operator()(size_t idx)
	-> std::enable_if_t<isVector && std::is_same_v<T, TAlias>, TAlias>& {
	assert(idx < vectorDim);
	const auto [i, j] = RowCount() == 1 ? std::tuple(size_t(0), idx) : std::tuple(idx, size_t(0));
	return (*this)(i, j);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class TAlias>
auto Matrix<T, Rows, Columns, Order, Layout, Packed>::operator()(size_t idx) const
	-> const std::enable_if_t<isVector && std::is_same_v<T, TAlias>, TAlias>& {
	assert(idx < vectorDim);
	const auto [i, j] = RowCount() == 1 ? std::tuple(size_t(0), idx) : std::tuple(idx, size_t(0));
	return (*this)(i, j);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Matrix<T, Rows, Columns, Order, Layout, Packed>::Column(size_t colIdx) const {
	if constexpr (Layout == eMatrixLayout::COLUMN_MAJOR) {
		return stripes[colIdx];
	}
	else {
		Vector<T, Rows, Packed> column;
		for (size_t i = 0; i < RowCount(); ++i) {
			column(i) = (*this)(i, colIdx);
		}
		return column;
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Matrix<T, Rows, Columns, Order, Layout, Packed>::Row(size_t rowIdx) const {
	if constexpr (Layout == eMatrixLayout::ROW_MAJOR) {
		return stripes[rowIdx];
	}
	else {
		Vector<T, Columns, Packed> row;
		for (size_t i = 0; i < ColumnCount(); ++i) {
			row(i) = (*this)(rowIdx, i);
		}
		return row;
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
void Matrix<T, Rows, Columns, Order, Layout, Packed>::Column(size_t colIdx, const Vector<T, Rows, Packed>& column) {
	if constexpr (Layout == eMatrixLayout::COLUMN_MAJOR) {
		stripes[colIdx] = column;
	}
	else {
		for (size_t i = 0; i < RowCount(); ++i) {
			(*this)(i, colIdx) = column(i);
		}
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
void Matrix<T, Rows, Columns, Order, Layout, Packed>::Row(size_t rowIdx, const Vector<T, Columns, Packed>& row) {
	if constexpr (Layout == eMatrixLayout::ROW_MAJOR) {
		stripes[rowIdx] = row;
	}
	else {
		for (size_t i = 0; i < ColumnCount(); ++i) {
			(*this)(rowIdx, i) = row(i);
		}
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class TOther, bool PackedOther, class>
Matrix<T, Rows, Columns, Order, Layout, Packed>::operator Vector<TOther, Matrix::vectorDim, PackedOther>() const {
	using Vec = Vector<TOther, vectorDim, PackedOther>;
	if constexpr (Rows == 1) {
		return Vec(Row(0));
	}
	else {
		return Vec(Column(0));
	}
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <int i, int j, class Head, class... Args>
void Matrix<T, Rows, Columns, Order, Layout, Packed>::Assign(Head head, Args... args) {
	(*this)(i, j) = static_cast<T>(head);
	Assign<((j != Columns - 1) ? i : (i + 1)), ((j + 1) % Columns)>(args...);
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <int ExtractedRows, int ExtractedColumns>
Matrix<T, ExtractedRows, ExtractedColumns, Order, Layout, Packed> Matrix<T, Rows, Columns, Order, Layout, Packed>::Extract(int whereRow, int whereCol) const {
	assert(ExtractedRows + whereRow <= Rows);
	assert(ExtractedColumns + whereCol <= Columns);
	Matrix<T, ExtractedRows, ExtractedColumns, Order, Layout, Packed> submatrix;
	for (size_t row = 0; row < ExtractedRows; ++row) {
		for (size_t column = 0; column < ExtractedColumns; ++column) {
			submatrix(row, column) = (*this)(row + whereRow, column + whereCol);
		}
	}
	return submatrix;
}


template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <int InsertedRows, int InsertedColumns>
void Matrix<T, Rows, Columns, Order, Layout, Packed>::Insert(int whereRow, int whereCol, const Matrix<T, InsertedRows, InsertedColumns, Order, Layout, Packed>& submatrix) {
	assert(InsertedRows + whereRow <= Rows);
	assert(InsertedColumns + whereCol <= Columns);
	for (size_t row = 0; row < InsertedRows; ++row) {
		for (size_t column = 0; column < InsertedColumns; ++column) {
			(*this)(row + whereRow, column + whereCol) = submatrix(row, column);
		}
	}
}


} // namespace mathter