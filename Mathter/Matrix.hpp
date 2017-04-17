#pragma once

#include <type_traits>
#include <iostream> // debug only

#include "Vector.hpp"

namespace mathter {



//------------------------------------------------------------------------------
// Matrix base class only allocating the memory
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout Layout = eMatrixLayout::ROW_MAJOR, bool Packed = false>
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
	static constexpr int StripeDim = Layout == eMatrixLayout::ROW_MAJOR ? Columns : Rows;
	static constexpr int StripeCount = Layout == eMatrixLayout::ROW_MAJOR ? Rows : Columns;

	Vector<T, StripeDim, Packed> stripes[StripeCount];

	// Get element
	inline T& GetElement(int col, int row) {
		return GetElement(col, row, std::integral_constant<bool, Layout == eMatrixLayout::ROW_MAJOR>());
	}
	inline T GetElement(int col, int row) const {
		return GetElement(col, row, std::integral_constant<bool, Layout == eMatrixLayout::ROW_MAJOR>());
	}

private:
	inline T& GetElement(int col, int row, std::true_type) {
		return stripes[row][col];
	}
	inline T GetElement(int col, int row, std::true_type) const {
		return stripes[row][col];
	}
	inline T& GetElement(int col, int row, std::false_type) {
		return stripes[col][row];
	}
	inline T GetElement(int col, int row, std::false_type) const {
		return stripes[col][row];
	}
};


//------------------------------------------------------------------------------
// Matrix operations class extending MatrixData by special operations
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed, bool Square = (Columns == Rows)>
class MatrixOps;

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixOps<T, Columns, Rows, Order, Layout, Packed, false> : public MatrixData<T, Columns, Rows, Order, Layout, Packed> {
protected:

};

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixOps<T, Dim, Dim, Order, Layout, Packed, true> : public MatrixOps<T, Dim, Dim, Order, Layout, Packed, false> {
	using MyMatrixT = Matrix<T, Dim, Dim, Order>;
public:
	template <class T2, eMatrixOrder Order2, eMatrixLayout Layout2>
	MyMatrixT& operator*=(const Matrix<T2, Dim, Dim, Order2, Layout2, Packed>& rhs) {
		auto& thisMat = static_cast<MyMatrixT&>(*this);
		thisMat = operator*<T, T2, Dim, Dim, Dim, Order, Order2, Layout, Layout2, Packed, T>(thisMat, rhs);
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

template <
	class T,
	class U,
	int Columns,
	int Rows,
	eMatrixOrder Order1,
	eMatrixOrder Order2,
	bool Packed,
	class V = decltype(T() + U())>
	Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed> operator+(
		const Matrix<T, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>&,
		const Matrix<U, Columns, Rows, Order2, eMatrixLayout::ROW_MAJOR, Packed>&);

template <
	class T,
	class U,
	int Columns,
	int Rows,
	eMatrixOrder Order1,
	eMatrixOrder Order2,
	bool Packed,
	class V = decltype(T() + U())>
	Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed> operator-(
		const Matrix<T, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>&,
		const Matrix<U, Columns, Rows, Order2, eMatrixLayout::ROW_MAJOR, Packed>&);



//------------------------------------------------------------------------------
// Matrix class providing the common interface for all matrices
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class Matrix : public MatrixOps<T, Columns, Rows, Order, Layout, Packed> {
	static_assert(Columns >= 1 && Rows >= 1, "Dimensions must be positive integers.");
protected:
	using MatrixData<T, Columns, Rows, Order, Layout, Packed>::GetElement;
	using MatrixData::stripes;
public:
	//--------------------------------------------
	// Constructors
	//--------------------------------------------

	Matrix() = default;
	
	template <class T2, eMatrixOrder Order2, eMatrixLayout Layout2>
	Matrix(const Matrix<T2, Columns, Rows, Order2, Layout2, Packed>& rhs) {
		for (int x = 0; x < Width(); ++x) {
			for (int y = 0; y < Width(); ++y) {
				(*this)(x, y) = rhs(x, y);
			}
		}
	}

	template <class T2, eMatrixOrder Order2, eMatrixLayout Layout2>
	explicit Matrix(const Matrix<T2, Columns, Rows, Order2, Layout2, !Packed>& rhs) {
		for (int x = 0; x < Width(); ++x) {
			for (int y = 0; y < Width(); ++y) {
				(*this)(x, y) = rhs(x, y);
			}
		}
	}

	template <class H, class... Args, typename std::enable_if<impl::All<impl::IsScalar, H, Args...>::value, int>::type = 0>
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
	// Compare
	//--------------------------------------------

	template <eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2>
	bool operator==(const Matrix<T, Columns, Rows, Order2, Layout2, Packed2>& rhs) const {
		bool equal = true;
		for (int i = 0; i < StripeCount; ++i) {
			equal = equal && stripes[i] == rhs.stripes[i];
		}
		return equal;
	}

	template <eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2>
	bool operator!=(const Matrix<T, Columns, Rows, Order2, Layout2, Packed2>& rhs) const {
		return !(*this == rhs);
	}

	template <eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2, class = typename std::enable_if<std::is_floating_point<T>::value>::type>
	bool AlmostEqual(const Matrix<T, Columns, Rows, Order2, Layout2, Packed2>& rhs) const {
		bool equal = true;
		for (int i = 0; i < StripeCount; ++i) {
			equal = equal && stripes[i].AlmostEqual(rhs.stripes[i]);
		}
		return equal;
	}

	//--------------------------------------------
	// Arithmetic
	//--------------------------------------------

	// Non-modifying external operators

	// Addition row-row
	template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	friend auto operator+(const Matrix<T, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
						  const Matrix<U, Columns, Rows, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
		->Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

	// Subtraction row-row
	template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	friend auto operator-(const Matrix<T, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
						  const Matrix<U, Columns, Rows, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
		->Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

	// Multiplication row-row
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
						  const Matrix<U, Columns2, Match, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
		->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

	// Multiplication row-col
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
						  const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
		->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

	// Multiplication col-col
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>& lhs,
				   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
		->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>;

	// Internal modifying operators
	template <class U, eMatrixOrder Order2, eMatrixLayout Layout2>
	Matrix& operator+=(const Matrix<U, Columns, Rows, Order2, Layout2, Packed>& rhs) {
		*this = operator+<T, U, Columns, Rows, Order, Order2, Layout2, Packed, T>(*this, rhs);
	}

	template <class U, eMatrixOrder Order2, eMatrixLayout Layout2>
	Matrix& operator-=(const Matrix<U, Columns, Rows, Order2, Layout2, Packed>& rhs) {
		*this = operator-<T, U, Columns, Rows, Order, Order2, Layout2, Packed, T>(*this, rhs);
	}

	// Scalar multiplication
	Matrix& operator*=(T s) {
		for (auto& stripe : stripes) {
			stripe *= s;
		}
		return *this;
	}
	Matrix& operator/=(T s) {
		*this *= (1 / s);
		return *this;
	}


	//--------------------------------------------
	// Matrix functions
	//--------------------------------------------
	auto Transposed() const -> Matrix<T, Rows, Columns, Order, Layout, Packed> {
		Matrix<T, Rows, Columns, Order, Layout, Packed> result;
		for (int y = 0; y < Height(); ++y) {
			for (int x = 0; x < Width(); ++x) {
				result(y, x) = (*this)(x, y);
			}
		}
		return result;
	}

	static Matrix Zero() {
		Matrix m;
		for (auto& stripe : m.stripes) {
			stripe.Spread(T(0));
		}
		return m;
	}

	Matrix& SetZero() {
		*this = Zero();
		return *this;
	}


	//--------------------------------------------
	// Geometry
	//--------------------------------------------

	// Scale
	template <class... Args, typename std::enable_if<(impl::All<impl::IsScalar, Args...>::value), int>::type = 0>
	static Matrix Scale(Args&&... args) {
		static_assert(sizeof...(Args) <= std::min(Rows, Columns), "You must provide scales for dimensions equal to matrix dimension");
		Matrix m;
		m.SetZero();
		T tableArgs[sizeof...(Args)] = { std::forward<Args>(args)... };
		for (int i = 0; i < sizeof...(Args); ++i) {
			m(i, i) = std::move(tableArgs[i]);
		}
		return m;
	}

	template <class U, int Dim>
	static void Scale(Vector<U, Dim> scale) {
		static_assert(Dim < std::min(Rows, Columns), "Vector dimension must be smaller than or equal to matrix dimension.");
		Matrix m;
		m.SetZero();
		for (int i = 0; i < scale.Dimension(); ++i) {
			m(i, i) = std::move(scale(i));
		}
		return m;
	}

	// Static rotation
	template <class = typename std::enable_if<(2 <= Rows && Rows <= 3 && 2 <= Columns && Columns <= 3)>::type>
	static Matrix Rotation(T angle);

	template <int Axis, class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix Rotation(T angle);
	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix RotationX(T angle);
	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix RotationY(T angle);
	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix RotationZ(T angle);

	template <int FirstAxis, int SecondAxis, int ThirdAxis, class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix Rotation(T angle1, T angle2, T angle3);

	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix RotationEuler(float z1, float x2, float z3);

	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix RotationRPY(float x1, float y2, float z3);

	template <class U, class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	static Matrix Rotation(Vector<U, 3> axis, T angle);


	// Member rotation
	template <class = typename std::enable_if<(2 <= Rows && Rows <= 3 && 2 <= Columns && Columns <= 3)>::type>
	Matrix& SetRotation(T angle) { *this = Matrix::Rotation(angle); return *this; }

	template <int Axis, class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotation(T angle) { *this = Matrix::Rotation<Axis>(angle); return *this; }
	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotationX(T angle) { *this = Matrix::RotationX(angle); return *this; }
	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotationY(T angle) { *this = Matrix::RotationY(angle); return *this; }
	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotationZ(T angle) { *this = Matrix::RotationZ(angle); return *this; }

	template <int FirstAxis, int SecondAxis, int ThirdAxis, class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotation(T angle1, T angle2, T angle3) { *this = Matrix::Rotation<FirstAxis, SecondAxis, ThirdAxis>(angle1, angle2, angle3); return *this; }

	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotationEuler(float z1, float x2, float z3) { *this = Matrix::RotationEuler(z1, x2, z3); return *this; }

	template <class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotationRPY(float x1, float y2, float z3) { *this = Matrix::RotationRPY(x1, y2, z3); return *this; }

	template <class U, class = typename std::enable_if<(3 <= Rows && Rows <= 4 && 3 <= Columns && Columns <= 4)>::type>
	Matrix& SetRotation(Vector<U, 3> axis, T angle) { *this = Matrix::Rotation(axis, angle); return *this; }


	// Translation
	template <class... Args, typename std::enable_if<(impl::All<impl::IsScalar, Args...>::value), int>::type = 0>
	Matrix Translation(Args&&... args) {
		constexpr int TranslationDim = Order == eMatrixOrder::FOLLOW_VECTOR ? Rows - 1 : Columns - 1;
		static_assert(sizeof...(Args) == TranslationDim, "Number of arguments must match the dimension of translation.");
		static_assert(TranslationDim > 0, "1x1 matrices cannot represent translation.");

		Matrix m;
		if (Order == eMatrixOrder::FOLLOW_VECTOR) {
			T tableArgs[sizeof...(Args)] = { std::forward<Args>(args)... };
			for (int i = 0; i < sizeof...(Args); ++i) {
				m(i, i) = std::move(tableArgs[i]);
			}
			return m;
		}
	}


	//--------------------------------------------
	// Matrix-vector arithmetic
	//--------------------------------------------
	template <class Vt, class Mt, int Vd, int Mcol, eMatrixOrder Morder, bool Packed, class Rt>
	friend Vector<Rt, Mcol, Packed> operator*(const Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Mcol, Vd, Morder, eMatrixLayout::ROW_MAJOR, Packed>& mat);

	template <class Vt, class Mt, int Vd, int Mrow, eMatrixOrder Morder, bool Packed, class Rt>
	friend Vector<Rt, Mrow, Packed> operator*(const Matrix<Mt, Vd, Mrow, Morder, eMatrixLayout::ROW_MAJOR, Packed>& mat, const Vector<Vt, Vd, Packed>& vec);


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

template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V = MatMulElemT<T, U>>
auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
	-> Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>
{
	Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed> result;

	// custom acceleration for sluggish 2x2 multiplication
	// 'if' optimized away at compile time
	if (Rows1 == 2 && Match == 2 && Columns2 == 2) {
		Vector<V, 4, Packed> aacc = { lhs(0,0), lhs(0,0), lhs(0,1), lhs(0,1) };
		Vector<V, 4, Packed> bbdd = { lhs(1,0), lhs(1,0), lhs(1,1), lhs(1,1) };
		Vector<V, 4, Packed> efef = { rhs(0,0), rhs(1,0), rhs(0,0), rhs(1,0) };
		Vector<V, 4, Packed> ghgh = { rhs(0,1), rhs(1,1), rhs(0,1), rhs(1,0) };
		Vector<V, 4, Packed> res = aacc*efef + bbdd*ghgh;
		result(0, 0) = res(0);
		result(1, 0) = res(1);
		result(0, 1) = res(2);
		result(1, 1) = res(3);
		return result;
	}

	// general algorithm
	for (int y = 0; y < Rows1; ++y) {
		result.stripes[y] = lhs(0, y) * rhs.stripes[0];
	}
	for (int x = 1; x < Match; ++x) {
		for (int y = 0; y < Rows1; ++y) {
			result.stripes[y] += lhs(x, y) * rhs.stripes[x];
		}
	}

	return result;
}

template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V = MatMulElemT<T, U>>
auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
	-> Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>
{
	Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed> result;

	for (int x = 0; x < Columns2; ++x) {
		for (int y = 0; y < Rows1; ++y) {
			result(x, y) = Vector<T, Match, Packed>::Dot(lhs.stripes[y], rhs.stripes[x]);
		}
	}

	return result;

	//return lhs*Matrix<U, Columns2, Match, Order2, eMatrixLayout::ROW_MAJOR, Packed>(rhs);
}


template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V = MatMulElemT<T, U>>
auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>& lhs,
			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
	->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>
{
	Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed> result;

	// custom acceleration for sluggish 2x2 multiplication
	// 'if' optimized away at compile time
	if (Rows1 == 2 && Match == 2 && Columns2 == 2) {
		Vector<V, 4, Packed> aacc = { lhs(0,0), lhs(0,0), lhs(0,1), lhs(0,1) };
		Vector<V, 4, Packed> bbdd = { lhs(1,0), lhs(1,0), lhs(1,1), lhs(1,1) };
		Vector<V, 4, Packed> efef = { rhs(0,0), rhs(1,0), rhs(0,0), rhs(1,0) };
		Vector<V, 4, Packed> ghgh = { rhs(0,1), rhs(1,1), rhs(0,1), rhs(1,0) };
		Vector<V, 4, Packed> res = aacc*efef + bbdd*ghgh;
		result(0, 0) = res(0);
		result(1, 0) = res(1);
		result(0, 1) = res(2);
		result(1, 1) = res(3);
		return result;
	}

	// general algorithm
	for (int x = 0; x < Columns2; ++x) {
		result.stripes[x] = lhs.stripes[0] * rhs(x, 0);
	}
	for (int y = 1; y < Match; ++y) {
		for (int x = 0; x < Columns2; ++x) {
			result.stripes[x] += lhs.stripes[y] * rhs(x, y);
		}
	}

	return result;
}




template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
auto operator+(const Matrix<T, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
			   const Matrix<U, Columns, Rows, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
	-> Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>
{
	Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed> result;
	for (int i = 0; i < Rows; ++i) {
		result.stripes[i] = lhs.stripes[i] + rhs.stripes[i];
	}
	return result;
}

template <class T, class U, int Columns, int Rows, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
auto operator-(const Matrix<T, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
			   const Matrix<U, Columns, Rows, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
	->Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed>
{
	Matrix<U, Columns, Rows, Order1, eMatrixLayout::ROW_MAJOR, Packed> result;
	for (int i = 0; i < Rows; ++i) {
		result.stripes[i] = lhs.stripes[i] - rhs.stripes[i];
	}
	return result;
}

//------------------------------------------------------------------------------
// Matrix-Scalar arithmetic
//------------------------------------------------------------------------------


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
Matrix<T, Columns, Rows, Order, Layout, Packed> operator*(T s, Matrix<T, Columns, Rows, Order, Layout, Packed> mat) {
	mat *= s;
	return mat;
}

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
Matrix<T, Columns, Rows, Order, Layout, Packed> operator/(T s, Matrix<T, Columns, Rows, Order, Layout, Packed> mat) {
	mat /= s;
	return mat;
}

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
Matrix<T, Columns, Rows, Order, Layout, Packed> operator*(Matrix<T, Columns, Rows, Order, Layout, Packed> mat, T s) {
	mat *= s;
	return mat;
}

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
Matrix<T, Columns, Rows, Order, Layout, Packed> operator/(Matrix<T, Columns, Rows, Order, Layout, Packed> mat, T s) {
	mat /= s;
	return mat;
}


//------------------------------------------------------------------------------
// MatrixOps : Square matrices
//------------------------------------------------------------------------------

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T MatrixOps<T, Dim, Dim, Order, Layout, Packed, true>::Trace() const {
	T sum = GetElement(0, 0);
	for (int i = 1; i < Dim; ++i) {
		sum += GetElement(i, i);
	}
	return sum;
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T MatrixOps<T, Dim, Dim, Order, Layout, Packed, true>::Determinant() const {
	//static_assert(false, "Determinant not implemented yet.");
	return T();
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixOps<T, Dim, Dim, Order, Layout, Packed, true>::Transpose() -> MyMatrixT& {
	static_cast<MyMatrixT&>(*this) = static_cast<MyMatrixT&>(*this).Transposed();
	return static_cast<MyMatrixT&>(*this);
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixOps<T, Dim, Dim, Order, Layout, Packed, true>::Invert() -> MyMatrixT& {
	//static_assert(false, "Inverse not implemented yet.");
	return static_cast<MyMatrixT&>(*this);
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixOps<T, Dim, Dim, Order, Layout, Packed, true>::Inverted() const -> MyMatrixT {
	//static_assert(false, "Inverse not implemented yet.");
	MyMatrixT copy = static_cast<const MyMatrixT&>(*this);
	return copy.Invert();
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixOps<T, Dim, Dim, Order, Layout, Packed, true>::Identity() -> MyMatrixT {
	MyMatrixT res;

	for (int i = 0; i < Dim; ++i) {
		static_cast<MatrixOps&>(res).stripes[i].Spread(T(0));
		res(i, i) = T(1);
	}

	return res;
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixOps<T, Dim, Dim, Order, Layout, Packed, true>::SetIdentity() -> MyMatrixT& {
	static_cast<MyMatrixT&>(*this) = Identity();
	return static_cast<MyMatrixT&>(*this);
}



//------------------------------------------------------------------------------
// Matrix-vector arithmetic
//------------------------------------------------------------------------------

// v*M
template <class Vt, class Mt, int Vd, int Mcol, eMatrixOrder Morder, bool Packed, class Rt = MatMulElemT<Vt, Mt>>
Vector<Rt, Mcol, Packed> operator*(const Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Mcol, Vd, Morder, eMatrixLayout::ROW_MAJOR, Packed>& mat) {
	Vector<Rt, Mcol, Packed> result;
	result = vec(0) * mat.stripes[0];
	for (int i = 1; i < Vd; ++i) {
		result += vec(i) * mat.stripes[i];
	}
	return result;
}

// M*v
template <class Vt, class Mt, int Vd, int Mrow, eMatrixOrder Morder, bool Packed, class Rt = MatMulElemT<Vt, Mt>>
Vector<Rt, Mrow, Packed> operator*(const Matrix<Mt, Vd, Mrow, Morder, eMatrixLayout::ROW_MAJOR, Packed>& mat, const Vector<Vt, Vd, Packed>& vec) {
	Vector<Rt, Mrow, Packed> result;
	for (int i = 0; i < Mrow; ++i) {
		result(i) = vec.Dot(vec, mat.stripes[i]);
	}
	return result;
}

// v*=M
template <class Vt, class Mt, int Vd, eMatrixOrder Morder, eMatrixLayout Layout, bool Packed>
Vector<Vt, Vd, Packed>& operator*=(Vector<Vt, Vd, Packed>& vec, const Matrix<Mt, Vd, Vd, Morder, Layout, Packed>& mat) {
	vec = operator*<Vt, Mt, Vd, Vd, Morder, Layout, Packed, Vt>(vec, mat);
	return vec;
}


//------------------------------------------------------------------------------
// Geometry
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::Rotation(T angle) 
-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	Matrix<T, Columns, Rows, Order, Layout, Packed> m;

	T C = cos(angle);
	T S = sin(angle);

	auto elem = [&m](int i, int j) -> T& {
		return Order == eMatrixOrder::FOLLOW_VECTOR ? m(i, j) : m(j, i);
	};

	// Indices according to follow vector order
	elem(0, 0) = C;		elem(1, 0) = S;
	elem(0, 1) = -S;	elem(1, 1) = C;

	// Rest
	for (int x = 2; x < m.Width(); ++x) {
		for (int y = 2; y < m.Height(); ++y) {
			m(x, y) = (x == y);
		}
	}

	return m;
}



template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <int Axis, class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::Rotation(T angle)
	-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	Matrix<T, Columns, Rows, Order, Layout, Packed> m;

	T C = cos(angle);
	T S = sin(angle);

	auto elem = [&m](int i, int j) -> T& {
		return Order == eMatrixOrder::FOLLOW_VECTOR ? m(i, j) : m(j, i);
	};

	static_assert(0 <= Axis && Axis < 3, "You may choose X=0, Y=1 or Z=2 axes.");

	// Indices according to follow vector order
	if (Axis == 0) {
		// Rotate around X
		elem(0, 0) = 1;		elem(1, 0) = 0;		elem(2, 0) = 0;
		elem(0, 1) = 0;		elem(1, 1) = C;		elem(2, 1) = S;
		elem(0, 2) = 0;		elem(1, 2) = -S;	elem(2, 2) = C;
	}
	else if (Axis == 1) {
		// Rotate around Y
		elem(0, 0) = C;		elem(1, 0) = 0;		elem(2, 0) = -S;
		elem(0, 1) = 0;		elem(1, 1) = 1;		elem(2, 1) = 0;
		elem(0, 2) = S;		elem(1, 2) = 0;		elem(2, 2) = C;
	}
	else {
		// Rotate around Z
		elem(0, 0) = C;		elem(1, 0) = S;		elem(2, 0) = 0;
		elem(0, 1) = -S;	elem(1, 1) = C;		elem(2, 1) = 0;
		elem(0, 2) = 0;		elem(1, 2) = 0;		elem(2, 2) = 1;
	}

	// Rest
	for (int x = 3; x < m.Width(); ++x) {
		for (int y = 3; y < m.Height(); ++y) {
			m(x, y) = (x == y);
		}
	}

	return m;
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::RotationX(T angle)
	-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	return Rotation<0>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::RotationY(T angle)
-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	return Rotation<1>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::RotationZ(T angle)
-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	return Rotation<2>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <int Axis1, int Axis2, int Axis3, class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::Rotation(T angle1, T angle2, T angle3)
-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	return Rotation<Axis1>(angle1) * Rotation<Axis2>(angle2) * Rotation<Axis3>(angle3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::RotationEuler(float z1, float x2, float z3)
-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	return Rotation<2, 0, 2>(z1, x2, z3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::RotationRPY(float x1, float y2, float z3)
	-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	return Rotation<0, 1, 2>(x1, y2, z3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class U, class>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::Rotation(Vector<U, 3> axis, T angle)
	-> Matrix<T, Columns, Rows, Order, Layout, Packed>
{
	Matrix<T, Columns, Rows, Order> m;

	T C = cos(angle);
	T S = sin(angle);

	// 3x3 rotation sub-matrix
	using RotMat = Matrix<T, 3, 3, eMatrixOrder::FOLLOW_VECTOR>;
	Matrix<T, 1, 3, eMatrixOrder::FOLLOW_VECTOR> u(axis(0), axis(1), axis(2));
	RotMat cross = {
		0, -u(2), u(1),
		u(2), 0, -u(0),
		-u(1), u(0), 0
	};
	RotMat rot = C*RotMat::Identity() + S*cross + (1 - C)*(u*u.Transposed());


	// Elements
	auto elem = [&m](int i, int j) -> T& {
		return Order == eMatrixOrder::PRECEDE_VECTOR ? m(i, j) : m(j, i);
	};
	for (int x = 0; x < 3; ++x) {
		for (int y = 0; y < 3; ++y) {
			elem(x, y) = rot(x, y);
		}
	}

	// Rest
	for (int x = 3; x < m.Width(); ++x) {
		for (int y = 3; y < m.Height(); ++y) {
			m(x, y) = (x == y);
		}
	}

	return m;
}



} // namespace mathter




template <class T, int Columns, int Rows, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool Packed>
std::ostream& operator<<(std::ostream& os, const mathter::Matrix<T, Columns, Rows, Order, Layout, Packed>& mat) {
	for (int y = 0; y < mat.Height(); ++y) {
		os << "[";
		for (int x = 0; x < mat.Width(); ++x) {
			os << mat(x, y) << (x == mat.Width() - 1 ? "" : "\t");
		}
		os << "]\n";
	}
	return os;
}