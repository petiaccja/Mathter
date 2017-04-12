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
	// Compare
	//--------------------------------------------

	template <eMatrixOrder Order2>
	bool operator==(const Matrix<T, Columns, Rows, Order2>& rhs) const {
		bool equal = true;
		for (int i = 0; i < StripeCount; ++i) {
			equal = equal && stripes[i] == rhs.stripes[i];
		}
		return equal;
	}

	template <eMatrixOrder Order2>
	bool operator!=(const Matrix<T, Columns, Rows, Order2>& rhs) const {
		return !(*this == rhs);
	}

	template <eMatrixOrder Order2, class = typename std::enable_if<std::is_floating_point<T>::value>::type>
	bool AlmostEqual(const Matrix<T, Columns, Rows, Order2>& rhs) const {
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
	auto Transposed() const -> Matrix<T, Rows, Columns, Order> {
		Matrix<T, Rows, Columns, Order> result;
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

	template <class U, int D>
	static void Scale(Vector<U, D> scale) {
		static_assert(D < std::min(Rows, Columns), "Vector dimension must be smaller than or equal to matrix dimension.");
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
	Matrix Translation(Args&& args) {
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
// Matrix-Scalar arithmetic
//------------------------------------------------------------------------------


template <class T, int Columns, int Rows, eMatrixOrder Order>
Matrix<T, Columns, Rows, Order> operator*(T s, Matrix<T, Columns, Rows, Order> mat) {
	mat *= s;
	return mat;
}

template <class T, int Columns, int Rows, eMatrixOrder Order>
Matrix<T, Columns, Rows, Order> operator/(T s, Matrix<T, Columns, Rows, Order> mat) {
	mat /= s;
	return mat;
}

template <class T, int Columns, int Rows, eMatrixOrder Order>
Matrix<T, Columns, Rows, Order> operator*(Matrix<T, Columns, Rows, Order> mat, T s) {
	mat *= s;
	return mat;
}

template <class T, int Columns, int Rows, eMatrixOrder Order>
Matrix<T, Columns, Rows, Order> operator/(Matrix<T, Columns, Rows, Order> mat, T s) {
	mat /= s;
	return mat;
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
	//static_assert(false, "Determinant not implemented yet.");
	return T();
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::Transpose() -> MyMatrixT& {
	static_cast<MyMatrixT&>(*this) = static_cast<MyMatrixT&>(*this).Transposed();
	return static_cast<MyMatrixT&>(*this);
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::Invert() -> MyMatrixT& {
	//static_assert(false, "Inverse not implemented yet.");
	return static_cast<MyMatrixT&>(*this);
}

template <class T, int Dim, eMatrixOrder Order>
auto MatrixOps<T, Dim, Dim, Order, true>::Inverted() const -> MyMatrixT {
	//static_assert(false, "Inverse not implemented yet.");
	MyMatrixT copy = static_cast<const MyMatrixT&>(*this);
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
	return static_cast<MyMatrixT&>(*this);
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


//------------------------------------------------------------------------------
// Geometry
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order>
template <class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::Rotation(T angle) {
	Matrix<T, Columns, Rows, Order> m;

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



template <class T, int Columns, int Rows, eMatrixOrder Order>
template <int Axis, class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::Rotation(T angle) {
	Matrix<T, Columns, Rows, Order> m;

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


template <class T, int Columns, int Rows, eMatrixOrder Order>
template <class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::RotationX(T angle) {
	return Rotation<0>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order>
template <class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::RotationY(T angle) {
	return Rotation<1>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order>
template <class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::RotationZ(T angle) {
	return Rotation<2>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order>
template <int Axis1, int Axis2, int Axis3, class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::Rotation(T angle1, T angle2, T angle3) {
	return Rotation<Axis1>(angle1) * Rotation<Axis2>(angle2) * Rotation<Axis3>(angle3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order>
template <class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::RotationEuler(float z1, float x2, float z3) {
	return Rotation<2, 0, 2>(z1, x2, z3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order>
template <class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::RotationRPY(float x1, float y2, float z3) {
	return Rotation<0, 1, 2>(x1, y2, z3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order>
template <class U, class>
Matrix<T, Columns, Rows, Order> Matrix<T, Columns, Rows, Order>::Rotation(Vector<U, 3> axis, T angle) {
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




template <class T, int Columns, int Rows, mathter::eMatrixOrder Order>
std::ostream& operator<<(std::ostream& os, const mathter::Matrix<T, Columns, Rows, Order>& mat) {
	for (int y = 0; y < mat.Height(); ++y) {
		os << "[";
		for (int x = 0; x < mat.Width(); ++x) {
			os << mat(x, y) << (x == mat.Width() - 1 ? "" : "\t");
		}
		os << "]\n";
	}
	return os;
}