#pragma once

// Remove goddamn fucking bullshit crapware winapi macros.
#if _MSC_VER && defined(min)
#pragma push_macro("min")
#pragma push_macro("max")
#undef min
#undef max
#define MATHTER_MINMAX
#endif


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
		return GetElementImpl(col, row, std::integral_constant<bool, Layout == eMatrixLayout::ROW_MAJOR>());
	}
	inline T GetElement(int col, int row) const {
		return GetElementImpl(col, row, std::integral_constant<bool, Layout == eMatrixLayout::ROW_MAJOR>());
	}
private:
	inline T& GetElementImpl(int col, int row, std::true_type) {
		return stripes[row][col];
	}
	inline T GetElementImpl(int col, int row, std::true_type) const {
		return stripes[row][col];
	}
	inline T& GetElementImpl(int col, int row, std::false_type) {
		return stripes[col][row];
	}
	inline T GetElementImpl(int col, int row, std::false_type) const {
		return stripes[col][row];
	}
};


//------------------------------------------------------------------------------
// Matrix operations
//------------------------------------------------------------------------------

// Empty
template <class T>
class Empty {};

template <bool Enable, class Module>
using MatrixModule = typename std::conditional<Enable, Module, Empty<Module>>::type;


// Square matrices
template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixSquare {
	using MatrixT = Matrix<T, Columns, Rows, Order, Layout, Packed>;
protected:
	friend class MatrixT;
	using Inherit = Empty<MatrixSquare>;
};

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixSquare<T, Dim, Dim, Order, Layout, Packed> {
	using MatrixT = Matrix<T, Dim, Dim, Order, Layout, Packed>;
	MatrixT& self() { return *static_cast<MatrixT*>(this); }
	const MatrixT& self() const { return *static_cast<const MatrixT*>(this); }
public:
	template <class T2, eMatrixOrder Order2, eMatrixLayout Layout2>
	MatrixT& operator*=(const Matrix<T2, Dim, Dim, Order2, Layout2, Packed>& rhs) {
		self() = operator*<T, T2, Dim, Dim, Dim, Order, Order2, Layout, Layout2, Packed, T>(self(), rhs);
		return self();
	}

	T Trace() const;
	T Determinant() const;
	MatrixT& Transpose();
	MatrixT& Invert();
	MatrixT Inverted() const;
protected:
	friend class MatrixT;
	using Inherit = MatrixSquare;
};



// Rotation 2D functions
template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixRotation2D {
	using MatrixT = Matrix<T, Columns, Rows, Order, Layout, Packed>;
	MatrixT& self() { return *static_cast<MatrixT*>(this); }
	const MatrixT& self() const { return *static_cast<const MatrixT*>(this); }
protected:
	static constexpr bool Enable2DRotation =
		(Columns == 3 && Rows == 3)
		|| (Columns == 2 && Rows == 3 && Order == eMatrixOrder::FOLLOW_VECTOR)
		|| (Columns == 3 && Rows == 2 && Order == eMatrixOrder::PRECEDE_VECTOR)
		|| (Columns == 2 && Rows == 2);
public:
	static MatrixT Rotation(T angle);
	MatrixT& SetRotation(T angle) { 
		*this = Rotation(angle);
		return *this;
	}
protected:
	friend class MatrixT;
	using Inherit = MatrixModule<Enable2DRotation, MatrixRotation2D>;
};



// Rotation 3D functions
template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixRotation3D {
	using MatrixT = Matrix<T, Columns, Rows, Order, Layout, Packed>;
	MatrixT& self() { return *static_cast<MatrixT*>(this); }
	const MatrixT& self() const { return *static_cast<const MatrixT*>(this); }
protected:
	static constexpr bool Enable3DRotation =
		(Columns == 4 && Rows == 4)
		|| (Columns == 3 && Rows == 4 && Order == eMatrixOrder::FOLLOW_VECTOR)
		|| (Columns == 4 && Rows == 3 && Order == eMatrixOrder::PRECEDE_VECTOR)
		|| (Columns == 3 && Rows == 3);
public:
	// Static rotation
	template <int Axis>
	static MatrixT RotationAxis(T angle);

	static MatrixT RotationX(T angle);
	static MatrixT RotationY(T angle);
	static MatrixT RotationZ(T angle);

	template <int FirstAxis, int SecondAxis, int ThirdAxis>
	static MatrixT RotationAxis3(T angle1, T angle2, T angle3);

	static MatrixT RotationEuler(float z1, float x2, float z3);
	static MatrixT RotationRPY(float x1, float y2, float z3);
	template <class U>
	static MatrixT RotationAxisAngle(Vector<U, 3> axis, T angle);

	// Member rotation
	template <int Axis>
	MatrixT& SetRotationAxis(T angle) {	self() = Rotation<Axis>(angle); return self(); }

	MatrixT& SetRotationX(T angle) { self() = RotationX(angle); return self(); }
	MatrixT& SetRotationY(T angle) { self() = RotationY(angle); return self(); }
	MatrixT& SetRotationZ(T angle) { self() = RotationZ(angle); return self(); }

	template <int FirstAxis, int SecondAxis, int ThirdAxis>
	MatrixT& SetRotationAxis3(T angle1, T angle2, T angle3) { self() = Rotation<FirstAxis, SecondAxis, ThirdAxis>(angle1, angle2, angle3); return self(); }

	MatrixT& SetRotationEuler(float z1, float x2, float z3) { self() = RotationEuler(z1, x3, z3); return self(); }
	MatrixT& SetRotationRPY(float x1, float y2, float z3) { self() = RotationRPY(x1, y2, z3); return self(); }
	template <class U>
	MatrixT& SetRotationAxisAngle(Vector<U, 3> axis, T angle) { self() = Rotation(axis, angle); return self(); }
protected:
	friend class MatrixT;
	using Inherit = MatrixModule<Enable3DRotation, MatrixRotation3D>;
};



// Translation functions
template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixTranslation {
	using MatrixT = Matrix<T, Columns, Rows, Order, Layout, Packed>;
	MatrixT& self() { return *static_cast<MatrixT*>(this); }
	const MatrixT& self() const { return *static_cast<const MatrixT*>(this); }
protected:
	static constexpr bool EnableTranslation =
		(Rows < Columns ? Rows : Columns) > 0 &&
		(Order == eMatrixOrder::FOLLOW_VECTOR && Rows - 1 <= Columns && Columns <= Rows)
		|| (Order == eMatrixOrder::PRECEDE_VECTOR && Columns - 1 <= Rows && Rows <= Columns);
	static constexpr int TranslationDim = Rows == Columns ? Rows - 1 : std::min(Rows, Columns);
public:
	template <class... Args, typename std::enable_if<(impl::All<impl::IsScalar, Args...>::value), int>::type = 0>
	static MatrixT Translation(Args&&... args) {
		static_assert(sizeof...(Args) == TranslationDim, "Number of arguments must match the dimension of translation.");

		MatrixT m;
		m.SetIdentity();
		T tableArgs[sizeof...(Args)] = { (T)std::forward<Args>(args)... };
		if (Order == eMatrixOrder::FOLLOW_VECTOR) {
			for (int i = 0; i < sizeof...(Args); ++i) {
				m(i, Rows - 1) = std::move(tableArgs[i]);
			}
		}
		else {
			for (int i = 0; i < sizeof...(Args); ++i) {
				m(Columns - 1, i) = std::move(tableArgs[i]);
			}
		}
		return m;
	}

	template <class Vt, bool Vpacked>
	static MatrixT Translation(const Vector<Vt, TranslationDim, Vpacked>& translation) {
		MatrixT m;
		m.SetIdentity();
		if (Order == eMatrixOrder::FOLLOW_VECTOR) {
			for (int i = 0; i < translation.Dimension(); ++i) {
				m(i, Rows - 1) = translation(i);
			}
		}
		else {
			for (int i = 0; i < translation.Dimension(); ++i) {
				m(Columns - 1, i) = translation(i);
			}
		}
		return m;
	}

	template <class... Args>
	MatrixT& SetTranslation(Args&&... args) { self() = Translation(std::forward<Args>(args)...); return self(); }

	template <class Vt, bool Vpacked>
	MatrixT& SetTranslation(const Vector<Vt, TranslationDim, Vpacked>& translation) { self() = Translation(translation); return self(); }
protected:
	friend class MatrixT;
	using Inherit = MatrixModule<EnableTranslation, MatrixTranslation>;
};


// Scale functions
template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
class MatrixScale {
	using MatrixT = Matrix<T, Columns, Rows, Order, Layout, Packed>;
	MatrixT& self() { return *static_cast<MatrixT*>(this); }
	const MatrixT& self() const { return *static_cast<const MatrixT*>(this); }
public:
	template <class... Args, typename std::enable_if<(impl::All<impl::IsScalar, Args...>::value), int>::type = 0>
	static MatrixT Scale(Args&&... args) {
		static_assert(sizeof...(Args) <= std::min(Rows, Columns), "You must provide scales for dimensions equal to matrix dimension");
		MatrixT m;
		m.SetZero();
		T tableArgs[sizeof...(Args)] = { (T)std::forward<Args>(args)... };
		int i;
		for (i = 0; i < sizeof...(Args); ++i) {
			m(i, i) = std::move(tableArgs[i]);
		}
		for (; i < std::min(Rows, Columns); ++i) {
			m(i, i) = T(1);
		}
		return m;
	}

	template <class Vt, int Vdim, bool Vpacked>
	static void Scale(const Vector<Vt, Vdim, Vpacked>& scale) {
		static_assert(Vdim < std::min(Rows, Columns), "Vector dimension must be smaller than or equal to matrix dimension.");
		MatrixT m;
		m.SetIdentity();
		int i;
		for (i = 0; i < scale.Dimension(); ++i) {
			m(i, i) = std::move(scale(i));
		}
		for (; i < std::min(Rows, Columns); ++i) {
			m(i, i) = T(1);
		}
		return m;
	}

	template <class... Args>
	MatrixT& SetScale(Args&&... args) { self() = Scale(std::forward<Args>(args)...); return self(); }

	template <class Vt, int Vdim, bool Vpacked>
	MatrixT& SetScale(const Vector<Vt, Vdim, Vpacked>& translation) { self() = Scale(translation); return self(); }
protected:
	friend class MatrixT;
	using Inherit = MatrixScale;
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
class __declspec(empty_bases) Matrix
	: public MatrixData<T, Columns, Rows, Order, Layout, Packed>,
	public MatrixSquare<T, Columns, Rows, Order, Layout, Packed>::Inherit,
	public MatrixRotation2D<T, Columns, Rows, Order, Layout, Packed>::Inherit,
	public MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::Inherit,
	public MatrixTranslation<T, Columns, Rows, Order, Layout, Packed>::Inherit,
	public MatrixScale<T, Columns, Rows, Order, Layout, Packed>::Inherit
{
	static_assert(Columns >= 1 && Rows >= 1, "Dimensions must be positive integers.");
protected:
	using MatrixData<T, Columns, Rows, Order, Layout, Packed>::GetElement;
	using MatrixData::stripes;
public:
	static void DumpLayout(std::ostream& os) {
		Matrix* ptr = reinterpret_cast<Matrix*>(1000);
		using T1 = MatrixData<T, Columns, Rows, Order, Layout, Packed>;
		using T2 = MatrixSquare<T, Columns, Rows, Order, Layout, Packed>::Inherit;
		using T3 = MatrixRotation2D<T, Columns, Rows, Order, Layout, Packed>::Inherit;
		using T4 = MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::Inherit;
		using T5 = MatrixTranslation<T, Columns, Rows, Order, Layout, Packed>::Inherit;
		using T6 = MatrixScale<T, Columns, Rows, Order, Layout, Packed>::Inherit;
		os << "MatrixData:        " << (intptr_t)static_cast<T1*>(ptr) - 1000 << " -> " << sizeof(T1) << std::endl;
		os << "MatrixSquare:      " << (intptr_t)static_cast<T2*>(ptr) - 1000 << " -> " << sizeof(T2) << std::endl;
		os << "MatrixRotation2D:  " << (intptr_t)static_cast<T3*>(ptr) - 1000 << " -> " << sizeof(T3) << std::endl;
		os << "MatrixRotation3D:  " << (intptr_t)static_cast<T4*>(ptr) - 1000 << " -> " << sizeof(T4) << std::endl;
		os << "MatrixTranslation: " << (intptr_t)static_cast<T5*>(ptr) - 1000 << " -> " << sizeof(T5) << std::endl;
		os << "MatrixScale:       " << (intptr_t)static_cast<T6*>(ptr) - 1000 << " -> " << sizeof(T6) << std::endl;
	}

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
		static_assert(sizeof(Matrix) == sizeof(stripes), "Compiler did not optimize matrix size.");

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

	static Matrix Identity();
	Matrix& SetIdentity();


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
		result.stripes[y] = rhs.stripes[0] * lhs(0, y);
	}
	for (int x = 1; x < Match; ++x) {
		for (int y = 0; y < Rows1; ++y) {
			result.stripes[y] += rhs.stripes[x] * lhs(x, y);
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
// General matrix function
//------------------------------------------------------------------------------
template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::Identity() -> Matrix<T, Columns, Rows, Order, Layout, Packed> {
	Matrix<T, Columns, Rows, Order, Layout, Packed> res;

	res.SetZero();
	for (int i = 0; i < std::min(Rows, Columns); ++i) {
		res(i, i) = T(1);
	}

	return res;
}

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto Matrix<T, Columns, Rows, Order, Layout, Packed>::SetIdentity() ->Matrix<T, Columns, Rows, Order, Layout, Packed>& {
	*this = Identity();
	return *this;
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
// Square matrices
//------------------------------------------------------------------------------

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T MatrixSquare<T, Dim, Dim, Order, Layout, Packed>::Trace() const {
	T sum = self()(0, 0);
	for (int i = 1; i < Dim; ++i) {
		sum += self()(i, i);
	}
	return sum;
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
T MatrixSquare<T, Dim, Dim, Order, Layout, Packed>::Determinant() const {
	//static_assert(false, "Determinant not implemented yet.");
	return T();
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixSquare<T, Dim, Dim, Order, Layout, Packed>::Transpose() -> MatrixT& {
	self() = self().Transposed();
	return self();
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixSquare<T, Dim, Dim, Order, Layout, Packed>::Invert() -> MatrixT& {
	//static_assert(false, "Inverse not implemented yet.");
	return self();
}

template <class T, int Dim, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixSquare<T, Dim, Dim, Order, Layout, Packed>::Inverted() const -> MatrixT {
	//static_assert(false, "Inverse not implemented yet.");
	MatrixT copy = self();
	return copy.Invert();
}



//------------------------------------------------------------------------------
// Rotation 2D
//------------------------------------------------------------------------------


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixRotation2D<T, Columns, Rows, Order, Layout, Packed>::Rotation(T angle) -> MatrixT {
	MatrixT m;

	T C = cos(angle);
	T S = sin(angle);

	auto elem = [&m](int i, int j) -> T& {
		return Order == eMatrixOrder::FOLLOW_VECTOR ? m(i, j) : m(j, i);
	};

	// Indices according to follow vector order
	elem(0, 0) = C;		elem(1, 0) = S;
	elem(0, 1) = -S;	elem(1, 1) = C;

	// Rest
	for (int x = 0; x < m.Width(); ++x) {
		for (int y = (x < 2 ? 2 : 0); y < m.Height(); ++y) {
			m(x, y) = T(x == y);
		}
	}

	return m;
}


//------------------------------------------------------------------------------
// Rotation 3D
//------------------------------------------------------------------------------

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <int Axis>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationAxis(T angle) -> MatrixT
{
	MatrixT m;

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
		for (int y = (x < 3 ? 3 : 0); y < m.Height(); ++y) {
			m(x, y) = T(x == y);
		}
	}

	return m;
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationX(T angle)-> MatrixT
{
	return RotationAxis<0>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationY(T angle) -> MatrixT
{
	return RotationAxis<1>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationZ(T angle) -> MatrixT
{
	return RotationAxis<2>(angle);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <int Axis1, int Axis2, int Axis3>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationAxis3(T angle1, T angle2, T angle3) -> MatrixT
{
	return RotationAxis<Axis1>(angle1) * RotationAxis<Axis2>(angle2) * RotationAxis<Axis3>(angle3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationEuler(float z1, float x2, float z3) -> MatrixT
{
	return RotationAxis3<2, 0, 2>(z1, x2, z3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationRPY(float x1, float y2, float z3) -> MatrixT
{
	return RotationAxis3<0, 1, 2>(x1, y2, z3);
}


template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
template <class U>
auto MatrixRotation3D<T, Columns, Rows, Order, Layout, Packed>::RotationAxisAngle(Vector<U, 3> axis, T angle) -> MatrixT 
{
	MatrixT m;

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
		for (int y = (x < 3 ? 3 : 0); y < m.Height(); ++y) {
			m(x, y) = T(x == y);
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



// Remove goddamn fucking bullshit crapware winapi macros.
#if defined(MATHTER_MINMAX)
#pragma pop_macro("min")
#pragma pop_macro("max")
#endif
