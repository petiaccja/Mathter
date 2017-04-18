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
#include <iostream>

#include "Simd.hpp"


namespace mathter {

//------------------------------------------------------------------------------
// Pre-declarations
//------------------------------------------------------------------------------

// Vector

template <class T, int Dim, bool Packed = false>
class Vector;


// Matrix

enum class eMatrixOrder {
	PRECEDE_VECTOR,
	FOLLOW_VECTOR,
};

enum class eMatrixLayout {
	ROW_MAJOR,
	COLUMN_MAJOR,
};


template <class T, class U>
using MatMulElemT = decltype(T() * U() + T() + U());

template <class T, int Columns, int Rows, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout Layout = eMatrixLayout::ROW_MAJOR, bool Packed = false>
class Matrix;

template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs, 
			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
	->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
	->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>& lhs,
			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
	->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>;


//------------------------------------------------------------------------------
// Vector data containers
//--------------------------------------

// Function they must have:
// mul, div, add, sum - vec x vec
// mul, div, add, sum - vec x scalar
// spread
// dot

//------------------------------------------------------------------------------
// Vector general data container
//------------------------------------------------------------------------------

template <class T, int Dim, bool Packed>
class VectorSpec {
public:
	T data[Dim];

	VectorSpec() = default;
protected:
	// Assignment
	inline void spread(T all) {
		for (int i = 0; i < Dim; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec& rhs) {
		for (int i = 0; i < Dim; ++i) {
			data[i] *= rhs.data[i];
		}
	}
	inline void div(const VectorSpec& rhs) {
		for (int i = 0; i < Dim; ++i) {
			data[i] /= rhs.data[i];
		};
	}
	inline void add(const VectorSpec& rhs) {
		for (int i = 0; i < Dim; ++i) {
			data[i] += rhs.data[i];
		}
	}
	inline void sub(const VectorSpec& rhs) {
		for (int i = 0; i < Dim; ++i) {
			data[i] -= rhs.data[i];
		}
	}

	// Scalar arithmetic
	inline void mul(T rhs) {
		for (int i = 0; i < Dim; ++i) {
			data[i] *= rhs;
		}
	}
	inline void div(T rhs) {
		T rcprhs = T(1) / rhs;
		for (int i = 0; i < Dim; ++i) {
			data[i] *= rcprhs;
		}
	}
	inline void add(T rhs) {
		for (int i = 0; i < Dim; ++i) {
			data[i] += rhs;
		}
	}
	inline void sub(T rhs) {
		for (int i = 0; i < Dim; ++i) {
			data[i] -= rhs;
		}
	}

	// Misc
	inline float dot(const VectorSpec& rhs) const {
		float sum = 0.0f;
		for (int i = 0; i < Dim; ++i) {
			sum += data[i] * rhs.data[i];
		}
		return sum;
	}
};


//------------------------------------------------------------------------------
// Vector 2D data container
//------------------------------------------------------------------------------

template <class T, bool Packed>
class VectorSpec<T, 2, Packed> {
public:
	union {
		struct {
			T x, y;
		};
		T data[2];
	};

	VectorSpec() = default;
protected:
	// Assignment
	inline void spread(T all) {
		for (int i = 0; i < 2; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec& rhs) {
		this->x *= rhs.x;
		this->y *= rhs.y;
	}
	inline void div(const VectorSpec& rhs) {
		this->x /= rhs.x;
		this->y /= rhs.y;
	}
	inline void add(const VectorSpec& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
	}
	inline void sub(const VectorSpec& rhs) {
		this->x -= rhs.x;
		this->y -= rhs.y;
	}

	// Scalar arithmetic
	inline void mul(T rhs) {
		this->x *= rhs;
		this->y *= rhs;
	}
	inline void div(T rhs) {
		this->x /= rhs;
		this->y /= rhs;
	}
	inline void add(T rhs) {
		this->x += rhs;
		this->y += rhs;
	}
	inline void sub(T rhs) {
		this->x -= rhs;
		this->y -= rhs;
	}

	// Misc
	inline float dot(const VectorSpec& rhs) const {
		float sum = this->x * rhs.x;
		sum += this->y * rhs.y;
		return sum;
	}
};

//------------------------------------------------------------------------------
// Vector 3D data container
//------------------------------------------------------------------------------

template <class T, bool Packed>
class VectorSpec<T, 3, Packed> {
public:
	union {
		struct {
			T x, y, z;
		};
		T data[3];
	};

	VectorSpec() = default;

	inline static Vector<T, 3, Packed> Cross(const Vector<T, 3, Packed>& lhs, const Vector<T, 3, Packed>& rhs);
protected:
	// Assignment
	inline void spread(T all) {
		for (int i = 0; i < 3; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec& rhs) {
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
	}
	inline void div(const VectorSpec& rhs) {
		this->x /= rhs.x;
		this->y /= rhs.y;
		this->z /= rhs.z;
	}
	inline void add(const VectorSpec& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
	}
	inline void sub(const VectorSpec& rhs) {
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
	}

	// Scalar arithmetic
	inline void mul(T rhs) {
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
	}
	inline void div(T rhs) {
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
	}
	inline void add(T rhs) {
		this->x += rhs;
		this->y += rhs;
		this->z += rhs;
	}
	inline void sub(T rhs) {
		this->x -= rhs;
		this->y -= rhs;
		this->z -= rhs;
	}

	// Misc
	float dot(const VectorSpec& rhs) const {
		float sum = this->x * rhs.x;
		sum += this->y * rhs.y;
		sum += this->z * rhs.z;
		return sum;
	}
};


//------------------------------------------------------------------------------
// Vector 4D data container
//------------------------------------------------------------------------------

template <class T, bool Packed>
class VectorSpec<T, 4, Packed> {
public:
	union {
		struct {
			T x, y, z, w;
		};
		T data[4];
	};

	VectorSpec() = default;
protected:
	// Assigment
	inline void spread(T all) {
		for (int i = 0; i < 4; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec& rhs) {
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
		this->w *= rhs.w;
	}
	inline void div(const VectorSpec& rhs) {
		this->x /= rhs.x;
		this->y /= rhs.y;
		this->z /= rhs.z;
		this->w /= rhs.w;
	}
	inline void add(const VectorSpec& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		this->w += rhs.w;
	}
	inline void sub(const VectorSpec& rhs) {
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		this->w -= rhs.w;
	}

	// Scalar arithmetic
	inline void mul(T rhs) {
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		this->w *= rhs;
	}
	inline void div(T rhs) {
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
		this->w /= rhs;
	}
	inline void add(T rhs) {
		this->x += rhs;
		this->y += rhs;
		this->z += rhs;
		this->w += rhs;
	}
	inline void sub(T rhs) {
		this->x -= rhs;
		this->y -= rhs;
		this->z -= rhs;
		this->w -= rhs;
	}

	// Misc
	inline float dot(const VectorSpec& rhs) const {
		float sum = this->x * rhs.x;
		sum += this->y * rhs.y;
		sum += this->z * rhs.z;
		sum += this->w * rhs.w;
		return sum;
	}
};


//------------------------------------------------------------------------------
// Vector 3D FLOAT SIMD data container
//------------------------------------------------------------------------------

template <>
class VectorSpec<float, 3, false> {
public:
	VectorSpec() { simd.v[3] = 0.0f; }

	union {
		Simd4f simd;
		struct {
			float x, y, z;
		};
		float data[3];
	};

	inline static Vector<float, 3> Cross(const Vector<float, 3>& lhs, const Vector<float, 3>& rhs);
protected:
	// Assignment
	inline void spread(float all) {
		simd = Simd4f::spread(all);
	}

	// Vector arithmetic
	inline void mul(const VectorSpec<float, 3, false>& rhs) {
		simd = Simd4f::mul(simd, rhs.simd);
	}
	inline void div(const VectorSpec<float, 3, false>& rhs) {
		simd = Simd4f::div(simd, rhs.simd);
	}
	inline void add(const VectorSpec<float, 3, false>& rhs) {
		simd = Simd4f::add(simd, rhs.simd);
	}
	inline void sub(const VectorSpec<float, 3, false>& rhs) {
		simd = Simd4f::sub(simd, rhs.simd);
	}

	// Scalar arithmetic
	inline void mul(float rhs) {
		simd = Simd4f::mul(simd, rhs);
	}
	inline void div(float rhs) {
		simd = Simd4f::div(simd, rhs);
	}
	inline void add(float rhs) {
		simd = Simd4f::add(simd, rhs);
	}
	inline void sub(float rhs) {
		simd = Simd4f::sub(simd, rhs);
	}

	// Misc
	inline float dot(const VectorSpec<float, 3, false>& rhs) const {
		return Simd4f::dot(simd, rhs.simd); // both w should be 0
	}
};


//------------------------------------------------------------------------------
// Vector 4D FLOAT SIMD data container
//------------------------------------------------------------------------------

template <>
class VectorSpec<float, 4, false> {
public:
	VectorSpec() = default;

	union {
		Simd4f simd;
		struct {
			float x, y, z, w;
		};
		float data[4];
	};
protected:
	// Assignment
	inline void spread(float all) {
		simd = Simd4f::spread(all);
	}

	// Vector arithmetic
	inline void mul(const VectorSpec<float, 4, false>& rhs) {
		simd = Simd4f::mul(simd, rhs.simd);
	}
	inline void div(const VectorSpec<float, 4, false>& rhs) {
		simd = Simd4f::div(simd, rhs.simd);
	}
	inline void add(const VectorSpec<float, 4, false>& rhs) {
		simd = Simd4f::add(simd, rhs.simd);
	}
	inline void sub(const VectorSpec<float, 4, false>& rhs) {
		simd = Simd4f::sub(simd, rhs.simd);
	}

	// Scalar arithmetic
	inline void mul(float rhs) {
		simd = Simd4f::mul(simd, rhs);
	}
	inline void div(float rhs) {
		simd = Simd4f::div(simd, rhs);
	}
	inline void add(float rhs) {
		simd = Simd4f::add(simd, rhs);
	}
	inline void sub(float rhs) {
		simd = Simd4f::sub(simd, rhs);
	}

	// Misc
	inline float dot(const VectorSpec<float, 4, false>& rhs) const {
		return Simd4f::dot(simd, rhs.simd);
	}
};


//------------------------------------------------------------------------------
// Template magic helper classes
//------------------------------------------------------------------------------

namespace impl {

template <template <class> class Cond, class... T>
struct All;

template <template <class> class Cond, class Head, class... Rest>
struct All<Cond, Head, Rest...> {
	static constexpr bool value = Cond<Head>::value && All<Cond, Rest...>::value;
};

template <template <class> class Cond>
struct All<Cond> {
	static constexpr bool value = true;
};


template <template <class> class Cond, class... T>
struct Any;

template <template <class> class Cond, class Head, class... Rest>
struct Any<Cond, Head, Rest...> {
	static constexpr bool value = Cond<Head>::value || Any<Cond, Rest...>::value;
};

template <template <class> class Cond>
struct Any<Cond> {
	static constexpr bool value = false;
};



template <class... T>
struct TypeList {};

template <class Tl1, class Tl2>
struct ConcatTypeList;

template <class... T, class... U>
struct ConcatTypeList<TypeList<T...>, TypeList<U...>> {
	using type = TypeList<T..., U...>;
};

template <class T, int N>
struct RepeatType {
	using type = typename std::conditional<N <= 0, TypeList<>, typename ConcatTypeList<TypeList<T>, typename RepeatType<T, N - 1>::type>::type>::type;
};


// Decide if type is Scalar, Vector or Matrix
template <class Arg>
struct IsVector {
	static constexpr bool value = false;
};
template <class T, int Dim, bool Packed>
struct IsVector<Vector<T, Dim, Packed>> {
	static constexpr bool value = true;
};
template <class Arg>
struct NotVector {
	static constexpr bool value = !IsVector<Arg>::value;
};

template <class T>
struct IsMatrix {
	static constexpr bool value = false;
};

template <class T, int Columns, int Rows, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct IsMatrix<Matrix<T, Columns, Rows, Order, Layout, Packed>> {
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

// Dimension of an argument (add dynamically sized vectors later)
template <class U, int Along = 0>
struct DimensionOf {
	static constexpr int value = 1;
};
template <class T, int Dim, bool Packed>
struct DimensionOf<Vector<T, Dim, Packed>, 0> {
	static constexpr int value = Dim;
};

// Sum dimensions of arguments
template <class... Rest>
struct SumDimensions;

template <class Head, class... Rest>
struct SumDimensions<Head, Rest...> {
	static constexpr int value = DimensionOf<Head>::value > 0 ? DimensionOf<Head>::value + SumDimensions<Rest...>::value : -1;
};

template <>
struct SumDimensions<> {
	static constexpr int value = 0;
};


template <class T>
bool AlmostEqual(T d1, T d2) {
	if (d1 < 1e-38 && d2 < 1e-38) {
		return true;
	}
	d1 /= pow(T(10), floor(log10(abs(d1))));
	d2 /= pow(T(10), floor(log10(abs(d2))));
	d1 *= 1000.f;
	d2 *= 1000.f;
	return floor(d1) == floor(d2);
}


} // namespace impl




//------------------------------------------------------------------------------
// General Vector class - no specializations needed
//------------------------------------------------------------------------------

template <class T, int Dim, bool Packed>
class Vector : public VectorSpec<T, Dim, Packed> {
	static_assert(Dim >= 1, "Dimension must be positive integer.");

	//template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	//friend auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
	//			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::ROW_MAJOR, Packed>& rhs)
	//	->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

	//template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	//friend auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>& lhs,
	//					  const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
	//	->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::ROW_MAJOR, Packed>;

	//template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, bool Packed, class V>
	//friend auto operator*(const Matrix<T, Match, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>& lhs,
	//			   const Matrix<U, Columns2, Match, Order2, eMatrixLayout::COLUMN_MAJOR, Packed>& rhs)
	//	->Matrix<V, Columns2, Rows1, Order1, eMatrixLayout::COLUMN_MAJOR, Packed>;
public:
	//--------------------------------------------
	// Data constructors
	//--------------------------------------------

	// Default ctor
	Vector() = default;


	// All element same ctor
	explicit Vector(T all) {
		VectorSpec<T, Dim, Packed>::spread(all);
	}

	// T array ctor
	template <class U>
	explicit Vector(const U* data) {
		for (int i = 0; i < Dim; ++i) {
			this->data[i] = data[i];
		}
	}


	//--------------------------------------------
	// Copy, assignment, set
	//--------------------------------------------

	// NOTE: somehow msvc 2015 is buggy and cannot compile sizeof... checks for ctors as in Set

	// Scalar concat constructor
	template <class H1, class H2, class... Scalars, typename std::enable_if<impl::All<impl::IsScalar, H1, H2, Scalars...>::value, int>::type = 0>
	Vector(H1 h1, H2 h2, Scalars... scalars) {
		static_assert(impl::SumDimensions<H1, H2, Scalars...>::value == Dim, "Arguments must match vector dimension.");
		Assign(0, h1, h2, scalars...);
	}

	// Generalized concat constructor
	template <class H1, class... Mixed, typename std::enable_if<impl::Any<impl::IsVector, H1, Mixed...>::value, int>::type = 0>
	Vector(const H1& h1, const Mixed&... mixed) {
		static_assert(impl::SumDimensions<H1, Mixed...>::value == Dim, "Arguments must match vector dimension.");
		Assign(0, h1, mixed...);
	}

	// Scalar concat set
	template <class... Scalars, typename std::enable_if<((sizeof...(Scalars) > 1) && impl::All<impl::IsScalar, Scalars...>::value), int>::type = 0>
	Vector& Set(Scalars... scalars) {
		static_assert(impl::SumDimensions<Scalars...>::value == Dim, "Arguments must match vector dimension.");
		Assign(0, scalars...);
		return *this;
	}

	// Generalized concat set
	template <class... Mixed, typename std::enable_if<(sizeof...(Mixed) > 0) && impl::Any<impl::IsVector, Mixed...>::value, int>::type = 0>
	Vector& Set(const Mixed&... mixed) {
		static_assert(impl::SumDimensions<Mixed...>::value == Dim, "Arguments must match vector dimension.");
		Assign(0, mixed...);
		return *this;
	}

	// Set all members to certain type
	Vector& Spread(T all) {
		VectorSpec<T, Dim, Packed>::spread(all);
		return *this;
	}

	//--------------------------------------------
	// Properties
	//--------------------------------------------

	constexpr int Dimension() const {
		return Dim;
	}


	//--------------------------------------------
	// Accessors
	//--------------------------------------------

	T operator[](int idx) const {
		return data[idx];
	}

	T& operator[](int idx) {
		return data[idx];
	}

	T operator()(int idx) const {
		return data[idx];
	}

	T& operator()(int idx) {
		return data[idx];
	}

	const T* cbegin() const {
		return data;
	}
	const T* begin() const {
		return data;
	}
	T* begin() {
		return data;
	}
	const T* cend() const {
		return data + Dim;
	}
	const T* end() const {
		return data + Dim;
	}
	T* end() {
		return data + Dim;
	}


	const T* Data() const {
		return data;
	}
	T* Data() {
		return data;
	}


	//--------------------------------------------
	// Compare
	//--------------------------------------------

	bool operator==(const Vector& rhs) const {
		bool same = data[0] == rhs.data[0];
		for (int i = 1; i < Dim; ++i) {
			same = same && data[i] == rhs.data[i];
		}
		return same;
	}

	bool operator!=(const Vector& rhs) const {
		return !operator==(rhs);
	}

	template <class = typename std::enable_if<std::is_floating_point<T>::value>::type>
	bool AlmostEqual(const Vector& rhs) const {
		bool same = true;
		for (int i = 0; i < Dim; ++i) {
			T d1 = data[i], d2 = rhs.data[i];
			bool memberEqual = impl::AlmostEqual(d1, d2);
			same = same && memberEqual;
		}
		return same;
	}

	//--------------------------------------------
	// Arithmetic
	//--------------------------------------------


	// Vector assign arithmetic
	inline Vector& operator*=(const Vector& rhs) {
		mul(rhs);
		return *this;
	}

	inline Vector& operator/=(const Vector& rhs) {
		div(rhs);
		return *this;
	}

	inline Vector& operator+=(const Vector& rhs) {
		add(rhs);
		return *this;
	}

	inline Vector& operator-=(const Vector& rhs) {
		sub(rhs);
		return *this;
	}


	// Scalar assign arithmetic
	inline Vector& operator*=(T rhs) {
		mul(rhs);
		return *this;
	}

	inline Vector& operator/=(T rhs) {
		div(rhs);
		return *this;
	}

	inline Vector& operator+=(T rhs) {
		add(rhs);
		return *this;
	}

	inline Vector& operator-=(T rhs) {
		sub(rhs);
		return *this;
	}

	// Scalar arithmetic
	inline Vector operator*(T rhs) const { return Vector(*this) *= rhs; }
	inline Vector operator/(T rhs) const { return Vector(*this) /= rhs; }
	inline Vector operator+(T rhs) const { return Vector(*this) += rhs; }
	inline Vector operator-(T rhs) const { return Vector(*this) -= rhs; }


	//--------------------------------------------
	// Common functions
	//--------------------------------------------

	void Normalize() {
		T l = Length();
		operator/=(l);
	}

	Vector Normalized() const {
		Vector v = *this;
		v.Normalize();
		return v;
	}


	static T Dot(const Vector& lhs, const Vector& rhs) {
		return lhs.dot(rhs);
	}

	template <class U, typename std::enable_if<!std::is_same<T, U>::value, int>::type = 0>
	static auto Dot(const Vector& lhs, const Vector<U, Dim>& rhs) {
		auto s = lhs.data[0] * rhs.data[0];
		for (int i = 1; i < Dim; ++i) {
			s = lhs.data[i] * rhs.data[i] + s;
		}
	}


	T LengthSquared() const {
		return Dot(*this, *this);
	}

	T Length() const {
		return sqrt(LengthSquared());
	}

protected:
	//--------------------------------------------
	// Helpers
	//--------------------------------------------

	// Get nth element of an argument
	template <class U>
	struct GetVectorElement {
		static U Get(const U& u, int idx) { return u; }
	};
	template <class U, int E>
	struct GetVectorElement<Vector<U, E>> {
		static U Get(const Vector<U, E>& u, int idx) { return u.data[idx]; }
	};

	// Assign
	// Scalar concat assign
	template <class Head, class... Scalars, typename std::enable_if<impl::All<impl::IsScalar, Head, Scalars...>::value, int>::type = 0>
	void Assign(int idx, Head head, Scalars... scalars) {
		data[idx] = head;
		Assign(idx + 1, scalars...);
	}

	// Generalized concat assign
	template <class Head, class... Mixed, typename std::enable_if<impl::Any<impl::IsVector, Head, Mixed...>::value, int>::type = 0>
	void Assign(int idx, const Head& head, const Mixed&... mixed) {
		for (int i = 0; i < impl::DimensionOf<Head>::value; ++i) {
			data[idx] = GetVectorElement<Head>::Get(head, i);
			++idx;
		}
		Assign(idx, mixed...);
	}

	// Assign terminator, fill stuff with zeros
	void Assign(int idx) {
		for (; idx < Dim; idx++) {
			data[idx] = T(0);
		}
	}
};


//------------------------------------------------------------------------------
// External Vector function
//------------------------------------------------------------------------------

// Vector-Vector copy arithmetic
template <class T, int Dim, bool Packed>
inline Vector<T, Dim> operator*(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	auto tmp = lhs;
	tmp *= rhs;
	return tmp;
}

template <class T, int Dim, bool Packed>
inline Vector<T, Dim> operator/(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	auto tmp = lhs;
	tmp /= rhs;
	return tmp;
}

template <class T, int Dim, bool Packed>
inline Vector<T, Dim> operator+(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	auto tmp = lhs;
	tmp += rhs;
	return tmp;
}

template <class T, int Dim, bool Packed>
inline Vector<T, Dim> operator-(const Vector<T, Dim, Packed>& lhs, const Vector<T, Dim, Packed>& rhs) {
	auto tmp = lhs;
	tmp -= rhs;
	return tmp;
}


// Vector-Scalar arithmetic
template <class T, int Dim, bool Packed, class U>
inline Vector<T, Dim, Packed> operator*(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs*lhs;
}

template <class T, int Dim, bool Packed, class U>
inline Vector<T, Dim, Packed> operator/(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs / lhs;
}

template <class T, int Dim, bool Packed, class U>
inline Vector<T, Dim, Packed> operator+(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs + lhs;
}

template <class T, int Dim, bool Packed, class U>
inline Vector<T, Dim, Packed> operator-(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs - lhs;
}


//------------------------------------------------------------------------------
// Special functions
//------------------------------------------------------------------------------

template <class T, bool Packed>
Vector<T, 3, Packed> VectorSpec<T, 3, Packed>::Cross(const Vector<T, 3, Packed>& lhs, const Vector<T, 3, Packed>& rhs) {
	return Vector<T, 3>(lhs.y * rhs.z - lhs.z * rhs.y,
						lhs.z * rhs.x - lhs.x * rhs.z,
						lhs.x * rhs.y - lhs.y * rhs.x);
}


// Unpacked version
// TODO: accelerate with simd
inline Vector<float, 3> VectorSpec<float, 3, false>::Cross(const Vector<float, 3>& lhs, const Vector<float, 3>& rhs) {
	return Vector<float, 3>(lhs.y * rhs.z - lhs.z * rhs.y,
							lhs.z * rhs.x - lhs.x * rhs.z,
							lhs.x * rhs.y - lhs.y * rhs.x);
}


} // namespace mathter


//------------------------------------------------------------------------------
// Vector concatenation
//------------------------------------------------------------------------------
template <class T, int Dim, bool Packed, class U>
mathter::Vector<T, Dim + 1, Packed> operator|(const mathter::Vector<T, Dim, Packed>& lhs, U rhs) {
	mathter::Vector<T, Dim + 1, Packed> ret;
	for (int i = 0; i < Dim; ++i) {
		ret(i) = lhs(i);
	}
	ret(Dim) = rhs;
	return ret;
}

template <class T1, int Dim1, class T2, int Dim2, bool Packed>
mathter::Vector<T1, Dim1 + Dim2, Packed> operator|(const mathter::Vector<T1, Dim1, Packed>& lhs, const mathter::Vector<T2, Dim2, Packed>& rhs) {
	mathter::Vector<T1, Dim1 + Dim2, Packed> ret;
	for (int i = 0; i < Dim1; ++i) {
		ret(i) = lhs(i);
	}
	for (int i = 0; i < Dim2; ++i) {
		ret(Dim1 + i) = rhs(i);
	}
	return ret;
}


template <class T, int Dim, bool Packed, class U>
mathter::Vector<T, Dim + 1, Packed> operator|(U lhs, const mathter::Vector<T, Dim, Packed>& rhs) {
	mathter::Vector<T, Dim + 1, Packed> ret;
	ret(0) = lhs;
	for (int i = 0; i < Dim; ++i) {
		ret(i + 1) = rhs(i);
	}
	return ret;
}



template <class T, int Dim>
std::ostream& operator<<(std::ostream& os, const mathter::Vector<T, Dim>& v) {
	os << "[";
	for (int x = 0; x < Dim; ++x) {
		os << v(x) << (x == Dim - 1 ? "" : "\t");
	}
	os << "]";
	return os;
}


// Remove goddamn fucking bullshit crapware winapi macros.
#if defined(MATHTER_MINMAX)
#pragma pop_macro("min")
#pragma pop_macro("max")
#endif
