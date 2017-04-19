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
//------------------------------------------------------------------------------

// Function they must have:
// mul, div, add, sub | vec x vec
// mul, div, add, sub | vec x scalar
// spread
// dot

// General
template <class T, int Dim, bool Packed>
class VectorData {
public:
	T data[Dim];
};


// Small vectors with x,y,z,w members
template <class T, bool Packed>
class VectorData<T, 2, Packed> {
public:
	union {
		struct { T x, y; };
		T data[2];
	};
};

template <class T, bool Packed>
class VectorData<T, 3, Packed> {
public:
	union {
		struct { T x, y, z; };
		T data[3];
	};
};

template <class T, bool Packed>
class VectorData<T, 4, Packed> {
public:
	union {
		struct { T x, y, z, w; };
		T data[4];
	};
};


// Small SIMD fp32 vectors
template <>
class VectorData<float, 2, false> {
public:
	union {
		Simd<float, 2> simd;
		struct { float x, y; };
		float data[2];
	};
};

template <>
class VectorData<float, 3, false> {
public:
	union {
		Simd<float, 4> simd;
		struct { float x, y, z; };
		float data[3];
	};
};

template <>
class VectorData<float, 4, false> {
public:
	union {
		Simd<float, 4> simd;
		struct { float x, y, z, w; };
		float data[4];
	};
};


// Small SIMD fp64 vectors
template <>
class VectorData<double, 2, false> {
public:
	union {
		Simd<double, 2> simd;
		struct { double x, y; };
		double data[2];
	};
};

template <>
class VectorData<double, 3, false> {
public:
	union {
		Simd<double, 4> simd;
		struct { double x, y, z; };
		double data[3];
	};
};

template <>
class VectorData<double, 4, false> {
public:
	union {
		Simd<double, 4> simd;
		struct { double x, y, z, w; };
		double data[4];
	};
};



//------------------------------------------------------------------------------
// Vector basic operations
//------------------------------------------------------------------------------


template <class T>
struct has_simd {
	template <class U>
	static std::false_type test(...) { return {}; }

	template <class U>
	static decltype(U::simd) test(int) { return {}; }


	static constexpr bool value = !std::is_same<std::false_type, decltype(test<T>(0))>::value;
};


template <class T, int Dim, bool Packed, bool = has_simd<VectorData<T, Dim, Packed>>::value>
class VectorOps;


// General vector ops
template <class T, int Dim, bool Packed>
class VectorOps<T, Dim, Packed, false> {
	using VectorT = Vector<T, Dim, Packed>;

	VectorT& self() { return *reinterpret_cast<VectorT*>(this); }
	const VectorT& self() const { return *reinterpret_cast<const VectorT*>(this); }
public:
	inline VectorT operator*(T rhs) const { 
		VectorT copy = self();
		for (int i = 0; i < Dim; ++i) {
			copy.data[i] *= rhs;
		}
		return copy;
	};
protected:
	// Assignment
	static inline void spread(VectorT& lhs, T all) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] = all;
		}
	}

	// Vector arithmetic
	static inline void mul(VectorT& lhs, const VectorT& rhs) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] *= rhs.data[i];
		}
	}
	static inline void div(VectorT& lhs, const VectorT& rhs) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] /= rhs.data[i];
		};
	}
	static inline void add(VectorT& lhs, const VectorT& rhs) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] += rhs.data[i];
		}
	}
	static inline void sub(VectorT& lhs, const VectorT& rhs) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] -= rhs.data[i];
		}
	}

	// Scalar arithmetic
	static inline void mul(VectorT& lhs, T rhs) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] *= rhs;
		}
	}
	static inline void div(VectorT& lhs, T rhs) {
		T rcprhs = T(1) / rhs;
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] *= rcprhs;
		}
	}
	static inline void add(VectorT& lhs, T rhs) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] += rhs;
		}
	}
	static inline void sub(VectorT& lhs, T rhs) {
		for (int i = 0; i < Dim; ++i) {
			lhs.data[i] -= rhs;
		}
	}

	// Misc
	static inline T dot(const VectorT& lhs, const VectorT& rhs) {
		T sum = 0.0f;
		for (int i = 0; i < Dim; ++i) {
			sum += lhs.data[i] * rhs.data[i];
		}
		return sum;
	}
};


// Simd accelerated vector ops
template <class T, int Dim, bool Packed>
class VectorOps<T, Dim, Packed, true> {
	using SimdT = decltype(VectorData<T, Dim, Packed>::simd);
	using VectorT = Vector<T, Dim, Packed>;

	VectorT& self() { return *reinterpret_cast<VectorT*>(this); }
	const VectorT& self() const { return *reinterpret_cast<const VectorT*>(this); }
public:
	inline VectorT operator*(T rhs) const {
		VectorT copy;
		copy.simd = SimdT::mul(self().simd, rhs);
		return copy;
	};
protected:
	// Assignment
	static inline void spread(VectorT& lhs, T all) {
		lhs.simd = SimdT::spread(all);
	}

	// Vector arithmetic
	static inline void mul(VectorT& lhs, const VectorT& rhs) {
		lhs.simd = SimdT::mul(lhs.simd, rhs.simd);
	}
	static inline void div(VectorT& lhs, const VectorT& rhs) {
		lhs.simd = SimdT::div(lhs.simd, rhs.simd);
	}
	static inline void add(VectorT& lhs, const VectorT& rhs) {
		lhs.simd = SimdT::add(lhs.simd, rhs.simd);
	}
	static inline void sub(VectorT& lhs, const VectorT& rhs) {
		lhs.simd = SimdT::sub(lhs.simd, rhs.simd);
	}

	// Scalar arithmetic
	static inline void mul(VectorT& lhs, T rhs) {
		lhs.simd = SimdT::mul(lhs.simd, rhs);
	}
	static inline void div(VectorT& lhs, T rhs) {
		lhs.simd = SimdT::div(lhs.simd, rhs);
	}
	static inline void add(VectorT& lhs, T rhs) {
		lhs.simd = SimdT::add(lhs.simd, rhs);
	}
	static inline void sub(VectorT& lhs, T rhs) {
		lhs.simd = SimdT::sub(lhs.simd, rhs);
	}

	// Misc
	static inline T dot(const VectorT& lhs, const VectorT& rhs) {
		return SimdT::dot<Dim>(lhs.simd, rhs.simd);
	}
};


//------------------------------------------------------------------------------
// Vector special operations
//------------------------------------------------------------------------------

template <class T, int Dim, bool Packed>
class VectorSpecialOps {};

template <class T, bool Packed>
class VectorSpecialOps<T, 3, Packed> {
	using VectorT = Vector<T, 3, Packed>;
public:
	static VectorT Cross(const VectorT& lhs, const VectorT& rhs) {
		return VectorT(lhs.y * rhs.z - lhs.z * rhs.y,
					   lhs.z * rhs.x - lhs.x * rhs.z,
					   lhs.x * rhs.y - lhs.y * rhs.x);
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
class __declspec(empty_bases) Vector
	: public VectorData<T, Dim, Packed>,
	public VectorOps<T, Dim, Packed>, 
	public VectorSpecialOps<T, Dim, Packed>
{
	static_assert(Dim >= 1, "Dimension must be positive integer.");
public:
	static void DumpLayout(std::ostream& os) {
		Vector* ptr = reinterpret_cast<Vector*>(1000);
		os << "VectorData:       " << (intptr_t)static_cast<VectorData<float, 4, false>*>(ptr) - 1000 << " -> " << sizeof(VectorData<float, 4, false>) << endl;
		os << "VectorOps:        " << (intptr_t)static_cast<VectorOps<float, 4, false>*>(ptr) - 1000 << " -> " << sizeof(VectorOps<float, 4, false>) << endl;
		os << "VectorSpecialOps: " << (intptr_t)static_cast<VectorSpecialOps<float, 4, false>*>(ptr) - 1000 << " -> " << sizeof(VectorSpecialOps<float, 4, false>) << endl;
		os << "Vector:           " << (intptr_t)static_cast<Vector<float, 4, false>*>(ptr) - 1000 << " -> " << sizeof(Vector<float, 4, false>) << endl;
	}

	using VectorOps::operator*;
	//--------------------------------------------
	// Data constructors
	//--------------------------------------------

	// Default ctor
	Vector() = default;


	// All element same ctor
	explicit Vector(T all) {
		VectorOps<T, Dim, Packed>::spread(*this, all);
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
		VectorOps<T, Dim, Packed>::spread(*this, all);

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
		mul(*this, rhs);
		return *this;
	}

	inline Vector& operator/=(const Vector& rhs) {
		div(*this, rhs);
		return *this;
	}

	inline Vector& operator+=(const Vector& rhs) {
		add(*this, rhs);
		return *this;
	}

	inline Vector& operator-=(const Vector& rhs) {
		sub(*this, rhs);
		return *this;
	}

	// Vector arithmetic
	inline Vector operator*(const Vector& rhs) const { return Vector(*this) *= rhs; }
	inline Vector operator/(const Vector& rhs) const { return Vector(*this) /= rhs; }
	inline Vector operator+(const Vector& rhs) const { return Vector(*this) += rhs; }
	inline Vector operator-(const Vector& rhs) const { return Vector(*this) -= rhs; }

	// Scalar assign arithmetic
	inline Vector& operator*=(T rhs) {
		mul(*this, rhs);
		return *this;
	}

	inline Vector& operator/=(T rhs) {
		div(*this, rhs);
		return *this;
	}

	inline Vector& operator+=(T rhs) {
		add(*this, rhs);
		return *this;
	}

	inline Vector& operator-=(T rhs) {
		sub(*this, rhs);
		return *this;
	}

	// Scalar arithmetic
	//inline Vector operator*(T rhs) const { return Vector(*this) *= rhs; }
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
		return dot(lhs, rhs);
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


// Vector-Scalar arithmetic
template <class T, int Dim, bool Packed, class U, class = typename std::enable_if<std::is_convertible<U, T>::value>::type>
inline Vector<T, Dim, Packed> operator*(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs*(T)lhs;
}

template <class T, int Dim, bool Packed, class U, class = typename std::enable_if<std::is_convertible<U, T>::value>::type>
inline Vector<T, Dim, Packed> operator/(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs / (T)lhs;
}

template <class T, int Dim, bool Packed, class U, class = typename std::enable_if<std::is_convertible<U, T>::value>::type>
inline Vector<T, Dim, Packed> operator+(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs + (T)lhs;
}

template <class T, int Dim, bool Packed, class U, class = typename std::enable_if<std::is_convertible<U, T>::value>::type>
inline Vector<T, Dim, Packed> operator-(U lhs, const Vector<T, Dim, Packed>& rhs) {
	return rhs - (T)lhs;
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


//------------------------------------------------------------------------------
// IO
//------------------------------------------------------------------------------

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
