#pragma once

#include <type_traits>
#include <iostream> // debug only

#include "Simd.hpp"


static void message(const char* str) {
	std::cout << str << std::endl;
}


//------------------------------------------------------------------------------
// Pre-declarations
//------------------------------------------------------------------------------

// Vector

template <class T, int D>
class Vector;


// Matrix

enum class eMatrixOrder {
	PRECEDE_VECTOR,
	FOLLOW_VECTOR,
};

template <class T, class U>
using MatMulElemT = decltype(T() * U() + T() + U());

template <class T, int Columns, int Rows, eMatrixOrder Order = eMatrixOrder::FOLLOW_VECTOR>
class Matrix;

template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
->Matrix<V, Columns2, Rows1, Order1>;


//------------------------------------------------------------------------------
// Vector data containers
//--------------------------------------

// Function they must have:
// Assign(x, [y, [z, [w]]])
// mul, div, add, sum - vec x vec
// mul, div, add, sum - vec x scalar
// spread
// dot

//------------------------------------------------------------------------------
// Vector 2D data container
//------------------------------------------------------------------------------

template <class T, int D>
class VectorSpec {
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
		->Matrix<V, Columns2, Rows1, Order1>;
public:
	T data[D];

	VectorSpec() = default;
protected:
	// Assigment
	inline void Assign(std::false_type) {} // not in overload resultion, useless
	inline void spread(T all) {
		for (int i = 0; i < D; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] *= rhs.data[i];
		}
	}
	inline void div(const VectorSpec& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] /= rhs.data[i];
		};
	}
	inline void add(const VectorSpec& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] += rhs.data[i];
		}
	}
	inline void sub(const VectorSpec& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] -= rhs.data[i];
		}
	}

	// Scalar arithmetic
	inline void mul(T rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] *= rhs;
		}
	}
	inline void div(T rhs) {
		T rcprhs = T(1) / rhs;
		for (int i = 0; i < D; ++i) {
			data[i] *= rcprhs;
		}
	}
	inline void add(T rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] += rhs;
		}
	}
	inline void sub(T rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] -= rhs;
		}
	}

	// Misc
	inline float dot(const VectorSpec& rhs) const {
		float sum = this->x * rhs.x;
		sum += this->y * rhs.y;
		return sum;
	}
};


//------------------------------------------------------------------------------
// Vector 2D data container
//------------------------------------------------------------------------------

template <class T>
class VectorSpec<T, 2> {
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
		->Matrix<V, Columns2, Rows1, Order1>;
public:
	union {
		struct {
			T x, y;
		};
		T data[2];
	};

	VectorSpec() = default;
	VectorSpec(T x, T y) {
		this->x = x;
		this->y = y;
	}
protected:
	// Assigment
	inline void Assign(T x, T y) {
		this->x = x;
		this->y = y;
	}
	inline void spread(T all) {
		for (int i = 0; i < 2; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec<T, 2>& rhs) {
		this->x *= rhs.x;
		this->y *= rhs.y;
	}
	inline void div(const VectorSpec<T, 2>& rhs) {
		this->x /= rhs.x;
		this->y /= rhs.y;
	}
	inline void add(const VectorSpec<T, 2>& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
	}
	inline void sub(const VectorSpec<T, 2>& rhs) {
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
	inline float dot(const VectorSpec<T, 2>& rhs) const {
		float sum = this->x * rhs.x;
		sum += this->y * rhs.y;
		return sum;
	}
};

//------------------------------------------------------------------------------
// Vector 3D data container
//------------------------------------------------------------------------------

template <class T>
class VectorSpec<T, 3> {
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
		->Matrix<V, Columns2, Rows1, Order1>;
public:
	union {
		struct {
			T x, y, z;
		};
		T data[3];
	};

	VectorSpec() = default;
	VectorSpec(T x, T y, T z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	inline static Vector<T, 3> Cross(const Vector<T, 3>& lhs, const Vector<T, 3>& rhs);
protected:
	// Assignment
	inline void Assign(T x, T y, T z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	inline void spread(T all) {
		for (int i = 0; i < 3; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec<T, 3>& rhs) {
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
	}
	inline void div(const VectorSpec<T, 3>& rhs) {
		this->x /= rhs.x;
		this->y /= rhs.y;
		this->z /= rhs.z;
	}
	inline void add(const VectorSpec<T, 3>& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
	}
	inline void sub(const VectorSpec<T, 3>& rhs) {
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
	float dot(const VectorSpec<T, 4>& rhs) const {
		float sum = this->x * rhs.x;
		sum += this->y * rhs.y;
		sum += this->z * rhs.z;
		return sum;
	}
};


//------------------------------------------------------------------------------
// Vector 4D data container
//------------------------------------------------------------------------------

template <class T>
class VectorSpec<T, 4> {
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
		->Matrix<V, Columns2, Rows1, Order1>;
public:
	union {
		struct {
			T x, y, z, w;
		};
		T data[4];
	};

	VectorSpec() = default;
	VectorSpec(T x, T y, T z = 0, T w = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}
protected:
	// Assigment
	inline void Assign(T x, T y, T z = 0, T w = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}
	inline void spread(T all) {
		for (int i = 0; i < 4; ++i) {
			data[i] = all;
		}
	}

	// Vector arithmetic
	inline void mul(const VectorSpec<T, 4>& rhs) {
		this->x *= rhs.x;
		this->y *= rhs.y;
		this->z *= rhs.z;
		this->w *= rhs.w;
	}
	inline void div(const VectorSpec<T, 4>& rhs) {
		this->x /= rhs.x;
		this->y /= rhs.y;
		this->z /= rhs.z;
		this->w /= rhs.w;
	}
	inline void add(const VectorSpec<T, 4>& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		this->w += rhs.w;
	}
	inline void sub(const VectorSpec<T, 4>& rhs) {
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
	inline float dot(const VectorSpec<T, 4>& rhs) const {
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
class VectorSpec<float, 3> {
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
		->Matrix<V, Columns2, Rows1, Order1>;
public:
	VectorSpec() = default;

	VectorSpec(float x, float y, float z = 0) {
		simd = Simd4f::set(x, y, z, 0);
	}

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
	inline void Assign(float x, float y, float z = 0) {
		simd = Simd4f::set(x, y, z, 0);
	}
	inline void spread(float all) {
		simd = Simd4f::spread(all);
	}

	// Vector arithmetic
	inline void mul(const VectorSpec<float, 3>& rhs) {
		simd = Simd4f::mul(simd, rhs.simd);
	}
	inline void div(const VectorSpec<float, 3>& rhs) {
		simd = Simd4f::div(simd, rhs.simd);
	}
	inline void add(const VectorSpec<float, 3>& rhs) {
		simd = Simd4f::add(simd, rhs.simd);
	}
	inline void sub(const VectorSpec<float, 3>& rhs) {
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
	inline float dot(const VectorSpec<float, 3>& rhs) const {
		return Simd4f::dot(simd, rhs.simd); // both w should be 0
	}
};


//------------------------------------------------------------------------------
// Vector 4D FLOAT SIMD data container
//------------------------------------------------------------------------------

template <>
class VectorSpec<float, 4> {
	template <class T, class U, int Match, int Rows1, int Columns2, eMatrixOrder Order1, eMatrixOrder Order2, class V, class>
	friend auto operator*(const Matrix<T, Match, Rows1, Order1>& lhs, const Matrix<U, Columns2, Match, Order2>& rhs)
		->Matrix<V, Columns2, Rows1, Order1>;
public:
	VectorSpec() = default;

	VectorSpec(float x, float y, float z = 0, float w = 0) {
		Assign(x, y, z, w);
	}

	union {
		Simd4f simd;
		struct {
			float x, y, z, w;
		};
		float data[4];
	};
protected:
	// Assignment
	inline void Assign(float x, float y, float z = 0, float w = 0) {
		simd = Simd4f::set(x, y, z, w);
	}
	inline void spread(float all) {
		simd = Simd4f::spread(all);
	}

	// Vector arithmetic
	inline void mul(const VectorSpec<float, 4>& rhs) {
		simd = Simd4f::mul(simd, rhs.simd);
	}
	inline void div(const VectorSpec<float, 4>& rhs) {
		simd = Simd4f::div(simd, rhs.simd);
	}
	inline void add(const VectorSpec<float, 4>& rhs) {
		simd = Simd4f::add(simd, rhs.simd);
	}
	inline void sub(const VectorSpec<float, 4>& rhs) {
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
	inline float dot(const VectorSpec<float, 4>& rhs) const {
		return Simd4f::dot(simd, rhs.simd);
	}
};


//------------------------------------------------------------------------------
// Template magic helper classes
//------------------------------------------------------------------------------


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
	static constexpr bool value = Cond<Head>::value || All<Cond, Rest...>::value;
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


// Decide if type is Vector
template <class Arg>
struct IsVector {
	static constexpr bool value = false;
};
template <class Ta, int Da>
struct IsVector<Vector<Ta, Da>> {
	static constexpr bool value = true;
};
template <class Arg>
struct NotVector {
	static constexpr bool value = !IsVector<Arg>::value;
};


//------------------------------------------------------------------------------
// General Vector class - no specializations needed
//------------------------------------------------------------------------------

template <class T, int D>
class Vector : public VectorSpec<T, D> {
protected:
	using VectorSpec<T, D>::Assign;

public:
	//--------------------------------------------
	// Data constructors
	//--------------------------------------------

	// Default ctor
	Vector() {
		// message("Vector = default");
	};

	// Special ctor
	using VectorSpec::VectorSpec;

	// All element same ctor
	explicit Vector(T all) {
		// message("Vector(all)");
		VectorSpec<T, D>::spread(all);
	}

	// T array ctor
	explicit Vector(T* data) {
		// message("Vector(a[])");
		for (int i = 0; i < D; ++i) {
			this->data[i] = data[i];
		}
	}


	//--------------------------------------------
	// Copy, assignment, set
	//--------------------------------------------

	// Copy constructor
	Vector(const Vector& rhs) {
		// message("Vector(copy)");
		Assign(rhs);
	}

	// Contruct from all scalars
	//template <class U, class... V, typename std::enable_if<All<NotVector, U, V...>::value, int>::type = 0>
	//Vector(U&& )


	// Concetaneting constructor
	template <class U, class... V, typename std::enable_if<Any<IsVector, U, V...>::value, int>::type = 0>
	Vector(const U& rhs1, const V&... rhs2) {
		// message("Vector(concat...)");
		Assign(rhs1, rhs2...);
	}


	// Assign from any vector of same dimension
	template <class U>
	Vector& operator=(const Vector<U, D>& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] = rhs.data[i];
		}
		return *this;
	}

	// Assign from all scalars
	template <class U, class... V, typename std::enable_if<All<NotVector, U, V...>::value, int>::type = 0>
	Vector& Set(U u, V... v) {
		// message("Set(a,b,c,d...)");
		Assign(u, v...);
		return *this;
	}

	// Concatenating assign
	template <class U, class... V, typename std::enable_if<Any<IsVector, U, V...>::value, int>::type = 0>
	Vector& Set(const U& u, const V&... v) {
		// message("Set(concat...)");
		Assign(u, v...);
		return *this;
	}

	Vector& Spread(T all) {
		VectorSpec<T, D>::spread(all);
		return *this;
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
		return data + D;
	}
	const T* end() const {
		return data + D;
	}
	T* end() {
		return data + D;
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
		for (int i = 1; i < D; ++i) {
			same = same && data[i] == rhs.data[i];
		}
		return same;
	}

	bool operator!=(const Vector& rhs) const {
		return !operator==(rhs);
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
	static auto Dot(const Vector& lhs, const Vector<U, D>& rhs) {
		auto s = lhs.data[0] * rhs.data[0];
		for (int i = 1; i < D; ++i) {
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

	// Dimension of an argument (add dynamically sized vectors later)
	template <class U>
	struct DimensionOf {
		static constexpr int value = 1;
	};
	template <class Ta, int Da>
	struct DimensionOf<Vector<Ta, Da>> {
		static constexpr int value = Da;
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

	// Get nth element of an argument
	template <class U>
	struct GetVectorElement {
		static U Get(const U& u, int idx) { return u; }
	};
	template <class U, int E>
	struct GetVectorElement<Vector<U, E>> {
		static U Get(const Vector<U, E>& u, int idx) { return u.data[idx]; }
	};

	// Assign helper
	template <int ArgDim>
	struct AssignHelper {
		template <int>
		static void Assign(Vector& target) {}

		template <int StartIdx, class Head, class... Rest>
		static void Assign(Vector& target, const Head& head, const Rest&... rest) {
			constexpr int HeadDim = DimensionOf<Head>::value;
			for (int i = 0; i < HeadDim; ++i) {
				target.data[StartIdx + i] = GetVectorElement<Head>::Get(head, i);
			}
			return Assign<StartIdx + HeadDim, Rest... >(target, rest...);
		}
	};

	// Assign from vector of same type
	void Assign(const Vector& rhs) {
		// message("Assign(copy)");
		for (int i = 0; i < D; ++i) {
			data[i] = rhs.data[i];
		}
	}

	// Filter for concatenating assign
	template <class... Args>
	struct AssignFilter {
		// value is true if general assign is allowed
		static constexpr bool value = 
			!std::is_same<TypeList<typename std::decay<Args>::type...>, TypeList<Vector>>::value // Args != Vector
			&& (!All<NotVector, Args...>::value || D > 4); // not scalars only || D > 4
	};

	// Concatenating assign
	template <class Head, class... Rest, typename std::enable_if<AssignFilter<Head, Rest...>::value, int>::type = 0>
	void Assign(const Head& head, const Rest&... rest) {
		// message("Assign(concat...)");

		constexpr int ArgDim = SumDimensions<Head, Rest...>::value;
		static_assert(ArgDim == D, "The sum of dimensions of arguments must match the dimension of the vector.");
		AssignHelper<ArgDim>::Assign<0, Head, Rest...>(*this, head, rest...);
	}
};


//------------------------------------------------------------------------------
// External Vector function
//------------------------------------------------------------------------------

// Vector copy arithmetic
template <class T, int D>
inline Vector<T, D> operator*(const Vector<T, D>& lhs, const Vector<T, D>& rhs) {
	auto tmp = lhs;
	tmp *= rhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator/(const Vector<T, D>& lhs, const Vector<T, D>& rhs) {
	auto tmp = lhs;
	tmp /= rhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator+(const Vector<T, D>& lhs, const Vector<T, D>& rhs) {
	auto tmp = lhs;
	tmp += rhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator-(const Vector<T, D>& lhs, const Vector<T, D>& rhs) {
	auto tmp = lhs;
	tmp -= rhs;
	return tmp;
}


// Scalar assign arithmetic
template <class T, int D>
inline Vector<T, D> operator*(const Vector<T, D>& lhs, T rhs) {
	auto tmp = lhs;
	tmp *= rhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator/(const Vector<T, D>& lhs, T rhs) {
	auto tmp = lhs;
	tmp /= rhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator+(const Vector<T, D>& lhs, T rhs) {
	auto tmp = lhs;
	tmp += rhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator-(const Vector<T, D>& lhs, T rhs) {
	auto tmp = lhs;
	tmp -= rhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator*(T lhs, const Vector<T, D>& rhs) {
	auto tmp = rhs;
	tmp *= lhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator/(T lhs, const Vector<T, D>& rhs) {
	auto tmp = rhs;
	tmp /= lhs;
	return tmp;;
}

template <class T, int D>
inline Vector<T, D> operator+(T lhs, const Vector<T, D>& rhs) {
	auto tmp = rhs;
	tmp += lhs;
	return tmp;
}

template <class T, int D>
inline Vector<T, D> operator-(T lhs, const Vector<T, D>& rhs) {
	auto tmp = rhs;
	tmp -= lhs;
	return tmp;
}


//------------------------------------------------------------------------------
// Special functions
//------------------------------------------------------------------------------

template <class T>
Vector<T, 3> VectorSpec<T, 3>::Cross(const Vector<T, 3>& lhs, const Vector<T, 3>& rhs) {
	return Vector<T, 3>(lhs.y * rhs.z - lhs.z * rhs.y,
						lhs.z * rhs.x - lhs.x * rhs.z,
						lhs.x * rhs.y - lhs.y * rhs.x);
}


// TODO: accelerate with simd
inline Vector<float, 3> VectorSpec<float, 3>::Cross(const Vector<float, 3>& lhs, const Vector<float, 3>& rhs) {
	return Vector<float, 3>(lhs.y * rhs.z - lhs.z * rhs.y,
							lhs.z * rhs.x - lhs.x * rhs.z,
							lhs.x * rhs.y - lhs.y * rhs.x);
}

