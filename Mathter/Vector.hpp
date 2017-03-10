#pragma once

#include <type_traits>
#include <iostream> // debug only


static void message(const char* str) {
	std::cout << str << std::endl;
}

template <class T, int D>
class Vector;


template <class T, int D>
class VectorSpec {
public:
	T data[D];
};

template <class T>
class VectorSpec<T, 2> {
public:
	VectorSpec() = default;

	VectorSpec(T x, T y) {
		this->x = x;
		this->y = y;
		message("Vector(a,b)");
	}

	union {
		struct {
			float x, y;
		};
		float data[2];
	};

protected:
	void Assign(T x, T y) {
		this->x = x;
		this->y = y;
		message("Assign(a,b)");
	}
};

template <class T>
class VectorSpec<T, 3> {
public:
	VectorSpec() = default;

	VectorSpec(T x, T y, T z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		message("Vector(a,b,c)");
	}

	union {
		struct {
			float x, y, z;
		};
		float data[3];
	};

	static Vector<T, 3> Cross(const Vector<T, 3>& lhs, const Vector<T, 3>& rhs);

protected:
	void Assign(T x, T y, T z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		message("Assign(a,b,c)");
	}
};

template <class T>
class VectorSpec<T, 4> {
public:
	VectorSpec() = default;

	VectorSpec(T x, T y, T z = 0, T w = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
		message("Vector(a,b,c,d)");
	}

	union {
		struct {
			float x, y, z, w;
		};
		float data[4];
	};
protected:
	void Assign(T x, T y, T z = 0, T w = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
		message("Assign(a,b,c,d)");
	}
};




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



template <class T, int D>
class Vector : public VectorSpec<T, D> {
protected:
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

protected:
	using VectorSpec<T, D>::Assign;

public:
	//--------------------------------------------
	// Data constructors
	//--------------------------------------------

	// Default ctor
	Vector() {
		message("Vector = default");
	};

	// Special ctor
	using VectorSpec::VectorSpec;

	// All element same ctor
	explicit Vector(T all) {
		message("Vector(all)");
		for (int i = 0; i < D; ++i) {
			data[i] = all;
		}
	}

	// T array ctor
	explicit Vector(T* data) {
		message("Vector(a[])");
		for (int i = 0; i < D; ++i) {
			this->data[i] = data[i];
		}
	}


	//--------------------------------------------
	// Copy, assignment, set
	//--------------------------------------------

	// Copy constructor
	Vector(const Vector& rhs) {
		message("Vector(copy)");
		Assign(rhs);
	}

	// Concetaneting constructor
	template <class U, class... V, typename std::enable_if<!All<NotVector, U, V...>::value, int>::type = 0>
	Vector(const U& rhs1, const V&... rhs2) {
		message("Vector(concat...)");
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
		message("Set(a,b,c,d...)");
		Assign(u, v...);
		return *this;
	}

	// Concatenating assign
	template <class U, class... V, typename std::enable_if<!All<NotVector, U, V...>::value, int>::type = 0>
	Vector& Set(const U& u, const V&... v) {
		message("Set(concat...)");
		Assign(u, v...);
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
	Vector& operator*=(const Vector& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] *= rhs.data[i];
		}
		return *this;
	}

	Vector& operator/=(const Vector& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] /= rhs.data[i];
		}
		return *this;
	}

	Vector& operator+=(const Vector& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] += rhs.data[i];
		}
		return *this;
	}

	Vector& operator-=(const Vector& rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] -= rhs.data[i];
		}
		return *this;
	}


	// Scalar assign arithmetic
	Vector& operator*=(T rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] *= rhs;
		}
		return *this;
	}

	Vector& operator/=(T rhs) {
		T tmp = T(1) / rhs;
		for (int i = 0; i < D; ++i) {
			data[i] *= tmp;
		}
		return *this;
	}

	Vector& operator+=(T rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] += rhs;
		}
		return *this;
	}

	Vector& operator-=(T rhs) {
		for (int i = 0; i < D; ++i) {
			data[i] -= rhs;
		}
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
		T s = lhs.data[0] * rhs.data[0];
		for (int i = 1; i < D; ++i) {
			s = lhs.data[i] * rhs.data[i] + s; // hope it goes MAD
		}
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
		message("Assign(copy)");
		for (int i = 0; i < D; ++i) {
			data[i] = rhs.data[i];
		}
	}

	// Filter for concatenating assign
	template <class... Args>
	struct AssignFilter {
		static constexpr bool value = !(
			std::is_same<TypeList<Args...>, TypeList<Vector>>::value ||
			All<NotVector, Args...>::value
			);
	};

	// Concatenating assign
	template <class Head, class... Rest, typename std::enable_if<AssignFilter<Head, Rest...>::value, int>::type = 0>
	void Assign(const Head& head, const Rest&... rest) {
		message("Assign(concat...)");

		constexpr int ArgDim = SumDimensions<Head, Rest...>::value;
		static_assert(ArgDim == D, "The sum of dimensions of arguments must match the dimension of the vector.");
		AssignHelper<ArgDim>::Assign<0, Head, Rest...>(*this, head, rest...);
	}
};




// Vector copy arithmetic
template <class T, int D>
Vector<T, D> operator*(Vector<T, D> lhs, const Vector<T, D>& rhs) {
	lhs *= rhs;
	return lhs;
}

template <class T, int D>
Vector<T, D> operator/(Vector<T, D> lhs, const Vector<T, D>& rhs) {
	lhs /= rhs;
	return lhs;
}

template <class T, int D>
Vector<T, D> operator+(Vector<T, D> lhs, const Vector<T, D>& rhs) {
	lhs += rhs;
	return lhs;
}

template <class T, int D>
Vector<T, D> operator-(Vector<T, D> lhs, const Vector<T, D>& rhs) {
	lhs -= rhs;
	return lhs;
}


// Scalar assign arithmetic
template <class T, int D>
Vector<T, D> operator*(Vector<T, D> lhs, T rhs) {
	lhs *= rhs;
	return lhs;
}

template <class T, int D>
Vector<T, D> operator/(Vector<T, D> lhs, T rhs) {
	lhs /= rhs;
	return lhs;
}

template <class T, int D>
Vector<T, D> operator+(Vector<T, D> lhs, T rhs) {
	lhs += rhs;
	return lhs;
}

template <class T, int D>
Vector<T, D> operator-(Vector<T, D> lhs, T rhs) {
	lhs -= rhs;
	return lhs;
}

template <class T, int D>
Vector<T, D> operator*(T lhs, Vector<T, D> rhs) {
	rhs *= lhs;
	return rhs;
}

template <class T, int D>
Vector<T, D> operator/(T lhs, Vector<T, D> rhs) {
	rhs /= lhs;
	return rhs;
}

template <class T, int D>
Vector<T, D> operator+(T lhs, Vector<T, D> rhs) {
	rhs += lhs;
	return rhs;
}

template <class T, int D>
Vector<T, D> operator-(T lhs, Vector<T, D> rhs) {
	rhs -= lhs;
	return rhs;
}


template <class T>
Vector<T, 3> VectorSpec<T, 3>::Cross(const Vector<T, 3>& lhs, const Vector<T, 3>& rhs) {
	return Vector<T, 3>(lhs.y * rhs.z - lhs.z * rhs.y,
		lhs.z * rhs.x - lhs.x * rhs.z,
		lhs.x * rhs.y - lhs.y * rhs.x);
}


