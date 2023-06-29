﻿// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/Definitions.hpp"
#include "../Common/DeterministicInitializer.hpp"
#include "../Common/MathUtil.hpp"
#include "../Common/Traits.hpp"
#include "../SIMD/Simd.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <type_traits>


namespace mathter {


//------------------------------------------------------------------------------
// Swizzle
//------------------------------------------------------------------------------


/// <summary> Enables element swizzling (reordering elements) for vectors. </summary>
/// <remarks>
/// To access swizzlers, use the xx, xy, xyz and similar elements of vectors.
/// Swizzlers can be used with assignments, concatenation, casting and constructors.
/// To perform arithmetic, cast swizzlers to corresponding vector type.
/// </remarks>
template <class VectorData, int... Indices>
class Swizzle {
	using T = typename traits::VectorTraits<VectorData>::Type;
	static constexpr int IndexTable[] = { Indices... };
	static constexpr int Dim = sizeof...(Indices);
	T* data() { return reinterpret_cast<T*>(this); }
	const T* data() const { return reinterpret_cast<const T*>(this); }

public:
	/// <summary> Builds the swizzled vector object. </summary>
	operator Vector<T, sizeof...(Indices), false>() const;
	/// <summary> Builds the swizzled vector object. </summary>
	operator Vector<T, sizeof...(Indices), true>() const;

	/// <summary> Sets the parent vector's elements from the right-side argument. </summary>
	/// <remarks>
	/// Example: b = {1,2,3}; a.yxz = b; -> a contains {2,1,3}.
	/// You don't have to worry about aliasing (a.xyz = a is totally fine).
	/// </remarks>
	Swizzle& operator=(const Vector<T, sizeof...(Indices), false>& rhs);
	/// <summary> Sets the parent vector's elements from the right-side argument. </summary>
	/// <remarks>
	/// Example: b = {1,2,3}; a.yxz = b; -> a contains {2,1,3}.
	/// You don't have to worry about aliasing (a.xyz = a is totally fine).
	/// </remarks>
	Swizzle& operator=(const Vector<T, sizeof...(Indices), true>& rhs);

	/// <summary> Sets the parent vector's elements from the right-side argument. </summary>
	/// <remarks>
	/// Example: b = {1,2,3}; a.yxz = b.xyz; -> a contains {2,1,3}.
	/// You don't have to worry about aliasing (a.xyz = a.zyx is totally fine).
	/// </remarks>
	template <class T2, int... Indices2, typename std::enable_if<sizeof...(Indices) == sizeof...(Indices2), int>::type = 0>
	Swizzle& operator=(const Swizzle<T2, Indices2...>& rhs) {
		*this = Vector<T, sizeof...(Indices2), false>(rhs);
		return *this;
	}

	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy[2] returns y. </summary>
	T& operator[](int idx) {
		assert(idx < Dim);
		return data()[IndexTable[idx]];
	}
	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy[2] returns y. </summary>
	T operator[](int idx) const {
		assert(idx < Dim);
		return data()[IndexTable[idx]];
	}

	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy(2) returns y. </summary>
	T& operator()(int idx) {
		assert(idx < Dim);
		return data()[IndexTable[idx]];
	}
	/// <summary> Returns the nth element of the swizzled vector. Example: v.zxy(2) returns y. </summary>
	T operator()(int idx) const {
		assert(idx < Dim);
		return data()[IndexTable[idx]];
	}

	/// <summary> Builds the swizzled vector object. </summary>
	template <bool Packed = false>
	const auto ToVector() const {
		return Vector<T, Dim, Packed>(*this);
	}

protected:
	template <int... Rest, class = typename std::enable_if<sizeof...(Rest) == 0>::type>
	void Assign(const T*) {}

	template <int Index, int... Rest>
	void Assign(const T* rhs) {
		data()[Index] = *rhs;
		return Assign<Rest...>(rhs + 1);
	}
};

template <class T, int... Indices>
constexpr int Swizzle<T, Indices...>::IndexTable[];


//------------------------------------------------------------------------------
// VectorData
//------------------------------------------------------------------------------

// General
template <class T, int Dim, bool Packed>
class VectorData {
public:
	/// <summary> Raw array containing the elements. </summary>
	T data[Dim];
};


// Small vectors with x,y,z,w members
template <class T, bool Packed>
class VectorData<T, 2, Packed> {
	using ST = T;

public:
	VectorData() {}
	VectorData(const VectorData& rhs) {
		for (int i = 0; i < 2; ++i) {
			data[i] = rhs.data[i];
		}
	}
	VectorData& operator=(const VectorData& rhs) {
		for (int i = 0; i < 2; ++i) {
			data[i] = rhs.data[i];
		}
		return *this;
	}
	union {
		struct {
			T x, y;
		};
		/// <summary> Raw array containing the elements. </summary>
		T data[2];
#include "../Swizzle/Swizzle_2.inc.hpp"
	};
};

template <class T, bool Packed>
class VectorData<T, 3, Packed> {
	using ST = T;

public:
	VectorData() {}
	VectorData(const VectorData& rhs) {
		for (int i = 0; i < 3; ++i) {
			data[i] = rhs.data[i];
		}
	}
	VectorData& operator=(const VectorData& rhs) {
		for (int i = 0; i < 3; ++i) {
			data[i] = rhs.data[i];
		}
		return *this;
	}
	union {
		struct {
			T x, y, z;
		};
		/// <summary> Raw array containing the elements. </summary>
		T data[3];
#include "../Swizzle/Swizzle_3.inc.hpp"
	};
};

template <class T, bool Packed>
class VectorData<T, 4, Packed> {
	using ST = T;

public:
	VectorData() {}
	VectorData(const VectorData& rhs) {
		for (int i = 0; i < 4; ++i) {
			data[i] = rhs.data[i];
		}
	}
	VectorData& operator=(const VectorData& rhs) {
		for (int i = 0; i < 4; ++i) {
			data[i] = rhs.data[i];
		}
		return *this;
	}
	union {
		struct {
			T x, y, z, w;
		};
		/// <summary> Raw array containing the elements. </summary>
		T data[4];
#include "../Swizzle/Swizzle_4.inc.hpp"
	};
};


// Small SIMD fp32 vectors
template <>
class VectorData<float, 2, false> {
	using ST = float;

public:
	using SimdT = Simd<ST, 2>;
	VectorData() {}
	VectorData(const VectorData& rhs) { simd = rhs.simd; }
	explicit VectorData(SimdT simd) : simd(simd) {}
	VectorData& operator=(const VectorData& rhs) {
		simd = rhs.simd;
		return *this;
	}
	union {
		/// <summary> Leave this member alone. You can't fuck it up though. </summary>
		SimdT simd;
		struct {
			ST x, y;
		};
		/// <summary> Raw array containing the elements. </summary>
		ST data[2];
#include "../Swizzle/Swizzle_2.inc.hpp"
	};
};

template <>
class VectorData<float, 3, false> {
	using ST = float;

public:
	using SimdT = Simd<ST, 4>;
	VectorData() {}
	VectorData(const VectorData& rhs) { simd = rhs.simd; }
	explicit VectorData(SimdT simd) : simd(simd) {}
	VectorData& operator=(const VectorData& rhs) {
		simd = rhs.simd;
		return *this;
	}
	union {
		/// <summary> Leave this member alone. You can't fuck it up though. </summary>
		SimdT simd;
		struct {
			ST x, y, z;
		};
		/// <summary> Raw array containing the elements. </summary>
		ST data[3];
#include "../Swizzle/Swizzle_3.inc.hpp"
	};
};

template <>
class VectorData<float, 4, false> {
	using ST = float;

public:
	using SimdT = Simd<ST, 4>;
	VectorData() {}
	VectorData(const VectorData& rhs) { simd = rhs.simd; }
	explicit VectorData(SimdT simd) : simd(simd) {}
	VectorData& operator=(const VectorData& rhs) {
		simd = rhs.simd;
		return *this;
	}
	union {
		/// <summary> Leave this member alone. You can't fuck it up though. </summary>
		SimdT simd;
		struct {
			ST x, y, z, w;
		};
		/// <summary> Raw array containing the elements. </summary>
		ST data[4];
#include "../Swizzle/Swizzle_4.inc.hpp"
	};
};

template <>
class VectorData<float, 8, false> {
public:
	using SimdT = Simd<float, 8>;
	VectorData() {}
	VectorData(const VectorData& rhs) { simd = rhs.simd; }
	explicit VectorData(SimdT simd) : simd(simd) {}
	VectorData& operator=(const VectorData& rhs) {
		simd = rhs.simd;
		return *this;
	}
	union {
		/// <summary> Leave this member alone. You can't fuck it up though. </summary>
		SimdT simd;
		/// <summary> Raw array containing the elements. </summary>
		float data[8];
	};
};


// Small SIMD fp64 vectors
template <>
class VectorData<double, 2, false> {
	using ST = double;

public:
	using SimdT = Simd<ST, 2>;
	VectorData() {}
	VectorData(const VectorData& rhs) { simd = rhs.simd; }
	explicit VectorData(SimdT simd) : simd(simd) {}
	VectorData& operator=(const VectorData& rhs) {
		simd = rhs.simd;
		return *this;
	}
	union {
		/// <summary> Leave this member alone. You can't fuck it up though. </summary>
		SimdT simd;
		struct {
			ST x, y;
		};
		/// <summary> Raw array containing the elements. </summary>
		ST data[2];
#include "../Swizzle/Swizzle_2.inc.hpp"
	};
};

template <>
class VectorData<double, 3, false> {
	using ST = double;

public:
	using SimdT = Simd<ST, 4>;
	VectorData() {}
	VectorData(const VectorData& rhs) { simd = rhs.simd; }
	explicit VectorData(SimdT simd) : simd(simd) {}
	VectorData& operator=(const VectorData& rhs) {
		simd = rhs.simd;
		return *this;
	}
	union {
		/// <summary> Leave this member alone. You can't fuck it up though. </summary>
		SimdT simd;
		struct {
			ST x, y, z;
		};
		/// <summary> Raw array containing the elements. </summary>
		ST data[3];
#include "../Swizzle/Swizzle_3.inc.hpp"
	};
};

template <>
class VectorData<double, 4, false> {
	using ST = double;

public:
	using SimdT = Simd<ST, 4>;
	VectorData() {}
	VectorData(const VectorData& rhs) { simd = rhs.simd; }
	explicit VectorData(SimdT simd) : simd(simd) {}
	VectorData& operator=(const VectorData& rhs) {
		simd = rhs.simd;
		return *this;
	}
	union {
		/// <summary> Leave this member alone. You can't fuck it up though. </summary>
		SimdT simd;
		struct {
			ST x, y, z, w;
		};
		/// <summary> Raw array containing the elements. </summary>
		ST data[4];
#include "../Swizzle/Swizzle_4.inc.hpp"
	};
};


//------------------------------------------------------------------------------
// General vector class
//------------------------------------------------------------------------------


/// <summary> Represents a vector in N-dimensional space. </summary>
/// <typeparam name="T"> The scalar type on which the vector is based.
///	You can use builtin floating point or integer types. User-defined types and std::complex
/// may also work, but are not yet officially supported. </typeparam>
/// <typeparam name="Dim"> The dimension of the vector-space. Must be a positive integer.
/// Dynamically sized vectors are not supported yet, but you'll have to use
/// <see cref="mathter::DYNAMIC"/> to define dynamically sized vectors. </typeparam>
///	<typeparam name="Packed"> Set to true to tightly pack vector elements and
/// avoid padding of the vector struct. Disables SIMD optimizations. </typeparam>
/// <remarks>
/// There is not much extraordinary to vectors, they work as you would expect.
/// - you can use common vector space airhtmetic
/// - you have common function like normalization
/// - you can multiply them with <see cref="mathter::Matrix"/> from either side
/// - you can concatenate and swizzle them.
/// </remarks>
template <class T, int Dim, bool Packed = false>
class Vector : public VectorData<T, Dim, Packed> {
	static_assert(Dim >= 1, "Dimension must be positive integer.");

public:
	using VectorData<T, Dim, Packed>::data;
	struct FromSimd_ {};
	static constexpr FromSimd_ FromSimd = {};

	//--------------------------------------------
	// Properties
	//--------------------------------------------

	/// <summary> Returns the number of dimensions of the vector. </summary>
	constexpr int Dimension() const {
		return Dim;
	}

	//--------------------------------------------
	// Data constructors
	//--------------------------------------------

	/// <summary> Constructs the vector. Does NOT zero-initialize elements. </summary>
	Vector() MATHTER_VECTOR_INITIALIZER(T) {}
	Vector(const Vector&) = default;
	Vector& operator=(const Vector&) = default;

	/// <summary> Constructs the vector by converting elements of <paramref name="other"/>. </summary>
	template <class U, bool UPacked, std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
	Vector(const Vector<U, Dim, UPacked>& other) {
		for (int i = 0; i < Dim; ++i) {
			this->data[i] = (T)other.data[i];
		}
	}

	/// <summary> Sets all elements to the same value. </summary>
	explicit Vector(T all) {
		if constexpr (!traits::HasSimd<Vector>::value) {
			for (auto& v : *this) {
				v = all;
			}
		}
		else {
			using SimdT = decltype(VectorData<T, Dim, Packed>::simd);
			this->simd = SimdT::spread(all);
		}
	}

	/// <summary> Constructs the vector from an array of elements. </summary>
	/// <remarks> The number of elements must be the same as the vector's dimension. </remarks>
	template <class U, std::enable_if_t<std::is_convertible_v<U, T>, int> = 0>
	explicit Vector(const U* data) {
		for (int i = 0; i < Dim; ++i) {
			this->data[i] = data[i];
		}
	}

	template <class SimdArgT>
	Vector(FromSimd_, SimdArgT simd) : VectorData<T, Dim, Packed>(simd) {
		static_assert(traits::HasSimd<Vector>::value, "To the developer of Mathter: don't call this unless has SIMD.");
	}

	//--------------------------------------------
	// Homogeneous up- and downcast
	//--------------------------------------------

	/// <summary> Creates a homogeneous vector by appending a 1. </summary>
	template <class T2, bool Packed2, class = typename std::enable_if<(Dim >= 2), T2>::type>
	explicit Vector(const Vector<T2, Dim - 1, Packed2>& rhs) : Vector(rhs, 1) {}

	/// <summary> Truncates last coordinate of homogenous vector to create non-homogeneous. </summary>
	template <class T2, bool Packed2>
	explicit Vector(const Vector<T2, Dim + 1, Packed2>& rhs) : Vector(rhs.data) {}


	//--------------------------------------------
	// Scalar set & concatenation
	//--------------------------------------------

	/// <summary> Initializes the vector to the given scalar elements. </summary>
	/// <remarks> Number of arguments must equal vector dimension.
	///		Types of arguments may differ from vector's underlying type, in which case cast is forced without a warning. </remarks>
	template <class... Scalars, typename std::enable_if<Dim >= 2
															&& std::conjunction_v<std::is_convertible<Scalars, T>...>
															&& sizeof...(Scalars) == Dim,
														int>::type = 0>
	Vector(Scalars... scalars) {
		if constexpr (traits::HasSimd<Vector>::value) {
			if constexpr (Dim == 3) {
				this->simd = VectorData<T, 3, Packed>::SimdT::set(T(scalars)..., T(0));
			}
			else {
				this->simd = VectorData<T, Dim, Packed>::SimdT::set(T(scalars)...);
			}
		}
		else {
			*reinterpret_cast<std::array<T, Dim>*>(this->data) = { T(scalars)... };
		}
	}

	/// <summary> Initializes the vector by concatenating given scalar, vector or swizzle arguments. </summary>
	/// <remarks> Sum of the dimension of arguments must equal vector dimension.
	///		Types of arguments may differ from vector's underlying type, in which case cast is forced without a warning. </remarks>
	template <class... Mixed, typename std::enable_if<sizeof...(Mixed) >= 1
														  && !std::conjunction_v<std::is_convertible<Mixed, T>...>
														  && traits::SumDimensions<Mixed...>::value == Dim,
													  int>::type = 0>
	Vector(const Mixed&... mixed) {
		Assign(0, mixed...);
	}


	/// <summary> Sets the vector the provided elements. </summary>
	template <class... Args>
	void Set(Args&&... args) {
		new (this) Vector(std::forward<Args>(args)...);
	}


	//--------------------------------------------
	// Accessors
	//--------------------------------------------

	/// <summary> Returns the nth element of the vector. </summary>
	T operator[](int idx) const {
		assert(idx < Dimension());
		return data[idx];
	}
	/// <summary> Returns the nth element of the vector. </summary>
	T& operator[](int idx) {
		assert(idx < Dimension());
		return data[idx];
	}

	/// <summary> Returns the nth element of the vector. </summary>
	T operator()(int idx) const {
		assert(idx < Dimension());
		return data[idx];
	}
	/// <summary> Returns the nth element of the vector. </summary>
	T& operator()(int idx) {
		assert(idx < Dimension());
		return data[idx];
	}

	/// <summary> Returns an iterator to the first element. </summary>
	const T* cbegin() const {
		return data;
	}
	/// <summary> Returns an iterator to the first element. </summary>
	const T* begin() const {
		return data;
	}
	/// <summary> Returns an iterator to the first element. </summary>
	T* begin() {
		return data;
	}
	/// <summary> Returns an iterator to the end of the vector (works like STL). </summary>
	const T* cend() const {
		return data + Dim;
	}
	/// <summary> Returns an iterator to the end of the vector (works like STL). </summary>
	const T* end() const {
		return data + Dim;
	}
	/// <summary> Returns an iterator to the end of the vector (works like STL). </summary>
	T* end() {
		return data + Dim;
	}

	/// <summary> Returns a pointer to the underlying array of elements. </summary>
	const T* Data() const {
		return data;
	}
	/// <summary> Returns a pointer to the underlying array of elements. </summary>
	T* Data() {
		return data;
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
	template <class T2, int D2, bool P2>
	struct GetVectorElement<Vector<T2, D2, P2>> {
		static T2 Get(const Vector<T2, D2, P2>& u, int idx) { return u.data[idx]; }
	};
	template <class VectorDataU, int... Indices>
	struct GetVectorElement<Swizzle<VectorDataU, Indices...>> {
		static auto Get(const Swizzle<VectorDataU, Indices...>& u, int idx) { return u[idx]; }
	};

	// Assign
	// Scalar concat assign
	template <class Head, class... Scalars, typename std::enable_if<traits::All<traits::IsScalar, Head, Scalars...>::value, int>::type = 0>
	void Assign(int idx, Head head, Scalars... scalars) {
		data[idx] = (T)head;
		Assign(idx + 1, scalars...);
	}

	// Generalized concat assign
	template <class Head, class... Mixed, typename std::enable_if<traits::Any<traits::IsVectorOrSwizzle, Head, Mixed...>::value, int>::type = 0>
	void Assign(int idx, const Head& head, const Mixed&... mixed) {
		for (int i = 0; i < traits::DimensionOf<Head>::value; ++i) {
			data[idx] = (T)GetVectorElement<Head>::Get(head, i);
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


template <class SimdT, int... Indices>
auto ShuffleReverse(SimdT arg, std::integer_sequence<int, Indices...>) {
	return SimdT::template shuffle<Indices...>(arg);
}

template <class VectorData, int... Indices>
Swizzle<VectorData, Indices...>::operator Vector<typename Swizzle<VectorData, Indices...>::T, sizeof...(Indices), false>() const {
	using DestVecT = Vector<T, sizeof...(Indices), false>;
#if !_MSVC_LANG || _MSVC_LANG <= 201703L // MSVC has some bug...
	constexpr bool SourceHasSimd = traits::HasSimd<VectorData>::value;
	constexpr bool DestHasSimd = traits::HasSimd<DestVecT>::value;
	if constexpr (SourceHasSimd && DestHasSimd) {
		using SourceSimdT = decltype(std::declval<VectorData>().simd);
		using DestSimdT = decltype(std::declval<DestVecT>().simd);
		constexpr int VectorDataDim = traits::VectorTraits<VectorData>::Dim;
		if constexpr (std::is_same_v<SourceSimdT, DestSimdT>) {
			const auto& sourceSimd = reinterpret_cast<const VectorData*>(this)->simd;
			if constexpr (sizeof...(Indices) == 3 && VectorDataDim == 3 && VectorDataDim == 4) {
				return { DestVecT::FromSimd, ShuffleReverse(sourceSimd, typename traits::ReverseIntegerSequence<std::integer_sequence<int, Indices..., 3>>::type{}) };
			}
			else if constexpr (sizeof...(Indices) == 4 && VectorDataDim == 3 && VectorDataDim == 4) {
				return { DestVecT::FromSimd, ShuffleReverse(sourceSimd, typename traits::ReverseIntegerSequence<std::integer_sequence<int, Indices...>>::type{}) };
			}
			else if constexpr (sizeof...(Indices) == 2 && VectorDataDim == 2) {
				return { DestVecT::FromSimd, ShuffleReverse(sourceSimd, typename traits::ReverseIntegerSequence<std::integer_sequence<int, Indices...>>::type{}) };
			}
		}
	}
#endif
	return DestVecT(data()[Indices]...);
}
template <class VectorData, int... Indices>
Swizzle<VectorData, Indices...>::operator Vector<typename Swizzle<VectorData, Indices...>::T, sizeof...(Indices), true>() const {
	return Vector<T, sizeof...(Indices), true>(data()[Indices]...);
}

template <class VectorData, int... Indices>
Swizzle<VectorData, Indices...>& Swizzle<VectorData, Indices...>::operator=(const Vector<T, sizeof...(Indices), false>& rhs) {
	if (data() != rhs.data) {
		Assign<Indices...>(rhs.data);
	}
	else {
		Vector<T, sizeof...(Indices), false> tmp = rhs;
		*this = tmp;
	}
	return *this;
}
template <class VectorData, int... Indices>
Swizzle<VectorData, Indices...>& Swizzle<VectorData, Indices...>::operator=(const Vector<T, sizeof...(Indices), true>& rhs) {
	if (data() != rhs.data) {
		Assign<Indices...>(rhs.data);
	}
	else {
		Vector<T, sizeof...(Indices), false> tmp = rhs;
		*this = tmp;
	}
	return *this;
}



} // namespace mathter