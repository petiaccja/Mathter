// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/DeterministicInitializer.hpp"
#include "../Common/OptimizationUtil.hpp"
#include "../Common/TypeTraits.hpp"
#include "SIMDUtil.hpp"
#include "Swizzle.hpp"

#include <algorithm>
#include <cassert>
#include <type_traits>


namespace mathter {


template <class T, int Dim, bool Packed>
struct Elements {
	using MyStorage = Storage<T, Dim, Packed>;
	using Batch = MakeBatch<T, Dim, Packed>;
	alignas(GetStorageAlignment<T, Dim, Packed>()) MyStorage array;

	Batch Load() const {
		static_assert(!std::is_void_v<Batch>, "storage is not batched");
		return Batch::load_aligned(array.data());
	}

	void Store(const std::conditional_t<!std::is_void_v<Batch>, Batch, std::nullptr_t>& value) {
		value.store_aligned(array.data());
	}

	static_assert(std::is_standard_layout_v<MyStorage>);
};


template <class T, int Dim, bool Packed>
struct VectorStorage {
	Elements<T, Dim, Packed> elements;
};


template <class T, bool Packed>
struct VectorStorage<T, 2, Packed> {
	static constexpr int Dim = 2;

	VectorStorage() {}

	union {
		Elements<T, Dim, Packed> elements;
#include "SwizzleInc/Swizzle2.hpp.inc"
	};
};


template <class T, bool Packed>
struct VectorStorage<T, 3, Packed> {
	static constexpr int Dim = 3;

	VectorStorage() {}

	union {
		Elements<T, Dim, Packed> elements;
#include "SwizzleInc/Swizzle3.hpp.inc"
	};
};


template <class T, bool Packed>
struct VectorStorage<T, 4, Packed> {
	static constexpr int Dim = 4;

	VectorStorage() {}

	union {
		Elements<T, Dim, Packed> elements;
#include "SwizzleInc/Swizzle4.hpp.inc"
	};
};


namespace check_standard_layout {

	// This is extremely important to allow safe access to unions of Swizzles.
	using Storage = VectorStorage<float, 3, false>;
	static_assert(std::is_standard_layout_v<decltype(Storage::elements)>);
	static_assert(std::is_standard_layout_v<decltype(Storage::x)>);
	static_assert(std::is_standard_layout_v<decltype(Storage::xxxx)>);
	static_assert(std::is_same_v<decltype(decltype(Storage::elements)::array), decltype(decltype(Storage::x)::array)>);
	static_assert(std::is_same_v<decltype(decltype(Storage::elements)::array), decltype(decltype(Storage::xxxx)::array)>);
	static_assert(offsetof(decltype(Storage::elements), array) == 0);
	static_assert(offsetof(decltype(Storage::x), array) == 0);
	static_assert(offsetof(decltype(Storage::xxxx), array) == 0);

} // namespace check_standard_layout


/// <summary> Represents a vector in N-dimensional space. </summary>
/// <typeparam name="T"> The scalar type on which the vector is based.
///	You can use builtin floating point or integer types. User-defined types and std::complex
/// may also work, but are not yet officially supported. </typeparam>
/// <typeparam name="Dim"> The dimension of the vector-space. Must be a positive integer. </typeparam>
///	<typeparam name="Packed"> Set to true to tightly pack vector elements and
/// avoid padding of the vector struct. Disables SIMD optimizations. </typeparam>
template <class T, int Dim, bool Packed = false>
class Vector : public VectorStorage<T, Dim, Packed> {
	static_assert(Dim >= 1, "Dimension must be positive integer.");

public:
	using MyStorage = typename decltype(VectorStorage<T, Dim, Packed>::elements)::MyStorage;
	static constexpr auto isBatched = IsBatched<T, Dim, Packed>();
	using Batch = MakeBatch<T, Dim, Packed>;
	using VectorStorage<T, Dim, Packed>::elements;

public:
	//--------------------------------------------
	// Constructors
	//--------------------------------------------

	/// <summary> Constructs the vector. Does NOT zero-initialize elements. </summary>
	Vector();

	/// <summary> Constructs the vector by converting elements of <paramref name="other"/>. </summary>
	template <class T2, bool Packed2>
	Vector(const Vector<T2, Dim, Packed2>& other);

	/// <summary> Creates a homogeneous vector by appending a 1. </summary>
	template <class T2, bool Packed2, class = std::enable_if_t<(Dim >= 2), T2>>
	explicit Vector(const Vector<T2, Dim - 1, Packed2>& rhs) : Vector(rhs, 1) {}

	/// <summary> Truncates last coordinate of homogenous vector to and performs perspective division. </summary>
	template <class T2, bool Packed2>
	explicit Vector(const Vector<T2, Dim + 1, Packed2>& rhs);

	/// <summary> Constructs the vector from an array of elements. </summary>
	/// <remarks> The number of elements must be the same as the vector's dimension. </remarks>
	explicit Vector(const T* elements);

	/// <summary> Sets all elements to the same value. </summary>
	explicit Vector(const T& all);

	/// <summary> Sets elements from a swizzle. </summary>
	template <class TOther, int DimOther, int... IndicesOther, bool PackedOther>
	Vector(const Swizzle<TOther, DimOther, PackedOther, IndicesOther...>& swizzle);

	/// <summary> Initializes the vector by concatenating given scalar, vector or swizzle arguments. </summary>
	/// <remarks> Sum of the dimension of arguments must equal vector dimension.
	///		Types of arguments may differ from vector's underlying type, in which case cast is forced without a warning. </remarks>
	template <class... Parts, class = std::enable_if_t<(sizeof...(Parts) > 1)>>
	Vector(const Parts&... parts);

	/// <summary> Construct the vector using the SIMD batch type as content. </summary>
	explicit Vector(const std::conditional_t<!std::is_void_v<Batch>, Batch, std::nullptr_t>& batch);

	//--------------------------------------------
	// Cast operators
	//--------------------------------------------
	template <class T2, class = std::enable_if_t<std::is_convertible_v<T, T2> && Dim == 1, T>>
	operator T2() const;


	//--------------------------------------------
	// Accessors
	//--------------------------------------------

	/// <summary> Returns the number of dimensions of the vector. </summary>
	constexpr int Dimension() const;

	/// <summary> Returns the nth element of the vector. </summary>
	const T& operator[](size_t idx) const;
	/// <summary> Returns the nth element of the vector. </summary>
	T& operator[](size_t idx);

	/// <summary> Returns the nth element of the vector. </summary>
	const T& operator()(size_t idx) const;
	/// <summary> Returns the nth element of the vector. </summary>
	T& operator()(size_t idx);

	/// <summary> Returns an iterator to the first element. </summary>
	auto cbegin() const;
	/// <summary> Returns an iterator to the first element. </summary>
	auto begin() const;
	/// <summary> Returns an iterator to the first element. </summary>
	auto begin();
	/// <summary> Returns an iterator to the end of the vector (works like STL). </summary>
	auto cend() const;
	/// <summary> Returns an iterator to the end of the vector (works like STL). </summary>
	auto end() const;
	/// <summary> Returns an iterator to the end of the vector (works like STL). </summary>
	auto end();

	/// <summary> Returns a pointer to the underlying array of elements. </summary>
	const T* Data() const;
	/// <summary> Returns a pointer to the underlying array of elements. </summary>
	T* Data();

	/// <summary> Returns a pointer to the underlying array of elements. </summary>
	const T* data() const;
	/// <summary> Returns a pointer to the underlying array of elements. </summary>
	T* data();

private:
	void ZeroPadding() {
		for (size_t i = Dim; i < elements.array.size(); ++i) {
			elements.array[i] = static_cast<T>(0);
		}
	}
};


namespace impl {

	template <class T>
	struct parts_scalar_type : std::conditional_t<is_scalar_v<T>, std::enable_if<true, T>, scalar_type<T>> {};

	template <class T>
	using parts_scalar_type_t = typename parts_scalar_type<T>::type;

	template <class... Parts>
	constexpr int GetConcatDim() {
		constexpr auto getDim = [](auto* arg) constexpr {
			using Arg = std::decay_t<std::remove_pointer_t<decltype(arg)>>;
			if constexpr (is_vector_v<Arg> || is_swizzle_v<Arg>) {
				return dimension_v<Arg>;
			}
			return 1;
		};
		return (... + getDim(static_cast<Parts*>(nullptr)));
	}

	template <class... Parts>
	constexpr bool GetConcatPacking() {
		constexpr auto getPacking = [](auto* arg) constexpr {
			using Arg = std::decay_t<std::remove_pointer_t<decltype(arg)>>;
			if constexpr (is_vector_v<Arg> || is_swizzle_v<Arg>) {
				return is_packed_v<Arg>;
			}
			return false;
		};
		return (... && getPacking(static_cast<Parts*>(nullptr)));
	}

} // namespace impl


template <class T, int Dim, bool Packed, int... Indices>
Vector(const Swizzle<T, Dim, Packed, Indices...>& swizzle) -> Vector<T, sizeof...(Indices), Packed>;


template <class... Parts, std::enable_if_t<(sizeof...(Parts) > 1), int> = 0>
Vector(const Parts&... parts) -> Vector<common_arithmetic_type_t<impl::parts_scalar_type_t<Parts>...>,
										impl::GetConcatDim<Parts...>(),
										impl::GetConcatPacking<Parts...>()>;


template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed>::Vector() MATHTER_VECTOR_INITIALIZER(T) {
	ZeroPadding();
}


template <class T, int Dim, bool Packed>
template <class T2, bool Packed2>
Vector<T, Dim, Packed>::Vector(const Vector<T2, Dim, Packed2>& other) {
	for (int i = 0; i < Dim; ++i) {
		elements.array[i] = static_cast<T>(other[i]);
	}
	ZeroPadding();
}


template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed>::Vector(const std::conditional_t<!std::is_void_v<Batch>, Batch, std::nullptr_t>& batch) {
	static_assert(!std::is_void_v<Batch>, "This vector is not using SIMD batches.");
	elements.Store(batch);
}


template <class T, int Dim, bool Packed>
template <class T2, class>
Vector<T, Dim, Packed>::operator T2() const {
	return static_cast<T2>(elements.array[0]);
}


template <class T, int Dim, bool Packed>
template <class T2, bool Packed2>
Vector<T, Dim, Packed>::Vector(const Vector<T2, Dim + 1, Packed2>& rhs) : Vector(rhs.data()) {
	const auto divisor = rhs[Dim];
	if constexpr (isBatched) {
		elements.Store(elements.Load() / divisor);
	}
	else {
		std::for_each(begin(), end(), [&](auto& v) { v /= divisor; });
	}
}


template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed>::Vector(const T* elements) {
	std::copy(elements, elements + Dimension(), begin());
	ZeroPadding();
}


template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed>::Vector(const T& all) {
	if constexpr (isBatched) {
		elements.Store(Batch(all));
	}
	else {
		std::fill(begin(), end(), all);
		ZeroPadding();
	}
}


template <class T, int Dim, bool Packed>
template <class TOther, int DimOther, int... IndicesOther, bool PackedOther>
Vector<T, Dim, Packed>::Vector(const Swizzle<TOther, DimOther, PackedOther, IndicesOther...>& swizzle) {
	static_assert(std::is_convertible_v<TOther, T> && Dim == sizeof...(IndicesOther));

	const auto linSwizzle = swizzle.Linearize();
	if constexpr (std::is_convertible_v<std::decay_t<decltype(linSwizzle)>, MyStorage>) {
		elements.array = linSwizzle;
	}
	else {
		for (size_t i = 0; i < Dim; ++i) {
			elements.array[i] = static_cast<T>(linSwizzle[i]);
		}
	}
}


namespace impl {

	template <size_t Offset, class T, int Dim, bool Packed, class Part, class... Parts>
	void Assign(Elements<T, Dim, Packed>& elements, const Part& part, const Parts&... parts) {
		constexpr auto partDim = size_t(dimension_v<std::decay_t<Part>>);
		if constexpr (partDim > 1) {
			for (size_t i = 0; i < partDim; ++i) {
				elements.array[Offset + i] = static_cast<T>(part[i]);
			}
		}
		else {
			elements.array[Offset] = static_cast<T>(part);
		}
		if constexpr (sizeof...(Parts) > 0) {
			Assign<Offset + partDim>(elements, parts...);
		}
	}

} // namespace impl


template <class T, int Dim, bool Packed>
template <class... Parts, class>
Vector<T, Dim, Packed>::Vector(const Parts&... parts) {
	constexpr int totalDim = (... + dimension_v<Parts>);
	static_assert(totalDim == Dim);
	impl::Assign<0>(elements, parts...);
	ZeroPadding();
}


template <class T, int Dim, bool Packed>
constexpr int Vector<T, Dim, Packed>::Dimension() const {
	return Dim;
}


template <class T, int Dim, bool Packed>
const T& Vector<T, Dim, Packed>::operator[](size_t idx) const {
	return elements.array[idx];
}


template <class T, int Dim, bool Packed>
T& Vector<T, Dim, Packed>::operator[](size_t idx) {
	return elements.array[idx];
}


template <class T, int Dim, bool Packed>
const T& Vector<T, Dim, Packed>::operator()(size_t idx) const {
	return elements.array[idx];
}


template <class T, int Dim, bool Packed>
T& Vector<T, Dim, Packed>::operator()(size_t idx) {
	return elements.array[idx];
}


template <class T, int Dim, bool Packed>
auto Vector<T, Dim, Packed>::cbegin() const {
	return elements.array.cbegin();
}


template <class T, int Dim, bool Packed>
auto Vector<T, Dim, Packed>::begin() const {
	return elements.array.begin();
}


template <class T, int Dim, bool Packed>
auto Vector<T, Dim, Packed>::begin() {
	return elements.array.begin();
}


template <class T, int Dim, bool Packed>
auto Vector<T, Dim, Packed>::cend() const {
	return cbegin() + Dimension();
}


template <class T, int Dim, bool Packed>
auto Vector<T, Dim, Packed>::end() const {
	return begin() + Dimension();
}


template <class T, int Dim, bool Packed>
auto Vector<T, Dim, Packed>::end() {
	return begin() + Dimension();
}


template <class T, int Dim, bool Packed>
const T* Vector<T, Dim, Packed>::Data() const {
	return data();
}


template <class T, int Dim, bool Packed>
T* Vector<T, Dim, Packed>::Data() {
	return data();
}


template <class T, int Dim, bool Packed>
const T* Vector<T, Dim, Packed>::data() const {
	return elements.array.data();
}


template <class T, int Dim, bool Packed>
T* Vector<T, Dim, Packed>::data() {
	return elements.array.data();
}

} // namespace mathter