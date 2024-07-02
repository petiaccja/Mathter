// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Definitions.hpp"

#include <cmath>
#include <complex>
#include <cstdint>
#include <tuple>
#include <type_traits>
#include <utility>


namespace mathter::traits {

// Vector properties
template <class VectorT>
class VectorTraitsHelper {};

template <class T_, int Dim_, bool Packed_>
class VectorTraitsHelper<Vector<T_, Dim_, Packed_>> {
public:
	using Type = T_;
	static constexpr int Dim = Dim_;
	static constexpr bool Packed = Packed_;
};

template <class T_, int Dim_, bool Packed_>
class VectorTraitsHelper<VectorData<T_, Dim_, Packed_>> {
public:
	using Type = T_;
	static constexpr int Dim = Dim_;
	static constexpr bool Packed = Packed_;
};

template <class VectorT>
class VectorTraits : public VectorTraitsHelper<typename std::decay<VectorT>::type> {};


// Matrix properties
template <class MatrixT>
class MatrixTraitsHelper {};

template <class T_, int Rows_, int Columns_, eMatrixOrder Order_, eMatrixLayout Layout_, bool Packed_>
class MatrixTraitsHelper<Matrix<T_, Rows_, Columns_, Order_, Layout_, Packed_>> {
public:
	using Type = T_;
	static constexpr int Rows = Rows_;
	static constexpr int Columns = Columns_;
	static constexpr eMatrixOrder Order = Order_;
	static constexpr eMatrixLayout Layout = Layout_;
	static constexpr bool Packed = Packed_;
};

template <class MatrixT>
class MatrixTraits : public MatrixTraitsHelper<typename std::decay<MatrixT>::type> {};


template <eMatrixOrder Order>
class OppositeOrder {
public:
	static constexpr eMatrixOrder value = (Order == eMatrixOrder::FOLLOW_VECTOR ? eMatrixOrder::PRECEDE_VECTOR : eMatrixOrder::FOLLOW_VECTOR);
};

template <eMatrixLayout Layout>
class OppositeLayout {
public:
	static constexpr eMatrixLayout value = (Layout == eMatrixLayout::ROW_MAJOR ? eMatrixLayout::COLUMN_MAJOR : eMatrixLayout::ROW_MAJOR);
};


// Common utility
template <class T, class U>
using MatMulElemT = decltype(T() * U() + T() * U());



// Template metaprogramming utilities
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


// Decide if type is Scalar, Vector or Matrix.
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

template <class Arg>
struct IsSwizzle {
	static constexpr bool value = false;
};
template <class T, int... Indices>
struct IsSwizzle<Swizzle<T, Indices...>> {
	static constexpr bool value = true;
};
template <class Arg>
struct NotSwizzle {
	static constexpr bool value = !IsSwizzle<Arg>::value;
};

template <class Arg>
struct IsVectorOrSwizzle {
	static constexpr bool value = IsVector<Arg>::value || IsSwizzle<Arg>::value;
};

template <class T>
struct IsMatrix {
	static constexpr bool value = false;
};
template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct IsMatrix<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	static constexpr bool value = true;
};
template <class T>
struct NotMatrix {
	static constexpr bool value = !IsMatrix<T>::value;
};

template <class T>
struct IsSubmatrix {
	static constexpr bool value = false;
};
template <class M, int Rows, int Columns>
struct IsSubmatrix<SubmatrixHelper<M, Rows, Columns>> {
	static constexpr bool value = true;
};
template <class T>
struct NotSubmatrix {
	static constexpr bool value = !IsSubmatrix<T>::value;
};

template <class Arg>
struct IsQuaternion {
	static constexpr bool value = false;
};
template <class T, bool Packed>
struct IsQuaternion<Quaternion<T, Packed>> {
	static constexpr bool value = true;
};
template <class Arg>
struct NotQuaternion {
	static constexpr bool value = !IsQuaternion<Arg>::value;
};


template <class T>
struct IsScalar {
	static constexpr bool value = !IsMatrix<T>::value && !IsVector<T>::value && !IsSwizzle<T>::value && !IsQuaternion<T>::value && !IsSubmatrix<T>::value;
};



template <class T>
struct same_size_int {
	template <class U>
	static constexpr auto GetSize(std::complex<U>*) {
		return sizeof(U);
	}
	template <class U>
	static constexpr auto GetSize(U*) {
		return sizeof(U);
	}
	static constexpr auto size = GetSize(static_cast<T*>(nullptr));
	using SelectType = std::tuple<std::uint8_t, std::uint16_t, void, std::uint32_t, void, void, void, std::uint64_t>;
	using type = std::tuple_element_t<size - 1, SelectType>;
};

template <class T>
using same_size_int_t = typename same_size_int<T>::type;

} // namespace mathter::traits