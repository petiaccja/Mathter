// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================


#pragma once

#include "Types.hpp"

#include <complex>
#include <cstdint>
#include <functional>
#include <type_traits>

#ifdef MATHTER_ENABLE_SIMD
#include <xsimd/xsimd.hpp>
#endif


namespace mathter {

//------------------------------------------------------------------------------
// Classification
//------------------------------------------------------------------------------

template <class>
struct is_vector : std::false_type {};

template <class T, int Dim, bool Packed>
struct is_vector<Vector<T, Dim, Packed>> : std::true_type {};

template <class T>
constexpr auto is_vector_v = is_vector<T>::value;


template <class>
struct is_swizzle : std::false_type {};

template <class T, int Dim, bool Packed, int... Indices>
struct is_swizzle<Swizzle<T, Dim, Packed, Indices...>> : std::true_type {};

template <class T>
constexpr auto is_swizzle_v = is_swizzle<T>::value;


template <class>
struct is_matrix : std::false_type {};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct is_matrix<Matrix<T, Rows, Columns, Order, Layout, Packed>> : std::true_type {};

template <class T>
constexpr auto is_matrix_v = is_matrix<T>::value;


template <class>
struct is_quaternion : std::false_type {};

template <class T, eQuaternionLayout Layout, bool Packed>
struct is_quaternion<Quaternion<T, Layout, Packed>> : std::true_type {};

template <class T>
constexpr auto is_quaternion_v = is_quaternion<T>::value;


template <class T>
struct is_scalar {
	static constexpr auto value = !(is_vector_v<T> || is_swizzle_v<T> || is_matrix_v<T> || is_quaternion_v<T>);
};

template <class T>
constexpr auto is_scalar_v = is_scalar<T>::value;


//------------------------------------------------------------------------------
// Vector type properties
//------------------------------------------------------------------------------

//--------------------------------------
// Scalar type
//--------------------------------------

template <class T>
struct scalar_type {
	using type = T;
};

template <class T, int Dim, bool Packed>
struct scalar_type<Vector<T, Dim, Packed>> {
	using type = T;
};

template <class T, int Dim, bool Packed, int... Indices>
struct scalar_type<Swizzle<T, Dim, Packed, Indices...>> {
	using type = T;
};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct scalar_type<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	using type = T;
};

template <class T, eQuaternionLayout Layout, bool Packed>
struct scalar_type<Quaternion<T, Layout, Packed>> {
	using type = T;
};

template <class T>
using scalar_type_t = typename scalar_type<T>::type;


//--------------------------------------
// Dimension
//--------------------------------------

template <class T>
struct dimension {
	static constexpr int value = 1;
};

template <class T, int Dim, bool Packed>
struct dimension<Vector<T, Dim, Packed>> {
	static constexpr int value = Dim;
};

template <class T, int Dim, bool Packed, int... Indices>
struct dimension<Swizzle<T, Dim, Packed, Indices...>> {
	static constexpr int value = sizeof...(Indices);
};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct dimension<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	template <int Rows_, int Columns_>
	static constexpr int Get() {
		static_assert(Rows_ == 1 || Columns_ == 1, "Must be a row or column matrix to determine dimension.");
		return Rows_ == 1 ? Columns_ : Rows_;
	}
	static constexpr int value = Get<Rows, Columns>();
};


template <class T>
constexpr auto dimension_v = dimension<T>::value;


template <class T>
struct source_dimension {};


template <class T, int Dim, bool Packed, int... Indices>
struct source_dimension<Swizzle<T, Dim, Packed, Indices...>> {
	static constexpr int value = Dim;
};


template <class T>
constexpr auto source_dimension_v = source_dimension<T>::value;


template <class T>
struct target_dimension {};


template <class T, int Dim, bool Packed, int... Indices>
struct target_dimension<Swizzle<T, Dim, Packed, Indices...>> {
	static constexpr int value = sizeof...(Indices);
};


template <class T>
constexpr auto target_dimension_v = target_dimension<T>::value;


//--------------------------------------
// Rows
//--------------------------------------

template <class T>
struct row_count {};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct row_count<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	static constexpr auto value = Rows;
};

template <class T>
constexpr auto row_count_v = row_count<T>::value;

//--------------------------------------
// Columns
//--------------------------------------

template <class T>
struct column_count {};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct column_count<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	static constexpr auto value = Columns;
};

template <class T>
constexpr auto column_count_v = column_count<T>::value;


//--------------------------------------
// Order
//--------------------------------------

template <class T>
struct order {};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct order<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	static constexpr auto value = Order;
};

template <class T>
constexpr auto order_v = order<T>::value;


//--------------------------------------
// Layout
//--------------------------------------

template <class T>
struct layout {};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct layout<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	static constexpr auto value = Layout;
};

template <class T, eQuaternionLayout Layout, bool Packed>
struct layout<Quaternion<T, Layout, Packed>> {
	static constexpr auto value = Layout;
};

template <class T>
constexpr auto layout_v = layout<T>::value;


//--------------------------------------
// Packed
//--------------------------------------

template <class T>
struct is_packed {};

template <class T, int Dim, bool Packed>
struct is_packed<Vector<T, Dim, Packed>> {
	static constexpr auto value = Packed;
};

template <class T, int Dim, bool Packed, int... Indices>
struct is_packed<Swizzle<T, Dim, Packed, Indices...>> {
	static constexpr auto value = Packed;
};

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
struct is_packed<Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	static constexpr auto value = Packed;
};

template <class T, eQuaternionLayout Layout, bool Packed>
struct is_packed<Quaternion<T, Layout, Packed>> {
	static constexpr auto value = Packed;
};

template <class T>
constexpr auto is_packed_v = is_packed<T>::value;


//------------------------------------------------------------------------------
// Matrix layout and order
//------------------------------------------------------------------------------

template <eMatrixOrder Order>
struct opposite_order {
	static constexpr auto value = Order == eMatrixOrder::FOLLOW_VECTOR ? eMatrixOrder::PRECEDE_VECTOR : eMatrixOrder::FOLLOW_VECTOR;
};

template <eMatrixOrder Order>
constexpr auto opposite_order_v = opposite_order<Order>::value;


template <eMatrixLayout Layout>
struct opposite_layout {
	static constexpr auto value = Layout == eMatrixLayout::ROW_MAJOR ? eMatrixLayout::COLUMN_MAJOR : eMatrixLayout::ROW_MAJOR;
};

template <eMatrixLayout Layout>
constexpr auto opposite_layout_v = opposite_layout<Layout>::value;


//------------------------------------------------------------------------------
// Complex
//------------------------------------------------------------------------------

template <class T>
struct remove_complex {
	using type = T;
};

template <class T>
struct remove_complex<std::complex<T>> {
	using type = T;
};

template <class T>
using remove_complex_t = typename remove_complex<T>::type;


template <class T>
struct is_complex : std::false_type {};

template <class T>
struct is_complex<std::complex<T>> : std::true_type {};

template <class T>
inline constexpr bool is_complex_v = is_complex<T>::value;


//------------------------------------------------------------------------------
// Common arithmetic type
//------------------------------------------------------------------------------

namespace impl {

	template <class... T>
	struct common_arithmetic_type;

	template <>
	struct common_arithmetic_type<> {
		using type = void;
		static constexpr bool valid = false;
	};

	template <class T>
	struct common_arithmetic_type<T> {
		using type = T;
		static constexpr bool valid = true;
	};

	template <class T1, class T2, class... Rest>
	struct common_arithmetic_type<T1, T2, Rest...> {
		template <class U1, class U2, class R = std::invoke_result_t<std::plus<>, U1, U2>>
		static constexpr R get_type(U1&&, U2&&);

		static constexpr void get_type(...);

		template <class U1, class U2, class R = std::invoke_result_t<std::plus<>, U1, U2>>
		static constexpr std::true_type is_valid(U1&&, U2&&);

		static constexpr std::false_type is_valid(...);

		using type = typename common_arithmetic_type<decltype(common_arithmetic_type::get_type(std::declval<T1>(), std::declval<T2>())), Rest...>::type;
		static constexpr bool valid = decltype(common_arithmetic_type::is_valid(std::declval<T1>(), std::declval<T2>()))::value && common_arithmetic_type<T2, Rest...>::valid;
	};

} // namespace impl


template <class... T>
struct common_arithmetic_type : std::enable_if<impl::common_arithmetic_type<T...>::valid, typename impl::common_arithmetic_type<T...>::type> {};

#ifdef MATHTER_ENABLE_SIMD
namespace impl {
	template <class Batch>
	struct common_arithmetic_type_xsimd : std::enable_if<!std::is_void_v<Batch>, Batch> {};
} // namespace impl

template <class... Ts, class... As>
struct common_arithmetic_type<xsimd::batch<Ts, As>...>
	: impl::common_arithmetic_type_xsimd<
		  xsimd::make_sized_batch_t<typename common_arithmetic_type<Ts...>::type,
									(... + xsimd::batch<Ts, As>::size) / sizeof...(Ts)>> {};
#endif


template <class... T>
using common_arithmetic_type_t = typename common_arithmetic_type<T...>::type;

// Tests
static_assert(std::is_same_v<common_arithmetic_type_t<float, double>, double>);
static_assert(std::is_same_v<common_arithmetic_type_t<double, float>, double>);


//------------------------------------------------------------------------------
// Sized integer
//------------------------------------------------------------------------------

template <size_t Bytes, bool Signed>
struct make_sized_integer {};

template <bool Signed>
struct make_sized_integer<1, Signed> {
	using type = std::conditional_t<Signed, int8_t, uint8_t>;
};

template <bool Signed>
struct make_sized_integer<2, Signed> {
	using type = std::conditional_t<Signed, int16_t, uint16_t>;
};

template <bool Signed>
struct make_sized_integer<4, Signed> {
	using type = std::conditional_t<Signed, int32_t, uint32_t>;
};

template <bool Signed>
struct make_sized_integer<8, Signed> {
	using type = std::conditional_t<Signed, int64_t, uint64_t>;
};

template <size_t Bytes, bool Signed>
using make_sized_integer_t = typename make_sized_integer<Bytes, Signed>::type;


} // namespace mathter