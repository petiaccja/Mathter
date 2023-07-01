// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include <Mathter/Common/Definitions.hpp>
#include <Mathter/Common/Traits.hpp>

#include <complex>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <typeinfo>



// Case combination elements
template <class... Types>
struct TypeCases {};

template <mathter::eMatrixOrder... Orders>
struct OrderCases {};

template <mathter::eMatrixLayout... Layouts>
struct LayoutCases {};

template <bool... Packed>
struct PackingCases {};



// Commonly used cases
using TypesAll = TypeCases<float, double, int, std::complex<float>>;
using TypesReal = TypeCases<float, double, int>;
using TypesFloating = TypeCases<float, double>;
using TypesInt = TypeCases<int>;

using OrdersAll = OrderCases<mathter::eMatrixOrder::FOLLOW_VECTOR, mathter::eMatrixOrder::PRECEDE_VECTOR>;
using OrdersFollow = OrderCases<mathter::eMatrixOrder::FOLLOW_VECTOR>;

using LayoutsAll = LayoutCases<mathter::eMatrixLayout::ROW_MAJOR, mathter::eMatrixLayout::COLUMN_MAJOR>;
using LayoutsRow = LayoutCases<mathter::eMatrixLayout::ROW_MAJOR>;
using LayoutsCol = LayoutCases<mathter::eMatrixLayout::COLUMN_MAJOR>;

using PackedAll = PackingCases<true, false>;
using PackedFalse = PackingCases<false>;
using PackedTrue = PackingCases<true>;



namespace impl {

template <class Type, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool Packed>
struct TypeMaker {
	template <int Dim>
	using Vector = mathter::Vector<Type, Dim, Packed>;

	template <int Rows, int Columns>
	using Matrix = mathter::Matrix<Type, Rows, Columns, Order, Layout, Packed>;

	using Quat = mathter::Quaternion<Type, Packed>;

	static const char* Name() {
		static const std::string name = []() {
			std::stringstream ss;

			ss << "(";
			ss << typeid(Type).name() << ", ";
			ss << (Packed ? "PACKED" : "VECTORIZED") << ", ";
			ss << (Order == mathter::eMatrixOrder::FOLLOW_VECTOR ? "FOLLOW_VECTOR" : "PRECEDE_VECTOR") << ", ";
			ss << (Layout == mathter::eMatrixLayout::ROW_MAJOR ? "ROW_MAJOR" : "COLUMN_MAJOR");
			ss << ")";

			return ss.str();
		}();
		return name.c_str();
	}
};


template <class Type, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, class PackingCasesT>
struct expand_packed;

template <class Type, mathter::eMatrixOrder Order, class LayoutCasesT, class PackingCasesT>
struct expand_layouts;

template <class Type, class OrderCasesT, class LayoutCasesT, class PackingCasesT>
struct expand_orders;

template <class TypeCasesT, class OrderCasesT, class LayoutCasesT, class PackingCasesT>
struct expand_types;


template <class Type, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool... Packed>
struct expand_packed<Type, Order, Layout, PackingCases<Packed...>> {
	using type = std::tuple<TypeMaker<Type, Order, Layout, Packed>...>;
};

template <class Type, mathter::eMatrixOrder Order, mathter::eMatrixLayout... Layouts, bool... Packed>
struct expand_layouts<Type, Order, LayoutCases<Layouts...>, PackingCases<Packed...>> {
	using type = decltype(std::tuple_cat(
		std::declval<typename expand_packed<Type,
											Order,
											Layouts,
											PackingCases<Packed...>>::type>()...));
};

template <class Type, mathter::eMatrixOrder... Orders, mathter::eMatrixLayout... Layouts, bool... Packed>
struct expand_orders<Type, OrderCases<Orders...>, LayoutCases<Layouts...>, PackingCases<Packed...>> {
	using type = decltype(std::tuple_cat(
		std::declval<typename expand_layouts<Type,
											 Orders,
											 LayoutCases<Layouts...>,
											 PackingCases<Packed...>>::type>()...));
};

template <class... Types, mathter::eMatrixOrder... Orders, mathter::eMatrixLayout... Layouts, bool... Packed>
struct expand_types<TypeCases<Types...>, OrderCases<Orders...>, LayoutCases<Layouts...>, PackingCases<Packed...>> {
	using type = decltype(std::tuple_cat(
		std::declval<typename expand_orders<Types,
											OrderCases<Orders...>,
											LayoutCases<Layouts...>,
											PackingCases<Packed...>>::type>()...));
};

template <class TypeCasesT, class PackingCasesT, class OrderCasesT, class LayoutCasesT>
struct TestTypeList;

template <class... Types, bool... Packed, mathter::eMatrixOrder... Orders, mathter::eMatrixLayout... Layouts>
struct TestTypeList<TypeCases<Types...>, PackingCases<Packed...>, OrderCases<Orders...>, LayoutCases<Layouts...>> {
	using type = typename expand_types<TypeCases<Types...>, OrderCases<Orders...>, LayoutCases<Layouts...>, PackingCases<Packed...>>::type;
};


} // namespace impl


template <class TypeCasesT,
		  class PackingCasesT,
		  class OrderCasesT = OrderCases<mathter::eMatrixOrder::FOLLOW_VECTOR>,
		  class LayoutCasesT = LayoutCases<mathter::eMatrixLayout::ROW_MAJOR>>
using TestTypeList = typename impl::TestTypeList<TypeCasesT, PackingCasesT, OrderCasesT, LayoutCasesT>::type;


// Invert order
template <class MatrixT>
struct invert_order;

template <class Type, int Rows, int Columns, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool Packed>
struct invert_order<mathter::Matrix<Type, Rows, Columns, Order, Layout, Packed>> {
	static constexpr auto order = Order == mathter::eMatrixOrder::FOLLOW_VECTOR ?
									  mathter::eMatrixOrder::PRECEDE_VECTOR :
									  mathter::eMatrixOrder::FOLLOW_VECTOR;
	using type = mathter::Matrix<Type, Columns, Rows, order, Layout, Packed>;
};

template <class MatrixT>
using invert_order_t = typename invert_order<MatrixT>::type;


// Invert layout
template <class MatrixT>
struct invert_layout;

template <class Type, int Rows, int Columns, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool Packed>
struct invert_layout<mathter::Matrix<Type, Rows, Columns, Order, Layout, Packed>> {
	static constexpr auto layout = Layout == mathter::eMatrixLayout::COLUMN_MAJOR ?
									   mathter::eMatrixLayout::ROW_MAJOR :
									   mathter::eMatrixLayout::COLUMN_MAJOR;
	using type = mathter::Matrix<Type, Rows, Columns, Order, layout, Packed>;
};

template <class MatrixT>
using invert_layout_t = typename invert_layout<MatrixT>::type;
