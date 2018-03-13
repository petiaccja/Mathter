#pragma once

#include "Mathter/DefinitionsUtil.hpp"

#include <type_traits>
#include <typeinfo>
#include <iostream>
#include <complex>



// Get n-th type of a typelist
template <unsigned Index, class... Rest>
struct GetType;

template <unsigned Index, class Head, class... Rest>
struct GetType<Index, Head, Rest...> {
	using type = typename GetType<Index - 1, Rest...>::type;
};

template <class Head, class... Rest>
struct GetType<0, Head, Rest...> {
	using type = Head;
};


// Get n-th value of a value list
template <unsigned Index, class Type, Type... Values>
struct GetValue;

template <unsigned Index, class Type, Type Head, Type... Rest>
struct GetValue<Index, Type, Head, Rest...> {
	static constexpr Type value = GetValue<Index - 1, Type, Rest...>::value;
};

template <class Type, Type Head, Type... Rest>
struct GetValue<0, Type, Head, Rest...> {
	static constexpr Type value = Head;
};


// Case combination elements
template <class... Types>
struct TypeCases {
public:
	template <int Index>
	using type = typename GetType<Index, Types...>::type;
};

template <mathter::eMatrixOrder... Orders>
struct OrderCases {
	template <int Index>
	static constexpr mathter::eMatrixOrder value = GetValue<Index, mathter::eMatrixOrder, Orders...>::value;
};

template <mathter::eMatrixLayout... Layouts>
struct LayoutCases {
	template <int Index>
	static constexpr mathter::eMatrixLayout value = GetValue<Index, mathter::eMatrixLayout, Layouts...>::value;
};

template <bool... Packed>
struct PackingCases {
	template <int Index>
	static constexpr bool value = GetValue<Index, bool, Packed...>::value;
};



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



// Declare a combination of simple cases
template <class TypeCasesT, class OrderCasesT, class LayoutCasesT, class PackingCasesT>
struct MatrixCases;

template <class... Types, mathter::eMatrixOrder... Orders, mathter::eMatrixLayout... Layouts, bool... Packed>
struct MatrixCases<TypeCases<Types...>, OrderCases<Orders...>, LayoutCases<Layouts...>, PackingCases<Packed...>> {
	static constexpr unsigned NumTypes = sizeof...(Types);
	static constexpr unsigned NumOrders = sizeof...(Orders);
	static constexpr unsigned NumLayouts = sizeof...(Layouts);
	static constexpr unsigned NumPackings = sizeof...(Packed);
	static constexpr unsigned NumCombinations = NumTypes * NumOrders * NumLayouts * NumPackings;

	static constexpr mathter::eMatrixOrder order = OrderCases<Orders...>::value<0>;
	static constexpr mathter::eMatrixLayout layout = LayoutCases<Layouts...>::value<0>;

	using OrderCasesT = OrderCases<Orders...>;

	template <int Index>
	using MatrixT = mathter::Matrix<
		typename TypeCases<Types...>::type<Index / (NumOrders * NumLayouts * NumPackings)>,
		1, 1,
		GetValue<Index % (NumOrders * NumLayouts * NumPackings) / (NumLayouts * NumPackings), mathter::eMatrixOrder, Orders...>::value,
		GetValue<Index % (NumOrders * NumLayouts * NumPackings) % (NumLayouts * NumPackings) / NumPackings, mathter::eMatrixLayout, Layouts...>::value,
		GetValue<Index % (NumOrders * NumLayouts * NumPackings) % (NumLayouts * NumPackings) % NumPackings, bool, Packed...>::value
	>;

	template <int Index>
	using Properties = mathter::impl::MatrixProperties<MatrixT<Index>>;
};



// Helper that iterates over combinations
template <int Index, class Cases, template <typename, mathter::eMatrixOrder, mathter::eMatrixLayout, bool> class Functor>
int RunCasesHelper() {
	if (Index > 0) {
		Functor<
			typename Cases::Properties<Index>::Type,
			Cases::Properties<Index>::Order,
			Cases::Properties<Index>::Layout,
			Cases::Properties<Index>::Packed
		>()();
		return 1 + RunCasesHelper<std::max(0, Index - 1), Cases, Functor>();
	}
	else {
		Functor<
			typename Cases::Properties<Index>::Type,
			Cases::Properties<Index>::Order,
			Cases::Properties<Index>::Layout,
			Cases::Properties<Index>::Packed
		>()();
		return 1;
	}
}



// Run all combinations
template <class Cases, template <class, mathter::eMatrixOrder, mathter::eMatrixLayout, bool> class Functor>
int RunCases() {
	return RunCasesHelper<(int)Cases::NumCombinations - 1, Cases, Functor>();
}



template <
	class Type1,
	mathter::eMatrixOrder Order1,
	mathter::eMatrixLayout Layout1,
	bool Packed1,
	template <class, mathter::eMatrixOrder, mathter::eMatrixLayout, bool, class, mathter::eMatrixOrder, mathter::eMatrixLayout, bool> class BiFunctor
>
struct UniFunctorify {
	template<
		class Type2,
		mathter::eMatrixOrder Order2,
		mathter::eMatrixLayout Layout2,
		bool Packed2
	>
	struct UniFunctor {
		void operator()() const {
			BiFunctor<Type1, Order1, Layout1, Packed1, Type2, Order2, Layout2, Packed2>()();
		}
	};
};

template <
	int Index,
	class Cases1,
	class Cases2,
	template <class, mathter::eMatrixOrder, mathter::eMatrixLayout, bool, class, mathter::eMatrixOrder, mathter::eMatrixLayout, bool> class BiFunctor
>
int RunCasesHelper() {
	if (Index > 0) {
		// Call unifunctor helper
		int count = RunCasesHelper<
			Cases2::NumCombinations - 1,
			Cases2,
			typename UniFunctorify<
				typename Cases1::Properties<Index>::Type,
				Cases1::Properties<Index>::Order,
				Cases1::Properties<Index>::Layout,
				Cases1::Properties<Index>::Packed,
				BiFunctor
			>::UniFunctor
		>();
		// Recurse
		return count + RunCasesHelper<std::max(0, Index - 1), Cases1, Cases2, BiFunctor>();
	}
	else {
		// Call unifunctor helper
		return RunCasesHelper<
			Cases2::NumCombinations - 1,
			Cases2,
			typename UniFunctorify<
				typename Cases1::Properties<Index>::Type,
				Cases1::Properties<Index>::Order,
				Cases1::Properties<Index>::Layout,
				Cases1::Properties<Index>::Packed,
				BiFunctor
			>::UniFunctor
		>();
	}
}


// Run all combinations x tensor product
template <
	class Cases1,
	class Cases2,
	template <class, mathter::eMatrixOrder, mathter::eMatrixLayout, bool, class, mathter::eMatrixOrder, mathter::eMatrixLayout, bool> class BiFunctor
>
int RunCases() {
	return RunCasesHelper<Cases1::NumCombinations-1, Cases1, Cases2, BiFunctor>();
}



// Macro helper for easily declaring templated test cases
#define GEN_FUNC_NAME_HELPER2(Name, Cntr) Name ## Cntr
#define GEN_FUNC_NAME_HELPER1(Name, Cntr) GEN_FUNC_NAME_HELPER2(Name, Cntr)
#define GEN_FUNC_NAME GEN_FUNC_NAME_HELPER1(_GenFunc_, __COUNTER__)


#define TEST_CASE_VARIANT_H(NAME, TAG, TYPES, ORDERS, LAYOUTS, PACKINGS, FUNNAME)	\
template <class Type, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>	\
struct FUNNAME {																\
	template<int Rows, int Columns>												\
	using MatrixT = Matrix<Type, Rows, Columns, Order, Layout, Packed>;			\
	void operator()() const;													\
};																				\
TEST_CASE(NAME, TAG) {															\
	using Cases = MatrixCases<													\
		TYPES,																	\
		ORDERS,																	\
		LAYOUTS,																\
		PACKINGS>;																\
																				\
	RunCases<Cases, FUNNAME>();													\
}																				\
template <class Type, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>	\
void FUNNAME<Type, Order, Layout, Packed>::operator()() const													

#define TEST_CASE_VARIANT(NAME, TAG, TYPES, ORDERS, LAYOUTS, PACKINGS) TEST_CASE_VARIANT_H(NAME, TAG, TYPES, ORDERS, LAYOUTS, PACKINGS, GEN_FUNC_NAME) 


#define TEST_CASE_VARIANT_H_2(NAME, TAG, TYPES1, ORDERS1, LAYOUTS1, PACKINGS1, TYPES2, ORDERS2, LAYOUTS2, PACKINGS2, FUNNAME)	\
template <class Type1, eMatrixOrder Order1, eMatrixLayout Layout1, bool Packed1,\
		  class Type2, eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2>\
struct FUNNAME {																\
	template<int Rows, int Columns>												\
	using MatrixT1 = Matrix<Type1, Rows, Columns, Order1, Layout1, Packed1>;	\
	template<int Rows, int Columns>												\
	using MatrixT2 = Matrix<Type2, Rows, Columns, Order2, Layout2, Packed2>;	\
	void operator()() const;													\
};																				\
TEST_CASE(NAME, TAG) {															\
	using Cases1 = MatrixCases<													\
		TYPES1,																	\
		ORDERS1,																\
		LAYOUTS1,																\
		PACKINGS1>;																\
	using Cases2 = MatrixCases<													\
		TYPES2,																	\
		ORDERS2,																\
		LAYOUTS2,																\
		PACKINGS2>;																\
																				\
	RunCases<Cases1, Cases2, FUNNAME>();										\
}																				\
template <class Type1, eMatrixOrder Order1, eMatrixLayout Layout1, bool Packed1,\
		  class Type2, eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2>\
void FUNNAME<Type1, Order1, Layout1, Packed1, Type2, Order2, Layout2, Packed2>::operator()() const						


#define TEST_CASE_VARIANT_2(NAME, TAG, TYPES1, ORDERS1, LAYOUTS1, PACKINGS1, TYPES2, ORDERS2, LAYOUTS2, PACKINGS2) \
TEST_CASE_VARIANT_H_2(NAME, TAG, TYPES1, ORDERS1, LAYOUTS1, PACKINGS1, TYPES2, ORDERS2, LAYOUTS2, PACKINGS2, GEN_FUNC_NAME)

template <class Type, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool Packed>
const std::string& SectionName() {
	std::stringstream ss;

	ss << typeid(Type).name() << ", ";
	ss << (Order == mathter::eMatrixOrder::FOLLOW_VECTOR ? "FOLLOW_VECTOR" : "PRECEDE_VECTOR") << ", ";
	ss << (Layout == mathter::eMatrixLayout::ROW_MAJOR ? "ROW_MAJOR" : "COLUMN_MAJOR") << ", ";
	ss << (Packed ? "true" : "false");

	thread_local std::string str; // Ugly hack but who cares?
	str = ss.str();
	return str;
}

template <
	class Type1, mathter::eMatrixOrder Order1, mathter::eMatrixLayout Layout1, bool Packed1,
	class Type2, mathter::eMatrixOrder Order2, mathter::eMatrixLayout Layout2, bool Packed2>
const std::string& SectionName2() {
	std::stringstream ss;

	ss << typeid(Type1).name() << ", ";
	ss << (Order1 == mathter::eMatrixOrder::FOLLOW_VECTOR ? "FOLLOW_VECTOR" : "PRECEDE_VECTOR") << ", ";
	ss << (Layout1 == mathter::eMatrixLayout::ROW_MAJOR ? "ROW_MAJOR" : "COLUMN_MAJOR") << ", ";
	ss << (Packed1 ? "true" : "false");

	ss << " ? ";

	ss << typeid(Type2).name() << ", ";
	ss << (Order2 == mathter::eMatrixOrder::FOLLOW_VECTOR ? "FOLLOW_VECTOR" : "PRECEDE_VECTOR") << ", ";
	ss << (Layout2 == mathter::eMatrixLayout::ROW_MAJOR ? "ROW_MAJOR" : "COLUMN_MAJOR") << ", ";
	ss << (Packed2 ? "true" : "false");

	thread_local std::string str; // Ugly hack but who cares?
	str = ss.str();
	return str;
}
#define SECTIONNAME SectionName<Type, Order, Layout, Packed>().c_str()
#define SECTIONNAME2 SectionName2<Type1, Order1, Layout1, Packed1, Type2, Order2, Layout2, Packed2>().c_str()



// Prints the type of a matrix in human-readable form
template <class MatrixT, std::enable_if_t<mathter::impl::IsMatrix<MatrixT>::value, int> = 0>
std::string PrintType() {
	std::stringstream ss;

	ss << "Matrix<";
	ss << typeid(typename impl::MatrixProperties<MatrixT>::Type).name() << ", ";
	ss << impl::MatrixProperties<MatrixT>::Rows << ", ";
	ss << impl::MatrixProperties<MatrixT>::Columns << ", ";
	ss << (impl::MatrixProperties<MatrixT>::Order == mathter::eMatrixOrder::FOLLOW_VECTOR ? "FOLLOW_VECTOR" : "PRECEDE_VECTOR") << ", ";
	ss << (impl::MatrixProperties<MatrixT>::Layout == mathter::eMatrixLayout::ROW_MAJOR ? "ROW_MAJOR" : "COLUMN_MAJOR") << ", ";
	ss << (impl::MatrixProperties<MatrixT>::Packed ? "true" : "false");
	ss << ">";

	return ss.str();
}
