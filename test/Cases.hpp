#pragma once


#include <Mathter/Common/Types.hpp>

#include <complex>
#include <type_traits>


#define CASE_LIST(X) decltype(X{})


namespace test_util {

//------------------------------------------------------------------------------
// Template argument list definitions
//------------------------------------------------------------------------------

template <class... Arguments>
struct TemplateArgumentList {};


template <class T, T... Arguments>
struct NTTemplateArgumentList {};


//------------------------------------------------------------------------------
// Common argument lists
//------------------------------------------------------------------------------

using ScalarsFloat32 = TemplateArgumentList<float>;
using ScalarsInt32 = TemplateArgumentList<int32_t>;
using ScalarsFloatAndInt32 = TemplateArgumentList<float, int32_t>;
using ScalarsFloatAndInt32AndComplex = TemplateArgumentList<float, int32_t, std::complex<float>>;
using ScalarsFloating = TemplateArgumentList<float, double>;
using ScalarsIntegral = TemplateArgumentList<int32_t, int64_t>;
using ScalarsComplex = TemplateArgumentList<std::complex<float>, std::complex<double>>;
using ScalarsComplex32 = TemplateArgumentList<std::complex<float>>;
using ScalarsFloatingAndComplex = TemplateArgumentList<float, double, std::complex<float>, std::complex<double>>;
using ScalarsFloatingAndComplex32 = TemplateArgumentList<float, std::complex<float>>;
using ScalarsAll = TemplateArgumentList<float, double, int32_t, int64_t, std::complex<float>, std::complex<double>>;

using OrdersFollow = NTTemplateArgumentList<mathter::eMatrixOrder, mathter::eMatrixOrder::FOLLOW_VECTOR>;
using OrdersPrecede = NTTemplateArgumentList<mathter::eMatrixOrder, mathter::eMatrixOrder::PRECEDE_VECTOR>;
using OrdersAll = NTTemplateArgumentList<mathter::eMatrixOrder, mathter::eMatrixOrder::PRECEDE_VECTOR, mathter::eMatrixOrder::FOLLOW_VECTOR>;

using LayoutsRM = NTTemplateArgumentList<mathter::eMatrixLayout, mathter::eMatrixLayout::ROW_MAJOR>;
using LayoutsCM = NTTemplateArgumentList<mathter::eMatrixLayout, mathter::eMatrixLayout::COLUMN_MAJOR>;
using LayoutsAll = NTTemplateArgumentList<mathter::eMatrixLayout, mathter::eMatrixLayout::COLUMN_MAJOR, mathter::eMatrixLayout::ROW_MAJOR>;

using PackingsYes = NTTemplateArgumentList<bool, true>;
using PackingsNo = NTTemplateArgumentList<bool, false>;
using PackingsAll = NTTemplateArgumentList<bool, false, true>;


//------------------------------------------------------------------------------
// Vector cases
//------------------------------------------------------------------------------

template <class Scalar, bool Packed>
struct VectorCase {
	template <int Dim>
	using Vector = mathter::Vector<Scalar, Dim, Packed>;
};


template <class ScalarOptions, class PackingOptions>
struct VectorCaseList_Builder;


template <class... Scalars, bool... Packings>
struct VectorCaseList_Builder<TemplateArgumentList<Scalars...>, NTTemplateArgumentList<bool, Packings...>> {
	template <class Scalar>
	static auto ExpandPackings() -> std::tuple<VectorCase<Scalar, Packings>...>;
	static auto ExpandScalars() -> decltype(std::tuple_cat(ExpandPackings<Scalars>()...));

	using Cases = decltype(ExpandScalars());
};


template <class ScalarOptions, class PackingOptions>
using VectorCaseList = typename VectorCaseList_Builder<ScalarOptions, PackingOptions>::Cases;


//------------------------------------------------------------------------------
// Matrix cases
//------------------------------------------------------------------------------

template <class Scalar, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool Packed>
struct MatrixCase {
	template <int Rows, int Columns>
	using Matrix = mathter::Matrix<Scalar, Rows, Columns, Order, Layout, Packed>;
};


template <class ScalarOptions, class OrderOptions, class LayoutOptions, class PackingOptions>
struct MatrixCaseList_Builder;


template <class... Scalars, mathter::eMatrixOrder... Orders, mathter::eMatrixLayout... Layouts, bool... Packings>
struct MatrixCaseList_Builder<TemplateArgumentList<Scalars...>,
							  NTTemplateArgumentList<mathter::eMatrixOrder, Orders...>,
							  NTTemplateArgumentList<mathter::eMatrixLayout, Layouts...>,
							  NTTemplateArgumentList<bool, Packings...>> {
	template <class Scalar, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout>
	static auto ExpandPackings() -> std::tuple<MatrixCase<Scalar, Order, Layout, Packings>...>;

	template <class Scalar, mathter::eMatrixOrder Order>
	static auto ExpandLayouts() -> decltype(std::tuple_cat(ExpandPackings<Scalar, Order, Layouts>()...));

	template <class Scalar>
	static auto ExpandOrders() -> decltype(std::tuple_cat(ExpandLayouts<Scalar, Orders>()...));

	static auto ExpandScalars() -> decltype(std::tuple_cat(ExpandOrders<Scalars>()...));

	using Cases = decltype(ExpandScalars());
};


template <class ScalarOptions, class OrderOptions, class LayoutOptions, class PackingOptions>
using MatrixCaseList = typename MatrixCaseList_Builder<ScalarOptions, OrderOptions, LayoutOptions, PackingOptions>::Cases;


//------------------------------------------------------------------------------
// Quaternion cases
//------------------------------------------------------------------------------

template <class Scalar, bool Packed>
struct QuaternionCase {
	using Quat = mathter::Quaternion<Scalar, Packed>;
};


template <class ScalarOptions, class PackingOptions>
struct QuaternionCaseList_Builder;


template <class... Scalars, bool... Packings>
struct QuaternionCaseList_Builder<TemplateArgumentList<Scalars...>, NTTemplateArgumentList<bool, Packings...>> {
	template <class Scalar>
	static auto ExpandPackings() -> std::tuple<QuaternionCase<Scalar, Packings>...>;
	static auto ExpandScalars() -> decltype(std::tuple_cat(ExpandPackings<Scalars>()...));

	using Cases = decltype(ExpandScalars());
};


template <class ScalarOptions, class PackingOptions>
using QuaternionCaseList = typename QuaternionCaseList_Builder<ScalarOptions, PackingOptions>::Cases;


//------------------------------------------------------------------------------
// Swizzle cases
//------------------------------------------------------------------------------

template <class Scalar, bool Packed>
struct SwizzleCase {
	template <int Dim, int... Indices>
	using Swizzle = mathter::Swizzle<Scalar, Dim, Packed, Indices...>;
};


template <class ScalarOptions, class PackingOptions>
struct SwizzleCaseList_Builder;


template <class... Scalars, bool... Packings>
struct SwizzleCaseList_Builder<TemplateArgumentList<Scalars...>, NTTemplateArgumentList<bool, Packings...>> {
	template <class Scalar>
	static auto ExpandPackings() -> std::tuple<SwizzleCase<Scalar, Packings>...>;
	static auto ExpandScalars() -> decltype(std::tuple_cat(ExpandPackings<Scalars>()...));

	using Cases = decltype(ExpandScalars());
};


template <class ScalarOptions, class PackingOptions>
using SwizzleCaseList = typename SwizzleCaseList_Builder<ScalarOptions, PackingOptions>::Cases;


//------------------------------------------------------------------------------
// Scalar cases
//------------------------------------------------------------------------------

template <class Scalar_>
struct ScalarCase {
	using Scalar = Scalar_;
};


template <class ScalarOptions>
struct ScalarCaseList_Builder;


template <class... Scalars>
struct ScalarCaseList_Builder<TemplateArgumentList<Scalars...>> {
	static auto ExpandScalars() -> std::tuple<ScalarCase<Scalars>...>;

	using Cases = decltype(ExpandScalars());
};


template <class ScalarOptions>
using ScalarCaseList = typename ScalarCaseList_Builder<ScalarOptions>::Cases;


//------------------------------------------------------------------------------
// Product cases
//------------------------------------------------------------------------------

template <class LhsCase, class RhsCase>
struct BinaryCase {
	using Lhs = LhsCase;
	using Rhs = RhsCase;
};


template <class LhsCaseList, class RhsCaseList>
struct BinaryCaseList_Builder;


template <class... LhsCases, class... RhsCases>
struct BinaryCaseList_Builder<std::tuple<LhsCases...>, std::tuple<RhsCases...>> {
	template <class LhsCase>
	static auto ExpandRhs() -> std::tuple<BinaryCase<LhsCase, RhsCases>...>;
	static auto ExpandLhs() -> decltype(std::tuple_cat(ExpandRhs<LhsCases>()...));

	using Cases = decltype(ExpandLhs());
};


template <class LhsCaseList, class RhsCaseList>
using BinaryCaseList = typename BinaryCaseList_Builder<LhsCaseList, RhsCaseList>::Cases;

} // namespace test_util