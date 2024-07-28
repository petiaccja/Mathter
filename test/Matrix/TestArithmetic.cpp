// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Cases.hpp"

#include <Mathter/Matrix/Arithmetic.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


//------------------------------------------------------------------------------
// Multiplication
//------------------------------------------------------------------------------

template <class MatLhs, class MatRhs>
static void TestMultiplication() {
	const MatLhs a = {
		1, 2, 3,
		4, 5, 6
	};

	const MatRhs b = {
		4, 3,
		3, 2,
		2, 1
	};

	const auto r = a * b;
	using Scalar = scalar_type_t<std::decay_t<decltype(r)>>;

	REQUIRE(r(0, 0) == static_cast<Scalar>(16));
	REQUIRE(r(0, 1) == static_cast<Scalar>(10));
	REQUIRE(r(1, 0) == static_cast<Scalar>(43));
	REQUIRE(r(1, 1) == static_cast<Scalar>(28));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiplication (follow)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, OrdersFollow, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloatingAndComplex32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {

	using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;
	using MatRhs = typename TestType::Rhs::template Matrix<3, 2>;

	TestMultiplication<MatLhs, MatRhs>();
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiplication (precede)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, OrdersPrecede, LayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloatingAndComplex32, OrdersPrecede, LayoutsAll, PackingsAll>>{})) {

	using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;
	using MatRhs = typename TestType::Rhs::template Matrix<3, 2>;

	TestMultiplication<MatLhs, MatRhs>();
}


//------------------------------------------------------------------------------
// Elementwise
//------------------------------------------------------------------------------

#define TEST_ELEMENTWISE(NAME, MAT_FUNC, SCALAR_FUNC, ORDERS)                                                                           \
	TEMPLATE_LIST_TEST_CASE("Matrix - " NAME " (" #ORDERS ")", "[Matrix]",                                                              \
							decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, ORDERS, LayoutsAll, PackingsAll>,       \
													MatrixCaseList<ScalarsFloatingAndComplex32, ORDERS, LayoutsAll, PackingsAll>>{})) { \
                                                                                                                                        \
		using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;                                                                   \
		using MatRhs = typename TestType::Rhs::template Matrix<2, 3>;                                                                   \
		const MatLhs a = {                                                                                                              \
			1, 2, 3,                                                                                                                    \
			4, 5, 6                                                                                                                     \
		};                                                                                                                              \
                                                                                                                                        \
		const MatRhs b = {                                                                                                              \
			4, 3, 3,                                                                                                                    \
			2, 2, 1                                                                                                                     \
		};                                                                                                                              \
                                                                                                                                        \
		const auto r = MAT_FUNC(a, b);                                                                                                  \
                                                                                                                                        \
		REQUIRE(r(0, 0) == SCALAR_FUNC(a(0, 0), b(0, 0)));                                                                              \
		REQUIRE(r(0, 1) == SCALAR_FUNC(a(0, 1), b(0, 1)));                                                                              \
		REQUIRE(r(0, 2) == SCALAR_FUNC(a(0, 2), b(0, 2)));                                                                              \
		REQUIRE(r(1, 0) == SCALAR_FUNC(a(1, 0), b(1, 0)));                                                                              \
		REQUIRE(r(1, 1) == SCALAR_FUNC(a(1, 1), b(1, 1)));                                                                              \
		REQUIRE(r(1, 2) == SCALAR_FUNC(a(1, 2), b(1, 2)));                                                                              \
	}


#define TEST_ELEMENTWISE_ASSIGN(NAME, OP, ORDERS)                                                                                       \
	TEMPLATE_LIST_TEST_CASE("Matrix - " NAME " (" #ORDERS ")", "[Matrix]",                                                              \
							decltype(BinaryCaseList<MatrixCaseList<ScalarsComplex32, ORDERS, LayoutsAll, PackingsAll>,                  \
													MatrixCaseList<ScalarsFloatingAndComplex32, ORDERS, LayoutsAll, PackingsAll>>{})) { \
                                                                                                                                        \
		using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;                                                                   \
		using MatRhs = typename TestType::Rhs::template Matrix<2, 3>;                                                                   \
                                                                                                                                        \
		const MatLhs a = {                                                                                                              \
			1, 2, 3,                                                                                                                    \
			4, 5, 6                                                                                                                     \
		};                                                                                                                              \
                                                                                                                                        \
		const MatRhs b = {                                                                                                              \
			4, 3, 3,                                                                                                                    \
			2, 2, 1                                                                                                                     \
		};                                                                                                                              \
                                                                                                                                        \
		auto copy = a;                                                                                                                  \
		REQUIRE(&(copy OP## = b) == &copy);                                                                                             \
                                                                                                                                        \
		REQUIRE(copy(0, 0) == a(0, 0) OP b(0, 0));                                                                                      \
		REQUIRE(copy(0, 1) == a(0, 1) OP b(0, 1));                                                                                      \
		REQUIRE(copy(0, 2) == a(0, 2) OP b(0, 2));                                                                                      \
		REQUIRE(copy(1, 0) == a(1, 0) OP b(1, 0));                                                                                      \
		REQUIRE(copy(1, 1) == a(1, 1) OP b(1, 1));                                                                                      \
		REQUIRE(copy(1, 2) == a(1, 2) OP b(1, 2));                                                                                      \
	}


TEST_ELEMENTWISE("Addition elementwise", std::plus{}, std::plus{}, OrdersFollow);
TEST_ELEMENTWISE("Addition elementwise", std::plus{}, std::plus{}, OrdersPrecede);
TEST_ELEMENTWISE("Subtraction elementwise", std::minus{}, std::minus{}, OrdersFollow);
TEST_ELEMENTWISE("Subtraction elementwise", std::minus{}, std::minus{}, OrdersPrecede);
TEST_ELEMENTWISE("Division elementwise", std::divides{}, std::divides{}, OrdersFollow);
TEST_ELEMENTWISE("Division elementwise", std::divides{}, std::divides{}, OrdersPrecede);
TEST_ELEMENTWISE("Hadamard product", Hadamard, std::multiplies{}, OrdersFollow);
TEST_ELEMENTWISE("Hadamard product", Hadamard, std::multiplies{}, OrdersPrecede);


TEST_ELEMENTWISE_ASSIGN("Addition elementwise assign", +, OrdersFollow);
TEST_ELEMENTWISE_ASSIGN("Addition elementwise assign", +, OrdersPrecede);
TEST_ELEMENTWISE_ASSIGN("Subtraction elementwise assign", -, OrdersFollow);
TEST_ELEMENTWISE_ASSIGN("Subtraction elementwise assign", -, OrdersPrecede);
TEST_ELEMENTWISE_ASSIGN("Division elementwise assign", /, OrdersFollow);
TEST_ELEMENTWISE_ASSIGN("Division elementwise assign", /, OrdersPrecede);


//------------------------------------------------------------------------------
// Scalar
//------------------------------------------------------------------------------


#define TEST_SCALAR(NAME, OP)                                                                                                        \
	TEMPLATE_LIST_TEST_CASE("Matrix - " NAME, "[Matrix]",                                                                            \
							decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, OrdersAll, LayoutsAll, PackingsAll>, \
													ScalarCaseList<ScalarsFloatingAndComplex32>>{})) {                               \
                                                                                                                                     \
		using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;                                                                \
		using ScalarRhs = typename TestType::Rhs::Scalar;                                                                            \
		const MatLhs a = {                                                                                                           \
			1, 2, 3,                                                                                                                 \
			4, 5, 6                                                                                                                  \
		};                                                                                                                           \
                                                                                                                                     \
		const ScalarRhs b = 4;                                                                                                       \
                                                                                                                                     \
		const auto r = a OP b;                                                                                                       \
                                                                                                                                     \
		REQUIRE(r(0, 0) == a(0, 0) OP b);                                                                                            \
		REQUIRE(r(0, 1) == a(0, 1) OP b);                                                                                            \
		REQUIRE(r(0, 2) == a(0, 2) OP b);                                                                                            \
		REQUIRE(r(1, 0) == a(1, 0) OP b);                                                                                            \
		REQUIRE(r(1, 1) == a(1, 1) OP b);                                                                                            \
		REQUIRE(r(1, 2) == a(1, 2) OP b);                                                                                            \
	}


#define TEST_SCALAR_REVERSE(NAME, OP)                                                                                                \
	TEMPLATE_LIST_TEST_CASE("Matrix - " NAME, "[Matrix]",                                                                            \
							decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, OrdersAll, LayoutsAll, PackingsAll>, \
													ScalarCaseList<ScalarsFloatingAndComplex32>>{})) {                               \
                                                                                                                                     \
		using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;                                                                \
		using ScalarRhs = typename TestType::Rhs::Scalar;                                                                            \
		const MatLhs a = {                                                                                                           \
			1, 2, 3,                                                                                                                 \
			4, 5, 6                                                                                                                  \
		};                                                                                                                           \
                                                                                                                                     \
		const ScalarRhs b = 4;                                                                                                       \
                                                                                                                                     \
		const auto r = b OP a;                                                                                                       \
                                                                                                                                     \
		REQUIRE(r(0, 0) == b OP a(0, 0));                                                                                            \
		REQUIRE(r(0, 1) == b OP a(0, 1));                                                                                            \
		REQUIRE(r(0, 2) == b OP a(0, 2));                                                                                            \
		REQUIRE(r(1, 0) == b OP a(1, 0));                                                                                            \
		REQUIRE(r(1, 1) == b OP a(1, 1));                                                                                            \
		REQUIRE(r(1, 2) == b OP a(1, 2));                                                                                            \
	}


#define TEST_SCALAR_ASSIGN(NAME, OP)                                                                                      \
	TEMPLATE_LIST_TEST_CASE("Matrix - " NAME, "[Matrix]",                                                                 \
							decltype(BinaryCaseList<MatrixCaseList<ScalarsComplex32, OrdersAll, LayoutsAll, PackingsAll>, \
													ScalarCaseList<ScalarsFloatingAndComplex32>>{})) {                    \
                                                                                                                          \
		using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;                                                     \
		using ScalarRhs = typename TestType::Rhs::Scalar;                                                                 \
		const MatLhs a = {                                                                                                \
			1, 2, 3,                                                                                                      \
			4, 5, 6                                                                                                       \
		};                                                                                                                \
                                                                                                                          \
		const ScalarRhs b = 4;                                                                                            \
                                                                                                                          \
		auto copy = a;                                                                                                    \
		REQUIRE(&(copy OP## = b) == &copy);                                                                               \
                                                                                                                          \
		REQUIRE(copy(0, 0) == a(0, 0) OP b);                                                                              \
		REQUIRE(copy(0, 1) == a(0, 1) OP b);                                                                              \
		REQUIRE(copy(0, 2) == a(0, 2) OP b);                                                                              \
		REQUIRE(copy(1, 0) == a(1, 0) OP b);                                                                              \
		REQUIRE(copy(1, 1) == a(1, 1) OP b);                                                                              \
		REQUIRE(copy(1, 2) == a(1, 2) OP b);                                                                              \
	}


TEST_SCALAR("Addition scalar", +);
TEST_SCALAR("Subtraction scalar", -);
TEST_SCALAR("Division scalar", /);
TEST_SCALAR("Multiplication scalar", *);


TEST_SCALAR_REVERSE("Addition scalar reverse", +);
TEST_SCALAR_REVERSE("Subtraction scalar reverse", -);
TEST_SCALAR_REVERSE("Division scalar reverse", /);
TEST_SCALAR_REVERSE("Multiplication scalar reverse", *);


TEST_SCALAR_ASSIGN("Addition scalar assign", +);
TEST_SCALAR_ASSIGN("Subtraction scalar assign", -);
TEST_SCALAR_ASSIGN("Division scalar assign", /);
TEST_SCALAR_ASSIGN("Multiplication scalar assign", *);


//------------------------------------------------------------------------------
// Vector
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Matrix - Multiply vector (precede)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, OrdersPrecede, LayoutsAll, PackingsAll>,
												VectorCaseList<ScalarsFloatingAndComplex32, PackingsAll>>{})) {

	using MatLhs = typename TestType::Lhs::template Matrix<2, 3>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;
	const MatLhs a = {
		1, 2, 3,
		4, 5, 6
	};
	const VecRhs b = { 2, 3, 4 };
	const auto r = a * b;
	using Scalar = scalar_type_t<std::decay_t<decltype(r)>>;
	REQUIRE(r[0] == static_cast<Scalar>(20));
	REQUIRE(r[1] == static_cast<Scalar>(47));
}

TEMPLATE_LIST_TEST_CASE("Matrix - Multiply vector (follow)", "[Matrix]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatingAndComplex32, PackingsAll>,
												MatrixCaseList<ScalarsFloatingAndComplex32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {

	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using MatRhs = typename TestType::Rhs::template Matrix<3, 2>;
	const VecLhs a = { 2, 3, 4 };
	const MatRhs b = {
		1, 4,
		2, 5,
		3, 6
	};
	const auto r = a * b;
	using Scalar = scalar_type_t<std::decay_t<decltype(r)>>;
	REQUIRE(r[0] == static_cast<Scalar>(20));
	REQUIRE(r[1] == static_cast<Scalar>(47));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiply vector - homogeneous augmentation (precede)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, OrdersPrecede, LayoutsAll, PackingsAll>,
												VectorCaseList<ScalarsFloatingAndComplex32, PackingsAll>>{})) {
	using MatLhs = typename TestType::Lhs::template Matrix<4, 4>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;
	const MatLhs a = {
		1, 0, 0, 0,
		0, 0, -1, 0,
		0, 1, 0, 0,
		0, 0, 0, 2
	};
	const VecRhs b = { 2, 4, 6 };
	const auto r = a * b;
	static_assert(dimension_v<std::decay_t<decltype(r)>> == 3);
	using Scalar = scalar_type_t<std::decay_t<decltype(r)>>;
	REQUIRE(r[0] == static_cast<Scalar>(1));
	REQUIRE(r[1] == static_cast<Scalar>(-3));
	REQUIRE(r[2] == static_cast<Scalar>(2));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiply vector - homogeneous augmentation (follow)", "[Matrix]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatingAndComplex32, PackingsAll>,
												MatrixCaseList<ScalarsFloatingAndComplex32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using MatRhs = typename TestType::Rhs::template Matrix<4, 4>;
	const VecLhs a = { 2, 4, 6 };
	const MatRhs b = {
		1, 0, 0, 0,
		0, 0, 1, 0,
		0, -1, 0, 0,
		0, 0, 0, 2
	};
	const auto r = a * b;
	static_assert(dimension_v<std::decay_t<decltype(r)>> == 3);
	using Scalar = scalar_type_t<std::decay_t<decltype(r)>>;
	REQUIRE(r[0] == static_cast<Scalar>(1));
	REQUIRE(r[1] == static_cast<Scalar>(-3));
	REQUIRE(r[2] == static_cast<Scalar>(2));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiply vector - affine augmentation (precede)", "[Matrix]",
						decltype(BinaryCaseList<MatrixCaseList<ScalarsFloatingAndComplex32, OrdersPrecede, LayoutsAll, PackingsAll>,
												VectorCaseList<ScalarsFloatingAndComplex32, PackingsAll>>{})) {
	using MatLhs = typename TestType::Lhs::template Matrix<3, 4>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;
	const MatLhs a = {
		1, 0, 0, 3,
		0, 0, -1, 4,
		0, 1, 0, 5
	};
	const VecRhs b = { 1, 2, 3 };
	const auto r = a * b;
	static_assert(dimension_v<std::decay_t<decltype(r)>> == 3);
	using Scalar = scalar_type_t<std::decay_t<decltype(r)>>;
	REQUIRE(r[0] == static_cast<Scalar>(4));
	REQUIRE(r[1] == static_cast<Scalar>(1));
	REQUIRE(r[2] == static_cast<Scalar>(7));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiply vector - affine augmentation (follow)", "[Matrix]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatingAndComplex32, PackingsAll>,
												MatrixCaseList<ScalarsFloatingAndComplex32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using MatRhs = typename TestType::Rhs::template Matrix<4, 3>;
	const VecLhs a = { 1, 2, 3 };
	const MatRhs b = {
		1, 0, 0,
		0, 0, 1,
		0, -1, 0,
		3, 4, 5
	};
	const auto r = a * b;
	static_assert(dimension_v<std::decay_t<decltype(r)>> == 3);
	using Scalar = scalar_type_t<std::decay_t<decltype(r)>>;
	REQUIRE(r[0] == static_cast<Scalar>(4));
	REQUIRE(r[1] == static_cast<Scalar>(1));
	REQUIRE(r[2] == static_cast<Scalar>(7));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiply-assign vector (only follow)", "[Matrix]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloat32, PackingsAll>,
												MatrixCaseList<ScalarsFloat32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using MatRhs = typename TestType::Rhs::template Matrix<3, 3>;
	const VecLhs a = { 1, 2, 3 };
	const MatRhs b = {
		1, 0, 0,
		0, 0, 1,
		0, -1, 0
	};
	auto copy = a;
	REQUIRE(&(copy *= b) == &copy);
	using Scalar = scalar_type_t<std::decay_t<decltype(copy)>>;
	REQUIRE(copy[0] == static_cast<Scalar>(1));
	REQUIRE(copy[1] == static_cast<Scalar>(-3));
	REQUIRE(copy[2] == static_cast<Scalar>(2));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiply-assign vector - homogeneous augmentation (only follow)", "[Matrix]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloat32, PackingsAll>,
												MatrixCaseList<ScalarsFloat32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using MatRhs = typename TestType::Rhs::template Matrix<4, 4>;
	const VecLhs a = { 2, 4, 6 };
	const MatRhs b = {
		1, 0, 0, 0,
		0, 0, 1, 0,
		0, -1, 0, 0,
		0, 0, 0, 2
	};
	auto copy = a;
	REQUIRE(&(copy *= b) == &copy);
	using Scalar = scalar_type_t<std::decay_t<decltype(copy)>>;
	REQUIRE(copy[0] == static_cast<Scalar>(1));
	REQUIRE(copy[1] == static_cast<Scalar>(-3));
	REQUIRE(copy[2] == static_cast<Scalar>(2));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Multiply-assign vector - affine augmentation (only follow)", "[Matrix]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloat32, PackingsAll>,
												MatrixCaseList<ScalarsFloat32, OrdersFollow, LayoutsAll, PackingsAll>>{})) {
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using MatRhs = typename TestType::Rhs::template Matrix<4, 3>;
	const VecLhs a = { 1, 2, 3 };
	const MatRhs b = {
		1, 0, 0,
		0, 0, 1,
		0, -1, 0,
		3, 4, 5
	};
	auto copy = a;
	REQUIRE(&(copy *= b) == &copy);
	using Scalar = scalar_type_t<std::decay_t<decltype(copy)>>;
	REQUIRE(copy[0] == static_cast<Scalar>(4));
	REQUIRE(copy[1] == static_cast<Scalar>(1));
	REQUIRE(copy[2] == static_cast<Scalar>(7));
}


//------------------------------------------------------------------------------
// Plus & minus
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Matrix - Plus", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {

	using Mat = typename TestType::template Matrix<2, 3>;
	const Mat value = {
		1, 2, 3,
		4, 5, 6
	};
	const auto r = +value;

	REQUIRE(r(0, 0) == value(0, 0));
	REQUIRE(r(0, 1) == value(0, 1));
	REQUIRE(r(0, 2) == value(0, 2));
	REQUIRE(r(1, 0) == value(1, 0));
	REQUIRE(r(1, 1) == value(1, 1));
	REQUIRE(r(1, 2) == value(1, 2));
}


TEMPLATE_LIST_TEST_CASE("Matrix - Minus", "[Matrix]",
						decltype(MatrixCaseList<ScalarsAll, OrdersAll, LayoutsAll, PackingsAll>{})) {

	using Mat = typename TestType::template Matrix<2, 3>;
	const Mat value = {
		1, 2, 3,
		4, 5, 6
	};
	const auto r = -value;

	REQUIRE(r(0, 0) == -value(0, 0));
	REQUIRE(r(0, 1) == -value(0, 1));
	REQUIRE(r(0, 2) == -value(0, 2));
	REQUIRE(r(1, 0) == -value(1, 0));
	REQUIRE(r(1, 1) == -value(1, 1));
	REQUIRE(r(1, 2) == -value(1, 2));
}