// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Quaternion/Arithmetic.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


//------------------------------------------------------------------------------
// Quat x quat
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Quaternion - Multiplication (quat x quat)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using QuatRhs = typename TestType::Rhs::Quat;

	const QuatLhs lhs(1, 2, 3, 4);
	const QuatRhs rhs(5, 6, 7, 8);
	const auto result = lhs * rhs;
	REQUIRE(double(result.s) == Catch::Approx(-60));
	REQUIRE(double(result.i) == Catch::Approx(12));
	REQUIRE(double(result.j) == Catch::Approx(30));
	REQUIRE(double(result.k) == Catch::Approx(24));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Addition (quat x quat)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using QuatRhs = typename TestType::Rhs::Quat;

	const QuatLhs lhs(1, 2, 3, 4);
	const QuatRhs rhs(5, 6, 7, 8);
	const auto result = lhs + rhs;
	REQUIRE(double(result.s) == Catch::Approx(6));
	REQUIRE(double(result.i) == Catch::Approx(8));
	REQUIRE(double(result.j) == Catch::Approx(10));
	REQUIRE(double(result.k) == Catch::Approx(12));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Subtraction (quat x quat)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using QuatRhs = typename TestType::Rhs::Quat;

	const QuatLhs lhs(1, 2, 3, 4);
	const QuatRhs rhs(5, 3, 7, 3);
	const auto result = lhs - rhs;
	REQUIRE(double(result.s) == Catch::Approx(-4));
	REQUIRE(double(result.i) == Catch::Approx(-1));
	REQUIRE(double(result.j) == Catch::Approx(-4));
	REQUIRE(double(result.k) == Catch::Approx(1));
}


//------------------------------------------------------------------------------
// Quat x quat assign
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Quaternion - Multiplication assign (quat x quat)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using QuatRhs = typename TestType::Rhs::Quat;

	const QuatLhs lhs(1, 2, 3, 4);
	const QuatRhs rhs(5, 6, 7, 8);
	auto copy = lhs;
	REQUIRE(&(copy *= rhs) == &copy);
	REQUIRE(double(copy.s) == Catch::Approx(-60));
	REQUIRE(double(copy.i) == Catch::Approx(12));
	REQUIRE(double(copy.j) == Catch::Approx(30));
	REQUIRE(double(copy.k) == Catch::Approx(24));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Addition assign (quat x quat)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using QuatRhs = typename TestType::Rhs::Quat;

	const QuatLhs lhs(1, 2, 3, 4);
	const QuatRhs rhs(5, 6, 7, 8);
	auto copy = lhs;
	REQUIRE(&(copy += rhs) == &copy);
	REQUIRE(double(copy.s) == Catch::Approx(6));
	REQUIRE(double(copy.i) == Catch::Approx(8));
	REQUIRE(double(copy.j) == Catch::Approx(10));
	REQUIRE(double(copy.k) == Catch::Approx(12));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Subtraction assign (quat x quat)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using QuatRhs = typename TestType::Rhs::Quat;

	const QuatLhs lhs(1, 2, 3, 4);
	const QuatRhs rhs(5, 3, 7, 3);
	auto copy = lhs;
	REQUIRE(&(copy -= rhs) == &copy);
	REQUIRE(double(copy.s) == Catch::Approx(-4));
	REQUIRE(double(copy.i) == Catch::Approx(-1));
	REQUIRE(double(copy.j) == Catch::Approx(-4));
	REQUIRE(double(copy.k) == Catch::Approx(1));
}


//------------------------------------------------------------------------------
// Quat x scalar
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Quaternion - Multiplication (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(1, 2, 3, 4);
	const ScalarRhs rhs(2);
	const auto result = lhs * rhs;
	REQUIRE(double(result.s) == Catch::Approx(2));
	REQUIRE(double(result.i) == Catch::Approx(4));
	REQUIRE(double(result.j) == Catch::Approx(6));
	REQUIRE(double(result.k) == Catch::Approx(8));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Division (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(2, 4, 6, 8);
	const ScalarRhs rhs(2);
	const auto result = lhs / rhs;
	REQUIRE(double(result.s) == Catch::Approx(1));
	REQUIRE(double(result.i) == Catch::Approx(2));
	REQUIRE(double(result.j) == Catch::Approx(3));
	REQUIRE(double(result.k) == Catch::Approx(4));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Addition (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(1, 2, 3, 4);
	const ScalarRhs rhs(5);
	const auto result = lhs + rhs;
	REQUIRE(double(result.s) == Catch::Approx(6));
	REQUIRE(double(result.i) == Catch::Approx(2));
	REQUIRE(double(result.j) == Catch::Approx(3));
	REQUIRE(double(result.k) == Catch::Approx(4));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Subtraction (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(1, 2, 3, 4);
	const ScalarRhs rhs(5);
	const auto result = lhs - rhs;
	REQUIRE(double(result.s) == Catch::Approx(-4));
	REQUIRE(double(result.i) == Catch::Approx(2));
	REQUIRE(double(result.j) == Catch::Approx(3));
	REQUIRE(double(result.k) == Catch::Approx(4));
}


//------------------------------------------------------------------------------
// Scalar x quat
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Quaternion - Multiplication (scalar x quat)", "[Quaternion]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using QuatRhs = typename TestType::Rhs::Quat;

	const ScalarLhs lhs(2);
	const QuatRhs rhs(1, 2, 3, 4);
	const auto result = lhs * rhs;
	REQUIRE(double(result.s) == Catch::Approx(2));
	REQUIRE(double(result.i) == Catch::Approx(4));
	REQUIRE(double(result.j) == Catch::Approx(6));
	REQUIRE(double(result.k) == Catch::Approx(8));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Addition (scalar x quat)", "[Quaternion]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using QuatRhs = typename TestType::Rhs::Quat;

	const ScalarLhs lhs(5);
	const QuatRhs rhs(1, 2, 3, 4);
	const auto result = lhs + rhs;
	REQUIRE(double(result.s) == Catch::Approx(6));
	REQUIRE(double(result.i) == Catch::Approx(2));
	REQUIRE(double(result.j) == Catch::Approx(3));
	REQUIRE(double(result.k) == Catch::Approx(4));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Subtraction (scalar x quat)", "[Quaternion]",
						decltype(BinaryCaseList<ScalarCaseList<ScalarsFloating>,
												QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>>{})) {
	using ScalarLhs = typename TestType::Lhs::Scalar;
	using QuatRhs = typename TestType::Rhs::Quat;

	const ScalarLhs lhs(5);
	const QuatRhs rhs(1, 2, 3, 4);
	const auto result = lhs - rhs;
	REQUIRE(double(result.s) == Catch::Approx(4));
	REQUIRE(double(result.i) == Catch::Approx(-2));
	REQUIRE(double(result.j) == Catch::Approx(-3));
	REQUIRE(double(result.k) == Catch::Approx(-4));
}

//------------------------------------------------------------------------------
// Quat x scalar assign
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Quaternion - Multiplication assign (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(1, 2, 3, 4);
	const ScalarRhs rhs(2);
	auto copy = lhs;
	REQUIRE(&(copy *= rhs) == &copy);
	REQUIRE(double(copy.s) == Catch::Approx(2));
	REQUIRE(double(copy.i) == Catch::Approx(4));
	REQUIRE(double(copy.j) == Catch::Approx(6));
	REQUIRE(double(copy.k) == Catch::Approx(8));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Division assign (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(2, 4, 6, 8);
	const ScalarRhs rhs(2);
	auto copy = lhs;
	REQUIRE(&(copy /= rhs) == &copy);
	REQUIRE(double(copy.s) == Catch::Approx(1));
	REQUIRE(double(copy.i) == Catch::Approx(2));
	REQUIRE(double(copy.j) == Catch::Approx(3));
	REQUIRE(double(copy.k) == Catch::Approx(4));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Addition assign (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(1, 2, 3, 4);
	const ScalarRhs rhs(5);
	auto copy = lhs;
	REQUIRE(&(copy += rhs) == &copy);
	REQUIRE(double(copy.s) == Catch::Approx(6));
	REQUIRE(double(copy.i) == Catch::Approx(2));
	REQUIRE(double(copy.j) == Catch::Approx(3));
	REQUIRE(double(copy.k) == Catch::Approx(4));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Subtraction assign (quat x scalar)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												ScalarCaseList<ScalarsFloating>>{})) {
	using QuatLhs = typename TestType::Lhs::Quat;
	using ScalarRhs = typename TestType::Rhs::Scalar;

	const QuatLhs lhs(1, 2, 3, 4);
	const ScalarRhs rhs(5);
	auto copy = lhs;
	REQUIRE(&(copy -= rhs) == &copy);
	REQUIRE(double(copy.s) == Catch::Approx(-4));
	REQUIRE(double(copy.i) == Catch::Approx(2));
	REQUIRE(double(copy.j) == Catch::Approx(3));
	REQUIRE(double(copy.k) == Catch::Approx(4));
}