// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Vector/Arithmetic.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace mathter;
using namespace test_util;


#define TEST_ARITHMETIC_VEC_x_VEC(OP, SCALARS)                                                  \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: vec " #OP " vec", "[Vector]",                 \
							decltype(BinaryCaseList<VectorCaseList<SCALARS, PackingsAll>,       \
													VectorCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                  \
                                                                                                \
		using VecLhs = typename TestType::Lhs::template Vector<Dim>;                            \
		using VecRhs = typename TestType::Rhs::template Vector<Dim>;                            \
                                                                                                \
		const VecLhs lhs{ 3, 6, 9 };                                                            \
		const VecRhs rhs{ 2, 3, 4 };                                                            \
		const auto result = lhs OP rhs;                                                         \
                                                                                                \
		for (int i = 0; i < Dim; ++i) {                                                         \
			REQUIRE(result[i] == (lhs[i] OP rhs[i]));                                           \
		}                                                                                       \
	}


#define TEST_ARITHMETIC_VEC_x_SCALAR(OP, SCALARS)                                         \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: vec " #OP " scalar", "[Vector]",        \
							decltype(BinaryCaseList<VectorCaseList<SCALARS, PackingsAll>, \
													ScalarCaseList<SCALARS>>{})) {        \
		constexpr int Dim = 3;                                                            \
                                                                                          \
		using VecLhs = typename TestType::Lhs::template Vector<Dim>;                      \
		using ScalarRhs = typename TestType::Rhs::Scalar;                                 \
                                                                                          \
		const VecLhs lhs{ 3, 6, 9 };                                                      \
		const ScalarRhs rhs(7);                                                           \
		const auto result = lhs OP rhs;                                                   \
                                                                                          \
		for (int i = 0; i < Dim; ++i) {                                                   \
			REQUIRE(result[i] == (lhs[i] OP rhs));                                        \
		}                                                                                 \
	}


#define TEST_ARITHMETIC_SCALAR_x_VEC(OP, SCALARS)                                               \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: scalar " #OP " vec", "[Vector]",              \
							decltype(BinaryCaseList<ScalarCaseList<SCALARS>,                    \
													VectorCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                  \
                                                                                                \
		using ScalarLhs = typename TestType::Lhs::Scalar;                                       \
		using VecRhs = typename TestType::Rhs::template Vector<Dim>;                            \
                                                                                                \
		const ScalarLhs lhs(7);                                                                 \
		const VecRhs rhs{ 3, 6, 9 };                                                            \
		const auto result = lhs OP rhs;                                                         \
                                                                                                \
		for (int i = 0; i < Dim; ++i) {                                                         \
			REQUIRE(result[i] == (lhs OP rhs[i]));                                              \
		}                                                                                       \
	}


#define TEST_ARITHMETIC_VEC_x_SWIZZLE(OP, SCALARS)                                               \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: vec " #OP " swizzle", "[Vector]",              \
							decltype(BinaryCaseList<VectorCaseList<SCALARS, PackingsAll>,        \
													SwizzleCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                   \
                                                                                                 \
		using VecLhs = typename TestType::Lhs::template Vector<Dim>;                             \
		using SwizRhs = typename TestType::Rhs::template Swizzle<Dim, 2, 0, 1>;                  \
                                                                                                 \
		const VecLhs lhs{ 3, 6, 9 };                                                             \
		const SwizRhs rhs{ 2, 3, 4 };                                                            \
		const auto result = lhs OP rhs;                                                          \
                                                                                                 \
		for (int i = 0; i < Dim; ++i) {                                                          \
			REQUIRE(result[i] == (lhs[i] OP rhs[i]));                                            \
		}                                                                                        \
	}


#define TEST_ARITHMETIC_SWIZZLE_x_VEC(OP, SCALARS)                                              \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: swizzle " #OP " vec", "[Vector]",             \
							decltype(BinaryCaseList<SwizzleCaseList<SCALARS, PackingsAll>,      \
													VectorCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                  \
                                                                                                \
		using SwizLhs = typename TestType::Lhs::template Swizzle<Dim, 2, 0, 1>;                 \
		using VecRhs = typename TestType::Rhs::template Vector<Dim>;                            \
                                                                                                \
		const SwizLhs lhs{ 2, 3, 4 };                                                           \
		const VecRhs rhs{ 3, 6, 9 };                                                            \
		const auto result = lhs OP rhs;                                                         \
                                                                                                \
		for (int i = 0; i < Dim; ++i) {                                                         \
			REQUIRE(result[i] == (lhs[i] OP rhs[i]));                                           \
		}                                                                                       \
	}


#define TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(OP, SCALARS)                                           \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: swizzle " #OP " swizzle", "[Vector]",          \
							decltype(BinaryCaseList<SwizzleCaseList<SCALARS, PackingsAll>,       \
													SwizzleCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                   \
                                                                                                 \
		using SwizLhs = typename TestType::Lhs::template Swizzle<Dim, 2, 0, 1>;                  \
		using SwizRhs = typename TestType::Rhs::template Swizzle<Dim, 1, 2, 0>;                  \
                                                                                                 \
		const SwizLhs lhs{ 2, 3, 4 };                                                            \
		const SwizRhs rhs{ 3, 6, 9 };                                                            \
		const auto result = lhs OP rhs;                                                          \
                                                                                                 \
		for (int i = 0; i < Dim; ++i) {                                                          \
			REQUIRE(result[i] == (lhs[i] OP rhs[i]));                                            \
		}                                                                                        \
	}


#define TEST_ARITHMETIC_SWIZZLE_x_SCALAR(OP, SCALARS)                                      \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: swizzle " #OP " scalar", "[Vector]",     \
							decltype(BinaryCaseList<SwizzleCaseList<SCALARS, PackingsAll>, \
													ScalarCaseList<SCALARS>>{})) {         \
		constexpr int Dim = 3;                                                             \
                                                                                           \
		using SwizLhs = typename TestType::Lhs::template Swizzle<Dim, 2, 0, 1>;            \
		using ScalarRhs = typename TestType::Rhs::Scalar;                                  \
                                                                                           \
		const SwizLhs lhs{ 3, 6, 9 };                                                      \
		const ScalarRhs rhs(7);                                                            \
		const auto result = lhs OP rhs;                                                    \
                                                                                           \
		for (int i = 0; i < Dim; ++i) {                                                    \
			REQUIRE(result[i] == (lhs[i] OP rhs));                                         \
		}                                                                                  \
	}


#define TEST_ARITHMETIC_SCALAR_x_SWIZZLE(OP, SCALARS)                                            \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: scalar " #OP " swizzle", "[Vector]",           \
							decltype(BinaryCaseList<ScalarCaseList<SCALARS>,                     \
													SwizzleCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                   \
                                                                                                 \
		using ScalarLhs = typename TestType::Lhs::Scalar;                                        \
		using SwizRhs = typename TestType::Rhs::template Swizzle<Dim, 2, 0, 1>;                  \
                                                                                                 \
		const ScalarLhs lhs(7);                                                                  \
		const SwizRhs rhs{ 3, 6, 9 };                                                            \
		const auto result = lhs OP rhs;                                                          \
                                                                                                 \
		for (int i = 0; i < Dim; ++i) {                                                          \
			REQUIRE(result[i] == (lhs OP rhs[i]));                                               \
		}                                                                                        \
	}


TEST_ARITHMETIC_VEC_x_VEC(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_VEC(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_VEC(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_VEC(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_VEC_x_VEC(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_VEC(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_VEC(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_VEC(-, ScalarsFloatingAndComplex32);


TEST_ARITHMETIC_VEC_x_SCALAR(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SCALAR(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SCALAR(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SCALAR(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_VEC_x_SCALAR(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_SCALAR(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_SCALAR(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_SCALAR(-, ScalarsFloatingAndComplex32);


TEST_ARITHMETIC_SCALAR_x_VEC(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SCALAR_x_VEC(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SCALAR_x_VEC(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SCALAR_x_VEC(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_SCALAR_x_VEC(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SCALAR_x_VEC(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SCALAR_x_VEC(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SCALAR_x_VEC(-, ScalarsFloatingAndComplex32);


TEST_ARITHMETIC_VEC_x_SWIZZLE(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SWIZZLE(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SWIZZLE(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SWIZZLE(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_VEC_x_SWIZZLE(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_SWIZZLE(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_SWIZZLE(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_VEC_x_SWIZZLE(-, ScalarsFloatingAndComplex32);


TEST_ARITHMETIC_SWIZZLE_x_VEC(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_VEC(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_VEC(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_VEC(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_SWIZZLE_x_VEC(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_VEC(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_VEC(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_VEC(-, ScalarsFloatingAndComplex32);


TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE(-, ScalarsFloatingAndComplex32);


TEST_ARITHMETIC_SWIZZLE_x_SCALAR(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_SWIZZLE_x_SCALAR(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR(-, ScalarsFloatingAndComplex32);


TEST_ARITHMETIC_SCALAR_x_SWIZZLE(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SCALAR_x_SWIZZLE(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SCALAR_x_SWIZZLE(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SCALAR_x_SWIZZLE(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_SCALAR_x_SWIZZLE(*, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SCALAR_x_SWIZZLE(/, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SCALAR_x_SWIZZLE(+, ScalarsFloatingAndComplex32);
TEST_ARITHMETIC_SCALAR_x_SWIZZLE(-, ScalarsFloatingAndComplex32);


#define TEST_ARITHMETIC_VEC_x_VEC_ASSIGN(OP, SCALARS)                                           \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: vec " #OP "= vec", "[Vector]",                \
							decltype(BinaryCaseList<VectorCaseList<SCALARS, PackingsAll>,       \
													VectorCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                  \
                                                                                                \
		using VecLhs = typename TestType::Lhs::template Vector<Dim>;                            \
		using VecRhs = typename TestType::Rhs::template Vector<Dim>;                            \
                                                                                                \
		const VecLhs lhs{ 3, 6, 9 };                                                            \
		const VecRhs rhs{ 2, 3, 4 };                                                            \
		auto result = lhs;                                                                      \
		result OP## = rhs;                                                                      \
                                                                                                \
		for (int i = 0; i < Dim; ++i) {                                                         \
			REQUIRE(result[i] == scalar_type_t<VecLhs>(lhs[i] OP rhs[i]));                      \
		}                                                                                       \
	}


#define TEST_ARITHMETIC_VEC_x_SCALAR_ASSIGN(OP, SCALARS)                                  \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: vec " #OP "= scalar", "[Vector]",       \
							decltype(BinaryCaseList<VectorCaseList<SCALARS, PackingsAll>, \
													ScalarCaseList<SCALARS>>{})) {        \
		constexpr int Dim = 3;                                                            \
                                                                                          \
		using VecLhs = typename TestType::Lhs::template Vector<Dim>;                      \
		using ScalarRhs = typename TestType::Rhs::Scalar;                                 \
                                                                                          \
		const VecLhs lhs{ 3, 6, 9 };                                                      \
		const ScalarRhs rhs(3);                                                           \
		auto result = lhs;                                                                \
		result OP## = rhs;                                                                \
                                                                                          \
		for (int i = 0; i < Dim; ++i) {                                                   \
			REQUIRE(result[i] == scalar_type_t<VecLhs>(lhs[i] OP rhs));                   \
		}                                                                                 \
	}


#define TEST_ARITHMETIC_VEC_x_SWIZZLE_ASSIGN(OP, SCALARS)                                        \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: vec " #OP "= swizzle", "[Vector]",             \
							decltype(BinaryCaseList<VectorCaseList<SCALARS, PackingsAll>,        \
													SwizzleCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                   \
                                                                                                 \
		using VecLhs = typename TestType::Lhs::template Vector<Dim>;                             \
		using SwizRhs = typename TestType::Rhs::template Swizzle<Dim, 2, 0, 1>;                  \
                                                                                                 \
		const VecLhs lhs{ 3, 6, 9 };                                                             \
		const SwizRhs rhs{ 2, 3, 4 };                                                            \
		auto result = lhs;                                                                       \
		result OP## = rhs;                                                                       \
                                                                                                 \
		for (int i = 0; i < Dim; ++i) {                                                          \
			REQUIRE(result[i] == scalar_type_t<VecLhs>(lhs[i] OP rhs[i]));                       \
		}                                                                                        \
	}


#define TEST_ARITHMETIC_SWIZZLE_x_VEC_ASSIGN(OP, SCALARS)                                       \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: swizzle " #OP "= vec", "[Vector]",            \
							decltype(BinaryCaseList<SwizzleCaseList<SCALARS, PackingsAll>,      \
													VectorCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                  \
                                                                                                \
		using SwizLhs = typename TestType::Lhs::template Swizzle<Dim, 1, 0, 2>;                 \
		using VecRhs = typename TestType::Rhs::template Vector<Dim>;                            \
                                                                                                \
		const SwizLhs lhs{ 3, 6, 9 };                                                           \
		const VecRhs rhs{ 2, 3, 4 };                                                            \
		auto result = lhs;                                                                      \
		result OP## = rhs;                                                                      \
                                                                                                \
		for (int i = 0; i < Dim; ++i) {                                                         \
			REQUIRE(result[i] == scalar_type_t<SwizLhs>(lhs[i] OP rhs[i]));                     \
		}                                                                                       \
	}


#define TEST_ARITHMETIC_SWIZZLE_x_SCALAR_ASSIGN(OP, SCALARS)                               \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: swizzle " #OP "= scalar", "[Vector]",    \
							decltype(BinaryCaseList<SwizzleCaseList<SCALARS, PackingsAll>, \
													ScalarCaseList<SCALARS>>{})) {         \
		constexpr int Dim = 3;                                                             \
                                                                                           \
		using SwizLhs = typename TestType::Lhs::template Swizzle<Dim, 1, 0, 2>;            \
		using ScalarRhs = typename TestType::Rhs::Scalar;                                  \
                                                                                           \
		const SwizLhs lhs{ 3, 6, 9 };                                                      \
		const ScalarRhs rhs(3);                                                            \
		auto result = lhs;                                                                 \
		result OP## = rhs;                                                                 \
                                                                                           \
		for (int i = 0; i < Dim; ++i) {                                                    \
			REQUIRE(result[i] == scalar_type_t<SwizLhs>(lhs[i] OP rhs));                   \
		}                                                                                  \
	}


#define TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE_ASSIGN(OP, SCALARS)                                    \
	TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: swizzle " #OP "= swizzle", "[Vector]",         \
							decltype(BinaryCaseList<SwizzleCaseList<SCALARS, PackingsAll>,       \
													SwizzleCaseList<SCALARS, PackingsAll>>{})) { \
		constexpr int Dim = 3;                                                                   \
                                                                                                 \
		using SwizLhs = typename TestType::Lhs::template Swizzle<Dim, 1, 0, 2>;                  \
		using SwizRhs = typename TestType::Rhs::template Swizzle<Dim, 2, 0, 1>;                  \
                                                                                                 \
		const SwizLhs lhs{ 3, 6, 9 };                                                            \
		const SwizRhs rhs{ 2, 3, 4 };                                                            \
		auto result = lhs;                                                                       \
		result OP## = rhs;                                                                       \
                                                                                                 \
		for (int i = 0; i < Dim; ++i) {                                                          \
			REQUIRE(result[i] == scalar_type_t<SwizLhs>(lhs[i] OP rhs[i]));                      \
		}                                                                                        \
	}


TEST_ARITHMETIC_VEC_x_VEC_ASSIGN(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_VEC_ASSIGN(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_VEC_ASSIGN(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_VEC_ASSIGN(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_VEC_x_SCALAR_ASSIGN(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SCALAR_ASSIGN(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SCALAR_ASSIGN(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SCALAR_ASSIGN(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_VEC_x_SWIZZLE_ASSIGN(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SWIZZLE_ASSIGN(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SWIZZLE_ASSIGN(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_VEC_x_SWIZZLE_ASSIGN(-, ScalarsFloatAndInt32);


TEST_ARITHMETIC_SWIZZLE_x_VEC_ASSIGN(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_VEC_ASSIGN(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_VEC_ASSIGN(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_VEC_ASSIGN(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_SWIZZLE_x_SCALAR_ASSIGN(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR_ASSIGN(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR_ASSIGN(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SCALAR_ASSIGN(-, ScalarsFloatAndInt32);

TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE_ASSIGN(*, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE_ASSIGN(/, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE_ASSIGN(+, ScalarsFloatAndInt32);
TEST_ARITHMETIC_SWIZZLE_x_SWIZZLE_ASSIGN(-, ScalarsFloatAndInt32);


TEMPLATE_LIST_TEST_CASE("Vector - Arithmetic: FMA vector", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {
	constexpr int Dim = 3;

	using Vec = typename TestType::template Vector<Dim>;

	const Vec a{ 3, 6, 9 };
	const Vec b{ 2, 3, 4 };
	const Vec c{ 2, 1, 3 };
	const auto result = MultiplyAdd(a, b, c);

	for (int i = 0; i < Dim; ++i) {
		const auto expected = a[i] * b[i] + c[i];
		REQUIRE(std::real(result[i]) == Catch::Approx(std::real(expected)));
		REQUIRE(std::imag(result[i]) == Catch::Approx(std::imag(expected)));
	}
}