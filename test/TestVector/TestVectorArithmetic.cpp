// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Vector.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace mathter;
using Catch::Approx;


using TypeListAll = TestTypeList<TypesAll, PackedAll>;
using TypeListFloating = TestTypeList<TypesFloating, PackedAll>;


//------------------------------------------------------------------------------
// Vector-vector
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector add", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(5, 7, 9);

		REQUIRE(a3 + b3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 b5(4, 5, 6, 7, 8);
		Vec5 c5(5, 7, 9, 11, 13);

		REQUIRE(a5 + b5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector sub", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(-3, -3, -3);

		REQUIRE(a3 - b3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 b5(4, 5, 6, 7, 8);
		Vec5 c5(-3, -3, -3, -3, -3);

		REQUIRE(a5 - b5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector multiply", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(4, 10, 18);

		REQUIRE(a3 * b3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 b5(4, 5, 6, 7, 8);
		Vec5 c5(4, 10, 18, 28, 40);

		REQUIRE(a5 * b5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector div", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(0.25f, 0.4f, 0.5f);

		REQUIRE(a3 / b3 == ApproxVec(c3));

		Vec5 a5(2, 6, 6, 12, 10);
		Vec5 b5(1, 2, 3, 4, 5);
		Vec5 c5(2, 3, 2, 3, 2);

		REQUIRE(a5 / b5 == ApproxVec(c5));
	}
}



TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector compound add", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(5, 7, 9);
		a3 += b3;
		REQUIRE(a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 b5(4, 5, 6, 7, 8);
		Vec5 c5(5, 7, 9, 11, 13);
		a5 += b5;
		REQUIRE(a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector compound sub", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(-3, -3, -3);
		a3 -= b3;
		REQUIRE(a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 b5(4, 5, 6, 7, 8);
		Vec5 c5(-3, -3, -3, -3, -3);
		a5 -= b5;
		REQUIRE(a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector compound multiply", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(4, 10, 18);
		a3 *= b3;
		REQUIRE(a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 b5(4, 5, 6, 7, 8);
		Vec5 c5(4, 10, 18, 28, 40);
		a5 *= b5;
		REQUIRE(a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector compound div", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(4, 5, 6);
		Vec3 c3(0.25f, 0.4f, 0.5f);
		a3 /= b3;
		REQUIRE(a3 == ApproxVec(c3));

		Vec5 a5(2, 6, 6, 12, 10);
		Vec5 b5(1, 2, 3, 4, 5);
		Vec5 c5(2, 3, 2, 3, 2);
		a5 /= b5;
		REQUIRE(a5 == ApproxVec(c5));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Vector-vector fma", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, 2, 3);
		Vec3 b3(2, 2, 2);
		Vec3 c3(3, 2, 1);
		Vec3 e3(5, 6, 7);
		const auto r3 = MultiplyAdd(a3, b3, c3);
		REQUIRE(e3 == ApproxVec(r3));

		Vec5 a5(0, 0, 1, 2, 3);
		Vec5 b5(0, 0, 2, 2, 2);
		Vec5 c5(0, 0, 3, 2, 1);
		Vec5 e5(0, 0, 5, 6, 7);
		const auto r5 = MultiplyAdd(a5, b5, c5);
		REQUIRE(e5 == ApproxVec(r5));
	}
}



//------------------------------------------------------------------------------
// Vector-scalar
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar add", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(5, 6, 7);

		REQUIRE(a3 + b == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(5, 6, 7, 8, 9);

		REQUIRE(a5 + b == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar sub", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(-3, -2, -1);

		REQUIRE(a3 - b == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(-3, -2, -1, 0, 1);

		REQUIRE(a5 - b == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar multiply", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(4, 8, 12);

		REQUIRE(a3 * b == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(4, 8, 12, 16, 20);

		REQUIRE(a5 * b == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar div", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(4, 8, 12);
		Vec3 c3(1, 2, 3);

		REQUIRE(a3 / b == c3);

		Vec5 a5(4, 8, 12, 16, 20);
		Vec5 c5(1, 2, 3, 4, 5);

		REQUIRE(a5 / b == c5);
	}
}



TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar compound add", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(5, 6, 7);
		a3 += b;
		REQUIRE(a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(5, 6, 7, 8, 9);
		a5 += b;
		REQUIRE(a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar compound sub", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(-3, -2, -1);
		a3 -= b;
		REQUIRE(a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(-3, -2, -1, 0, 1);
		a5 -= b;
		REQUIRE(a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar compound multiply", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(4, 8, 12);
		a3 *= b;
		REQUIRE(a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(4, 8, 12, 16, 20);
		a5 *= b;
		REQUIRE(a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar compound div", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(4, 8, 12);
		Vec3 c3(1, 2, 3);
		a3 /= b;
		REQUIRE(a3 == c3);

		Vec5 a5(4, 8, 12, 16, 20);
		Vec5 c5(1, 2, 3, 4, 5);
		a5 /= b;
		REQUIRE(a5 == c5);
	}
}



TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar reverse add", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(5, 6, 7);

		REQUIRE(b + a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(5, 6, 7, 8, 9);

		REQUIRE(b + a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar reverse sub", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(-3, -2, -1);

		REQUIRE(b - a3 == -c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(-3, -2, -1, 0, 1);

		REQUIRE(b - a5 == -c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar reverse multiply", "[Vector]", TypeListAll) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(1, 2, 3);
		Vec3 c3(4, 8, 12);

		REQUIRE(b * a3 == c3);

		Vec5 a5(1, 2, 3, 4, 5);
		Vec5 c5(4, 8, 12, 16, 20);

		REQUIRE(b * a5 == c5);
	}
}

TEMPLATE_LIST_TEST_CASE("Vector - Vector-scalar reverse div", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Type b = 4;

		Vec3 a3(4, 8, 12);
		Vec3 c3(1, 1.0 / 2.0, 1.0 / 3.0);

		REQUIRE(b / a3 == c3);

		Vec5 a5(4, 8, 12, 16, 20);
		Vec5 c5(1, 1.0f / 2.0, 1.0f / 3.0, 1.0f / 4.0, 1.0f / 5.0);

		REQUIRE(b / a5 == c5);
	}
}

//------------------------------------------------------------------------------
// Swizzle
//------------------------------------------------------------------------------

#define TEST_CASE_SWIZZLE_VECTOR_OP(NAME, OPERATOR)                                          \
	TEMPLATE_LIST_TEST_CASE("Vector - Swizzle-vector " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                          \
			using Vec3 = typename TestType::template Vector<3>;                              \
			using Vec5 = typename TestType::template Vector<5>;                              \
                                                                                             \
			Vec3 v1 = { 1, 2, 3 };                                                           \
			Vec3 v2 = { 1, 4, -2 };                                                          \
			Vec3 r = v1.xyz OPERATOR v2;                                                     \
			Vec3 e = v1 OPERATOR v2;                                                         \
			REQUIRE(r == e);                                                                 \
		}                                                                                    \
	}

#define TEST_CASE_VECTOR_SWIZZLE_OP(NAME, OPERATOR)                                          \
	TEMPLATE_LIST_TEST_CASE("Vector - Vector-swizzle " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                          \
			using Vec3 = typename TestType::template Vector<3>;                              \
			using Vec5 = typename TestType::template Vector<5>;                              \
                                                                                             \
			Vec3 v1 = { 1, 2, 3 };                                                           \
			Vec3 v2 = { 1, 4, -2 };                                                          \
			Vec3 r = v1 OPERATOR v2.xyz;                                                     \
			Vec3 e = v1 OPERATOR v2;                                                         \
			REQUIRE(r == e);                                                                 \
		}                                                                                    \
	}

#define TEST_CASE_VECTOR_SWIZZLE_COMPOUND_OP(NAME, OPERATOR)                                 \
	TEMPLATE_LIST_TEST_CASE("Vector - Vector-swizzle " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                          \
			using Vec3 = typename TestType::template Vector<3>;                              \
			using Vec5 = typename TestType::template Vector<5>;                              \
                                                                                             \
			Vec3 v1 = { 1, 2, 3 };                                                           \
			auto v1c = v1;                                                                   \
			Vec3 v2 = { 1, 4, -2 };                                                          \
			v1 OPERATOR v2.xyz;                                                              \
			v1c OPERATOR v2;                                                                 \
			REQUIRE(v1 == v1c);                                                              \
		}                                                                                    \
	}

#define TEST_CASE_SWIZZLE_VECTOR_COMPOUND_OP(NAME, OPERATOR)                                 \
	TEMPLATE_LIST_TEST_CASE("Vector - Swizzle-vector " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                          \
			using Vec3 = typename TestType::template Vector<3>;                              \
			using Vec5 = typename TestType::template Vector<5>;                              \
                                                                                             \
			Vec3 v1 = { 1, 2, 3 };                                                           \
			auto v1c = v1;                                                                   \
			Vec3 v2 = { 1, 4, -2 };                                                          \
			v1.xyz OPERATOR v2;                                                              \
			v1c OPERATOR v2;                                                                 \
			REQUIRE(v1 == v1c);                                                              \
		}                                                                                    \
	}

#define TEST_CASE_SWIZZLE_SWIZZLE_OP(NAME, OPERATOR)                                          \
	TEMPLATE_LIST_TEST_CASE("Vector - Swizzle-swizzle " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                           \
			using Vec3 = typename TestType::template Vector<3>;                               \
			using Vec5 = typename TestType::template Vector<5>;                               \
                                                                                              \
			const Vec3 v1 = { 1, 2, 3 };                                                      \
			const Vec3 v2 = { 1, 4, -2 };                                                     \
			const Vec3 r = v1.xyz OPERATOR v2.xyz;                                            \
			const Vec3 e = v1 OPERATOR v2;                                                    \
			REQUIRE(r == e);                                                                  \
		}                                                                                     \
	}

#define TEST_CASE_SWIZZLE_SWIZZLE_COMPOUND_OP(NAME, OPERATOR)                                 \
	TEMPLATE_LIST_TEST_CASE("Vector - Swizzle-swizzle " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                           \
			using Vec3 = typename TestType::template Vector<3>;                               \
			using Vec5 = typename TestType::template Vector<5>;                               \
                                                                                              \
			Vec3 v1 = { 1, 2, 3 };                                                            \
			auto v1c = v1;                                                                    \
			Vec3 v2 = { 1, 4, -2 };                                                           \
			v1.xyz OPERATOR v2.xyz;                                                           \
			v1c OPERATOR v2;                                                                  \
			REQUIRE(v1 == v1c);                                                               \
		}                                                                                     \
	}


#define TEST_CASE_SWIZZLE_SCALAR_COMPOUND_OP(NAME, OPERATOR)                                         \
	TEMPLATE_LIST_TEST_CASE("Vector - Swizzle-scalar compound" NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                                  \
			using Vec3 = typename TestType::template Vector<3>;                                      \
			using Vec5 = typename TestType::template Vector<5>;                                      \
			using Type = scalar_type_t<Vec3>;                                                        \
                                                                                                     \
			Vec3 v1 = { 1, 2, 3 };                                                                   \
			auto v1c = v1;                                                                           \
			Type b = 6;                                                                              \
			v1.xyz OPERATOR b;                                                                       \
			v1c OPERATOR b;                                                                          \
			REQUIRE(v1 == v1c);                                                                      \
		}                                                                                            \
	}

#define TEST_CASE_SWIZZLE_SCALAR_OP(NAME, OPERATOR)                                          \
	TEMPLATE_LIST_TEST_CASE("Vector - Swizzle-scalar " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                          \
			using Vec3 = typename TestType::template Vector<3>;                              \
			using Vec5 = typename TestType::template Vector<5>;                              \
			using Type = scalar_type_t<Vec3>;                                                \
                                                                                             \
			Vec3 v1 = { 1, 2, 3 };                                                           \
			Type b = 6;                                                                      \
			Vec3 r = v1.xyz OPERATOR b;                                                      \
			Vec3 e = v1 OPERATOR b;                                                          \
			REQUIRE(r == e);                                                                 \
		}                                                                                    \
	}

#define TEST_CASE_SCALAR_SWIZZLE_OP(NAME, OPERATOR)                                          \
	TEMPLATE_LIST_TEST_CASE("Vector - Scalar-swizzle " NAME, "[Vector]", TypeListFloating) { \
		SECTION(TestType::Name()) {                                                          \
			using Vec3 = typename TestType::template Vector<3>;                              \
			using Vec5 = typename TestType::template Vector<5>;                              \
			using Type = scalar_type_t<Vec3>;                                                \
                                                                                             \
			Vec3 v1 = { 1, 2, 3 };                                                           \
			Type b = 6;                                                                      \
			Vec3 r = b OPERATOR v1.xyz;                                                      \
			Vec3 e = b OPERATOR v1;                                                          \
			REQUIRE(r == e);                                                                 \
		}                                                                                    \
	}



TEST_CASE_SWIZZLE_VECTOR_OP("add", +)
TEST_CASE_SWIZZLE_VECTOR_OP("sub", -)
TEST_CASE_SWIZZLE_VECTOR_OP("mul", *)
TEST_CASE_SWIZZLE_VECTOR_OP("div", /)

TEST_CASE_VECTOR_SWIZZLE_OP("add", +)
TEST_CASE_VECTOR_SWIZZLE_OP("sub", -)
TEST_CASE_VECTOR_SWIZZLE_OP("mul", *)
TEST_CASE_VECTOR_SWIZZLE_OP("div", /)

TEST_CASE_VECTOR_SWIZZLE_COMPOUND_OP("compound add", +=)
TEST_CASE_VECTOR_SWIZZLE_COMPOUND_OP("compound sub", -=)
TEST_CASE_VECTOR_SWIZZLE_COMPOUND_OP("compound mul", *=)
TEST_CASE_VECTOR_SWIZZLE_COMPOUND_OP("compound div", /=)

TEST_CASE_SWIZZLE_VECTOR_COMPOUND_OP("compound add", +=)
TEST_CASE_SWIZZLE_VECTOR_COMPOUND_OP("compound sub", -=)
TEST_CASE_SWIZZLE_VECTOR_COMPOUND_OP("compound mul", *=)
TEST_CASE_SWIZZLE_VECTOR_COMPOUND_OP("compound div", /=)

TEST_CASE_SWIZZLE_SWIZZLE_OP("add", +)
TEST_CASE_SWIZZLE_SWIZZLE_OP("sub", -)
TEST_CASE_SWIZZLE_SWIZZLE_OP("mul", *)
TEST_CASE_SWIZZLE_SWIZZLE_OP("div", /)

TEST_CASE_SWIZZLE_SWIZZLE_COMPOUND_OP("compound add", +=)
TEST_CASE_SWIZZLE_SWIZZLE_COMPOUND_OP("compound sub", -=)
TEST_CASE_SWIZZLE_SWIZZLE_COMPOUND_OP("compound mul", *=)
TEST_CASE_SWIZZLE_SWIZZLE_COMPOUND_OP("compound div", /=)

TEST_CASE_SWIZZLE_SCALAR_OP("add", +)
TEST_CASE_SWIZZLE_SCALAR_OP("sub", -)
TEST_CASE_SWIZZLE_SCALAR_OP("mul", *)
TEST_CASE_SWIZZLE_SCALAR_OP("div", /)

TEST_CASE_SCALAR_SWIZZLE_OP("add", +)
TEST_CASE_SCALAR_SWIZZLE_OP("sub", -)
TEST_CASE_SCALAR_SWIZZLE_OP("mul", *)
TEST_CASE_SCALAR_SWIZZLE_OP("div", /)

TEST_CASE_SWIZZLE_SCALAR_COMPOUND_OP("compound add", +=)
TEST_CASE_SWIZZLE_SCALAR_COMPOUND_OP("compound sub", -=)
TEST_CASE_SWIZZLE_SCALAR_COMPOUND_OP("compound mul", *=)
TEST_CASE_SWIZZLE_SCALAR_COMPOUND_OP("compound div", /=)