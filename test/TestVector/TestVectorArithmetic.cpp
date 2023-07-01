// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../TestGenerators.hpp"

#include <Mathter/Common/Approx.hpp>
#include <Mathter/Vector.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace mathter;
using Catch::Approx;


//------------------------------------------------------------------------------
// Vector-vector
//------------------------------------------------------------------------------

TEST_CASE_VEC_VARIANT("Vector - Vector-vector add", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(5, 7, 9);

		REQUIRE(a3 + b3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> b5(4, 5, 6, 7, 8);
		VectorT<5> c5(5, 7, 9, 11, 13);

		REQUIRE(a5 + b5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-vector sub", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(-3, -3, -3);

		REQUIRE(a3 - b3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> b5(4, 5, 6, 7, 8);
		VectorT<5> c5(-3, -3, -3, -3, -3);

		REQUIRE(a5 - b5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-vector multiply", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(4, 10, 18);

		REQUIRE(a3 * b3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> b5(4, 5, 6, 7, 8);
		VectorT<5> c5(4, 10, 18, 28, 40);

		REQUIRE(a5 * b5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-vector div", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(0.25f, 0.4f, 0.5f);

		REQUIRE(a3 / b3 == ApproxVec(c3));

		VectorT<5> a5(2, 6, 6, 12, 10);
		VectorT<5> b5(1, 2, 3, 4, 5);
		VectorT<5> c5(2, 3, 2, 3, 2);

		REQUIRE(a5 / b5 == ApproxVec(c5));
	}
}



TEST_CASE_VEC_VARIANT("Vector - Vector-vector compound add", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(5, 7, 9);
		a3 += b3;
		REQUIRE(a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> b5(4, 5, 6, 7, 8);
		VectorT<5> c5(5, 7, 9, 11, 13);
		a5 += b5;
		REQUIRE(a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-vector compound sub", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(-3, -3, -3);
		a3 -= b3;
		REQUIRE(a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> b5(4, 5, 6, 7, 8);
		VectorT<5> c5(-3, -3, -3, -3, -3);
		a5 -= b5;
		REQUIRE(a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-vector compound multiply", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(4, 10, 18);
		a3 *= b3;
		REQUIRE(a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> b5(4, 5, 6, 7, 8);
		VectorT<5> c5(4, 10, 18, 28, 40);
		a5 *= b5;
		REQUIRE(a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-vector compound div", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a3(1, 2, 3);
		VectorT<3> b3(4, 5, 6);
		VectorT<3> c3(0.25f, 0.4f, 0.5f);
		a3 /= b3;
		REQUIRE(a3 == ApproxVec(c3));

		VectorT<5> a5(2, 6, 6, 12, 10);
		VectorT<5> b5(1, 2, 3, 4, 5);
		VectorT<5> c5(2, 3, 2, 3, 2);
		a5 /= b5;
		REQUIRE(a5 == ApproxVec(c5));
	}
}



//------------------------------------------------------------------------------
// Vector-scalar
//------------------------------------------------------------------------------

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar add", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(5, 6, 7);

		REQUIRE(a3 + b == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(5, 6, 7, 8, 9);

		REQUIRE(a5 + b == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar sub", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(-3, -2, -1);

		REQUIRE(a3 - b == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(-3, -2, -1, 0, 1);

		REQUIRE(a5 - b == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar multiply", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(4, 8, 12);

		REQUIRE(a3 * b == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(4, 8, 12, 16, 20);

		REQUIRE(a5 * b == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar div", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(4, 8, 12);
		VectorT<3> c3(1, 2, 3);

		REQUIRE(a3 / b == c3);

		VectorT<5> a5(4, 8, 12, 16, 20);
		VectorT<5> c5(1, 2, 3, 4, 5);

		REQUIRE(a5 / b == c5);
	}
}



TEST_CASE_VEC_VARIANT("Vector - Vector-scalar compound add", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(5, 6, 7);
		a3 += b;
		REQUIRE(a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(5, 6, 7, 8, 9);
		a5 += b;
		REQUIRE(a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar compound sub", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(-3, -2, -1);
		a3 -= b;
		REQUIRE(a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(-3, -2, -1, 0, 1);
		a5 -= b;
		REQUIRE(a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar compound multiply", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(4, 8, 12);
		a3 *= b;
		REQUIRE(a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(4, 8, 12, 16, 20);
		a5 *= b;
		REQUIRE(a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar compound div", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(4, 8, 12);
		VectorT<3> c3(1, 2, 3);
		a3 /= b;
		REQUIRE(a3 == c3);

		VectorT<5> a5(4, 8, 12, 16, 20);
		VectorT<5> c5(1, 2, 3, 4, 5);
		a5 /= b;
		REQUIRE(a5 == c5);
	}
}



TEST_CASE_VEC_VARIANT("Vector - Vector-scalar reverse add", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(5, 6, 7);

		REQUIRE(b + a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(5, 6, 7, 8, 9);

		REQUIRE(b + a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar reverse sub", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(-3, -2, -1);

		REQUIRE(b - a3 == -c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(-3, -2, -1, 0, 1);

		REQUIRE(b - a5 == -c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar reverse multiply", "[Vector]", TypesAll, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(1, 2, 3);
		VectorT<3> c3(4, 8, 12);

		REQUIRE(b * a3 == c3);

		VectorT<5> a5(1, 2, 3, 4, 5);
		VectorT<5> c5(4, 8, 12, 16, 20);

		REQUIRE(b * a5 == c5);
	}
}

TEST_CASE_VEC_VARIANT("Vector - Vector-scalar reverse div", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		Type b = 4;

		VectorT<3> a3(4, 8, 12);
		VectorT<3> c3(1, 1.0 / 2.0, 1.0 / 3.0);

		REQUIRE(b / a3 == c3);

		VectorT<5> a5(4, 8, 12, 16, 20);
		VectorT<5> c5(1, 1.0f / 2.0, 1.0f / 3.0, 1.0f / 4.0, 1.0f / 5.0);

		REQUIRE(b / a5 == c5);
	}
}

//------------------------------------------------------------------------------
// Swizzle
//------------------------------------------------------------------------------

#define TEST_CASE_SWIZZLE_VECTOR_OP(NAME, OPERATOR)                                                \
	TEST_CASE_VEC_VARIANT("Vector - Swizzle-vector " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                               \
		VectorT<3> v2 = { 1, 4, -2 };                                                              \
		VectorT<3> r = v1.xyz OPERATOR v2;                                                         \
		VectorT<3> e = v1 OPERATOR v2;                                                             \
		REQUIRE(r == e);                                                                           \
	}

#define TEST_CASE_VECTOR_SWIZZLE_OP(NAME, OPERATOR)                                                \
	TEST_CASE_VEC_VARIANT("Vector - Vector-swizzle " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                               \
		VectorT<3> v2 = { 1, 4, -2 };                                                              \
		VectorT<3> r = v1 OPERATOR v2.xyz;                                                         \
		VectorT<3> e = v1 OPERATOR v2;                                                             \
		REQUIRE(r == e);                                                                           \
	}

#define TEST_CASE_VECTOR_SWIZZLE_COMPOUND_OP(NAME, OPERATOR)                                       \
	TEST_CASE_VEC_VARIANT("Vector - Vector-swizzle " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                               \
		auto v1c = v1;                                                                             \
		VectorT<3> v2 = { 1, 4, -2 };                                                              \
		v1 OPERATOR v2.xyz;                                                                        \
		v1c OPERATOR v2;                                                                           \
		REQUIRE(v1 == v1c);                                                                        \
	}

#define TEST_CASE_SWIZZLE_VECTOR_COMPOUND_OP(NAME, OPERATOR)                                       \
	TEST_CASE_VEC_VARIANT("Vector - Swizzle-vector " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                               \
		auto v1c = v1;                                                                             \
		VectorT<3> v2 = { 1, 4, -2 };                                                              \
		v1.xyz OPERATOR v2;                                                                        \
		v1c OPERATOR v2;                                                                           \
		REQUIRE(v1 == v1c);                                                                        \
	}

#define TEST_CASE_SWIZZLE_SWIZZLE_OP(NAME, OPERATOR)                                                \
	TEST_CASE_VEC_VARIANT("Vector - Swizzle-swizzle " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                                \
		VectorT<3> v2 = { 1, 4, -2 };                                                               \
		VectorT<3> r = v1.xyz OPERATOR v2.xyz;                                                      \
		VectorT<3> e = v1 OPERATOR v2;                                                              \
		REQUIRE(r == e);                                                                            \
	}

#define TEST_CASE_SWIZZLE_SWIZZLE_COMPOUND_OP(NAME, OPERATOR)                                       \
	TEST_CASE_VEC_VARIANT("Vector - Swizzle-swizzle " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                                \
		auto v1c = v1;                                                                              \
		VectorT<3> v2 = { 1, 4, -2 };                                                               \
		v1.xyz OPERATOR v2.xyz;                                                                     \
		v1c OPERATOR v2;                                                                            \
		REQUIRE(v1 == v1c);                                                                         \
	}


#define TEST_CASE_SWIZZLE_SCALAR_COMPOUND_OP(NAME, OPERATOR)                                               \
	TEST_CASE_VEC_VARIANT("Vector - Swizzle-scalar compound" NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                                       \
		auto v1c = v1;                                                                                     \
		Type b = 6;                                                                                        \
		v1.xyz OPERATOR b;                                                                                 \
		v1c OPERATOR b;                                                                                    \
		REQUIRE(v1 == v1c);                                                                                \
	}

#define TEST_CASE_SWIZZLE_SCALAR_OP(NAME, OPERATOR)                                                \
	TEST_CASE_VEC_VARIANT("Vector - Swizzle-scalar " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                               \
		Type b = 6;                                                                                \
		VectorT<3> r = v1.xyz OPERATOR b;                                                          \
		VectorT<3> e = v1 OPERATOR b;                                                              \
		REQUIRE(r == e);                                                                           \
	}

#define TEST_CASE_SCALAR_SWIZZLE_OP(NAME, OPERATOR)                                                \
	TEST_CASE_VEC_VARIANT("Vector - Scalar-swizzle " NAME, "[Vector]", TypesFloating, PackedAll) { \
		VectorT<3> v1 = { 1, 2, 3 };                                                               \
		Type b = 6;                                                                                \
		VectorT<3> r = b OPERATOR v1.xyz;                                                          \
		VectorT<3> e = b OPERATOR v1;                                                              \
		REQUIRE(r == e);                                                                           \
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