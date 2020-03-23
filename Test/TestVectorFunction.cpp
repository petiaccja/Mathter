//==============================================================================
// This software is distributed under The Unlicense.
// For more information, please refer to <http://unlicense.org/>
//==============================================================================

#pragma warning(disable : 4244)

#include "../Mathter/Common/Approx.hpp"
#include "../Mathter/Vector.hpp"
#include "TestGenerators.hpp"

#include <Catch2/catch.hpp>

using namespace mathter;


TEST_CASE_VEC_VARIANT("Vector - Is nullvector", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1, 2, 3);
		REQUIRE(IsNullvector(a) == false);
		VectorT<3> b(0, 0, 0);
		REQUIRE(IsNullvector(b) == true);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Length", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1, 2, 3);
		REQUIRE(Length(a) == Approx(3.7416573867));


		VectorT<5> b(1, 0, 2, 0, 3);
		REQUIRE(Length(b) == Approx(3.7416573867));
	}
}


TEST_CASE_VEC_VARIANT("Vector - LengthPrecise", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1e-38f, 2e-38f, 3e-38f);
		REQUIRE(LengthPrecise(a) == Approx(3.7416573867e-38f));

		VectorT<5> b(1e+37f, 0, 2e+37f, 0, 3e+37f);
		REQUIRE(LengthPrecise(b) == Approx(3.7416573867e+37f));
	}
}


TEST_CASE_VEC_VARIANT("Vector - Normalize", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1, 2, 3);
		a = Normalize(a);
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(2 * a[0] == Approx(a[1]));
		REQUIRE(3 * a[0] == Approx(a[2]));

		VectorT<5> b(1, 0, 2, 0, 3);
		b = Normalize(b);
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(2 * b[0] == Approx(b[2]));
		REQUIRE(3 * b[0] == Approx(b[4]));
	}
}


TEST_CASE_VEC_VARIANT("Vector - SafeNormalize denorm", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(0, 1e-40, 0);
		a = SafeNormalize(a);
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(a[1] == Approx(1));

		VectorT<5> b(0, 0, 1e-40, 0, 0);
		b = SafeNormalize(b);
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(b[2] == Approx(1));
	}
}


TEST_CASE_VEC_VARIANT("Vector - SafeNormalize null", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(0, 0, 0);
		a = SafeNormalize(a);
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(a[0] == Approx(1));

		VectorT<5> b(0, 0, 0, 0, 0);
		b = SafeNormalize(b);
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(a[0] == Approx(1));
	}
}


TEST_CASE_VEC_VARIANT("Vector - SafeNormalize specific proper", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1, 2, 3);
		REQUIRE(ApproxVec(Normalize(a)) == SafeNormalize(a, VectorT<3>(0, 1, 0)));

		VectorT<5> b(1, 0, 2, 0, 3);
		REQUIRE(ApproxVec(Normalize(b)) == SafeNormalize(b, VectorT<5>(0, 1, 0, 0, 0)));
	}
}


TEST_CASE_VEC_VARIANT("Vector - SafeNormalize specific null", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(0, 0, 0);
		a = SafeNormalize(a, VectorT<3>(0,1,0));
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(a[1] == Approx(1));

		VectorT<5> b(0, 0, 0, 0, 0);
		b = SafeNormalize(b, VectorT<5>(0,1,0,0,0));
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(a[1] == Approx(1));
	}
}


TEST_CASE_VEC_VARIANT("Vector - Fill", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a;
		VectorT<3> b(4);
		Fill(a, 4);

		REQUIRE(a == b);

		VectorT<5> c;
		VectorT<5> d(4);
		Fill(c, 4);
		
		REQUIRE(c == d);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Min & Max", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1, 2, 3);
		VectorT<3> b(3, 2, 1);

		REQUIRE(Min(a,b) == VectorT<3>(1, 2, 1));
		REQUIRE(Max(a,b) == VectorT<3>(3, 2, 3));

		VectorT<5> c(1, 2, 3, 4, 5);
		VectorT<5> d(5, 4, 3, 2, 1);
		
		REQUIRE(Min(c, d) == VectorT<5>(1, 2, 3, 2, 1));
		REQUIRE(Max(c, d) == VectorT<5>(5, 4, 3, 4, 5));
	}
}


TEST_CASE_VEC_VARIANT("Vector - Dot", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1, 2, 3);
		VectorT<3> b(4, 5, 6);
		auto r1 = Dot(a, b);

		REQUIRE(r1 == 32);

		VectorT<5> c(1, 2, 3, 2, 1);
		VectorT<5> d(4, 5, 6, 5, 4);
		auto r2 = Dot(c, d);
		REQUIRE(r2 == 46);
	}
}


TEST_CASE_VEC_VARIANT("Vector - Cross", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		VectorT<3> a(1, 2, 3);
		VectorT<3> b(4, 5, 6);
		VectorT<3> r = Cross(a, b);
		VectorT<3> rexp(-3, 6, -3);

		REQUIRE(r == rexp);
	}
}


TEST_CASE_VEC_VARIANT("Vector - CrossND", "[Vector]", TypesFloating, PackedAll) {
	SECTION(SECTIONNAMEVEC) {
		// Simple 3D cross product
		VectorT<3> a(1, 2, 3);
		VectorT<3> b(4, 5, 6);
		VectorT<3> r = Cross(a, b);
		VectorT<3> rexp(-3, 6, -3);

		REQUIRE(r == rexp);

		// 2D cross product, that is, rotate by 90 degree
		VectorT<2> a2(1, 2);
		VectorT<2> r2 = Cross(a2);
		VectorT<2> r2exp(-2, 1);

		REQUIRE(ApproxVec(r2) == r2exp);

		// 4D cross product
		VectorT<4> a4(1, 2, 3, 4);
		VectorT<4> b4(4, 2, 6, 3);
		VectorT<4> c4(3, 6, 4, -9);
		VectorT<4> r4 = Cross(a4, b4, c4);

		auto dotprod = std::abs(Dot(a4, r4)) + std::abs(Dot(b4, r4)) + std::abs(Dot(c4, r4));
		REQUIRE(dotprod < 1e-5f);
	}
}