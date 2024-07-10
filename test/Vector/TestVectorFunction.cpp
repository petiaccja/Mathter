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

using namespace mathter;
using Catch::Approx;


using TypeListFloating = TestTypeList<TypesFloating, PackedAll>;


TEMPLATE_LIST_TEST_CASE("Vector - Is nullvector", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);
		REQUIRE(IsNullvector(a) == false);
		Vec3 b(0, 0, 0);
		REQUIRE(IsNullvector(b) == true);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Length", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);
		REQUIRE(Length(a) == Approx(3.7416573867));


		Vec5 b(1, 0, 2, 0, 3);
		REQUIRE(Length(b) == Approx(3.7416573867));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - LengthPrecise", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		SECTION("Underflow") {

			Vec3 a(1e-38f, 2e-38f, 3e-38f);
			REQUIRE(LengthPrecise(a) == Approx(3.7416573867e-38f));

			Vec5 b(1e+37f, 0, 2e+37f, 0, 3e+37f);
			REQUIRE(LengthPrecise(b) == Approx(3.7416573867e+37f));
		}
		SECTION("Denormal") {
			Vec3 a(1e-39f, 1e-39f, 1e-39f);
			REQUIRE(LengthPrecise(a) == Approx(1.7320508075689e-39f));
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Normalize", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);
		a = Normalize(a);
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(2 * a[0] == Approx(a[1]));
		REQUIRE(3 * a[0] == Approx(a[2]));

		Vec5 b(1, 0, 2, 0, 3);
		b = Normalize(b);
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(2 * b[0] == Approx(b[2]));
		REQUIRE(3 * b[0] == Approx(b[4]));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - SafeNormalize denorm", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(0, 1e-40, 0);
		a = SafeNormalize(a);
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(a[1] == Approx(1));

		Vec5 b(0, 0, 1e-40, 0, 0);
		b = SafeNormalize(b);
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(b[2] == Approx(1));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - SafeNormalize null", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(0, 0, 0);
		a = SafeNormalize(a);
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(a[0] == Approx(1));

		Vec5 b(0, 0, 0, 0, 0);
		b = SafeNormalize(b);
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(a[0] == Approx(1));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - SafeNormalize specific proper", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);
		REQUIRE(ApproxVec(Normalize(a)) == SafeNormalize(a, Vec3(0, 1, 0)));

		Vec5 b(1, 0, 2, 0, 3);
		REQUIRE(ApproxVec(Normalize(b)) == SafeNormalize(b, Vec5(0, 1, 0, 0, 0)));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - SafeNormalize specific null", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(0, 0, 0);
		a = SafeNormalize(a, Vec3(0, 1, 0));
		REQUIRE(Length(a) == Approx(1));
		REQUIRE(a[1] == Approx(1));

		Vec5 b(0, 0, 0, 0, 0);
		b = SafeNormalize(b, Vec5(0, 1, 0, 0, 0));
		REQUIRE(Length(b) == Approx(1));
		REQUIRE(a[1] == Approx(1));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Fill", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a;
		Vec3 b(4);
		Fill(a, 4);

		REQUIRE(a == b);

		Vec5 c;
		Vec5 d(4);
		Fill(c, 4);

		REQUIRE(c == d);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Min & Max rlementwise", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);
		Vec3 b(3, 2, 1);

		REQUIRE(Min(a, b) == Vec3(1, 2, 1));
		REQUIRE(Max(a, b) == Vec3(3, 2, 3));

		Vec5 c(1, 2, 3, 4, 5);
		Vec5 d(5, 4, 3, 2, 1);

		REQUIRE(Min(c, d) == Vec5(1, 2, 3, 2, 1));
		REQUIRE(Max(c, d) == Vec5(5, 4, 3, 4, 5));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Min reduce", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);

		REQUIRE(Min(a) == 1);

		Vec5 c(1, 2, 3, 4, 5);

		REQUIRE(Min(c) == 1);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Max reduce", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(-1, -2, -3);

		REQUIRE(Max(a) == -1);

		Vec5 c(-1, -2, -3, -4, -5);

		REQUIRE(Max(c) == -1);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Sum", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);

		REQUIRE(Sum(a) == 6);

		Vec5 c(1, 2, 3, 4, 5);

		REQUIRE(Sum(c) == 15);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Abs", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a3(1, -2, 3);
		Vec3 e3(1, 2, 3);

		REQUIRE(Abs(a3) == e3);

		Vec5 a5(1, 2, -3, 4, 5);
		Vec5 e5(1, 2, 3, 4, 5);

		REQUIRE(Abs(a5) == e5);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Dot", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);
		Vec3 b(4, 5, 6);
		auto r1 = Dot(a, b);

		REQUIRE(r1 == 32);

		Vec5 c(1, 2, 3, 2, 1);
		Vec5 d(4, 5, 6, 5, 4);
		auto r2 = Dot(c, d);
		REQUIRE(r2 == 46);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Cross", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;

		Vec3 a(1, 2, 3);
		Vec3 b(4, 5, 6);
		Vec3 r = Cross(a, b);
		Vec3 rexp(-3, 6, -3);

		REQUIRE(r == rexp);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - CrossND", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec2 = typename TestType::template Vector<2>;
		using Vec3 = typename TestType::template Vector<3>;
		using Vec4 = typename TestType::template Vector<4>;

		// Simple 3D cross product
		Vec3 a(1, 2, 3);
		Vec3 b(4, 5, 6);
		Vec3 r = Cross(a, b);
		Vec3 rexp(-3, 6, -3);

		REQUIRE(r == rexp);

		// 2D cross product, that is, rotate by 90 degree
		Vec2 a2(1, 2);
		Vec2 r2 = Cross(a2);
		Vec2 r2exp(-2, 1);

		REQUIRE(ApproxVec(r2) == r2exp);

		// 4D cross product
		Vec4 a4(1, 2, 3, 4);
		Vec4 b4(4, 2, 6, 3);
		Vec4 c4(3, 6, 4, -9);
		Vec4 r4 = Cross(a4, b4, c4);

		auto dotprod = std::abs(Dot(a4, r4)) + std::abs(Dot(b4, r4)) + std::abs(Dot(c4, r4));
		REQUIRE(dotprod < 1e-5f);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Distance", "[Vector]", TypeListFloating) {
	SECTION(TestType::Name()) {
		using Vec3 = typename TestType::template Vector<3>;
		using Vec5 = typename TestType::template Vector<5>;
		using Type = scalar_type_t<Vec3>;

		Vec3 a(1, 2, 3);
		Vec3 b(4, 5, 4);
		Type expected3 = Type(std::sqrt(3 * 3 + 3 * 3 + 1 * 1));
		REQUIRE(expected3 == Distance(a, b));

		Vec5 d(4, 5, 4, 2, 8);
		Vec5 c(1, 2, 3, 9, 2);
		Type expected5 = Type(std::sqrt(3 * 3 + 3 * 3 + 1 * 1 + 7 * 7 + 6 * 6));
		REQUIRE(expected5 == Distance(c, d));
	}
}