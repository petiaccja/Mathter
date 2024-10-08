// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Quaternion/Arithmetic.hpp>
#include <Mathter/Quaternion/Math.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


template <class T, eQuaternionLayout Layout, bool Packed>
auto ExpSeries(const Quaternion<T, Layout, Packed>& q) {
	uint64_t nFactorial = 1;
	Quaternion<T, Layout, Packed> qPowN = T(1);
	Quaternion<T, Layout, Packed> sum = T(0);
	for (uint64_t n = 0; n < 20; ++n) {
		nFactorial *= std::max(uint64_t(1), n);
		sum += qPowN / T(nFactorial);
		qPowN *= q;
	}
	return sum;
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Axis & angle", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;
	using Scalar = scalar_type_t<Quat>;

	SECTION("Proper (0<=theta<=pi") {
		const auto theta = Scalar(1.2);
		const auto axis = Normalize(Vector<Scalar, 3, is_packed_v<Quat>>(1, 2, 3));
		const auto scalar = std::cos(theta / 2);
		const auto vector = std::sin(theta / 2) * axis;

		const Quat q(scalar, vector);

		REQUIRE(Angle(q) == Catch::Approx(theta));
		REQUIRE((Axis(q) == test_util::Approx(axis)));
	}
	SECTION("Proper (-pi<=theta<=0") {
		const auto theta = Scalar(-1.2);
		const auto axis = Normalize(Vector<Scalar, 3, is_packed_v<Quat>>(1, 2, 3));
		const auto scalar = std::cos(theta / 2);
		const auto vector = std::sin(theta / 2) * axis;

		const Quat q(scalar, vector);

		REQUIRE(Angle(q) == Catch::Approx(-theta));
		REQUIRE((Axis(q) == test_util::Approx(-axis)));
	}
	SECTION("Improper (pi<=theta") {
		const auto theta = Scalar(-5.0);
		const auto axis = Normalize(Vector<Scalar, 3, is_packed_v<Quat>>(1, 2, 3));
		const auto scalar = std::cos(theta / 2);
		const auto vector = std::sin(theta / 2) * axis;

		const Quat q(scalar, vector);

		REQUIRE(Angle(q) == Catch::Approx(-theta));
		REQUIRE((Axis(q) == test_util::Approx(-axis)));
	}
	SECTION("Improper (theta<=-pi") {
		const auto theta = Scalar(5.0);
		const auto axis = Normalize(Vector<Scalar, 3, is_packed_v<Quat>>(1, 2, 3));
		const auto scalar = std::cos(theta / 2);
		const auto vector = std::sin(theta / 2) * axis;

		const Quat q(scalar, vector);

		REQUIRE(Angle(q) == Catch::Approx(theta));
		REQUIRE((Axis(q) == test_util::Approx(axis)));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - LengthSquared", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat q(2, 5, 14, 0);
	REQUIRE(LengthSquared(q) == Catch::Approx(225));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Length", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat q(2, 5, 14, 0);
	REQUIRE(Length(q) == Catch::Approx(15));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - LengthPrecise", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	SECTION("Underflow") {
		const Quat q(2e-30f, 5e-30f, 14e-30f, 0.0f);
		REQUIRE(LengthPrecise(q) == Catch::Approx(15e-30));
	}
	SECTION("Overflow") {
		const Quat q(2e+20f, 5e+20f, 14e+20f, 0.0f);
		REQUIRE(LengthPrecise(q) == Catch::Approx(15e+20));
	}
	SECTION("Zero") {
		const Quat q(0.0f, 0.0f, 0.0f, 0.0f);
		REQUIRE(LengthPrecise(q) == Catch::Approx(0));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Normalize", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat q(2, 5, 14, 0);
	REQUIRE(Normalize(q) == test_util::Approx(Quat(2.f / 15, 5.f / 15, 14.f / 15, 0), 1e-6f));
}


TEMPLATE_LIST_TEST_CASE("Quaternion - NormalizePrecise", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	SECTION("Underflow") {
		const Quat q(2e-30f, 5e-30f, 14e-30f, 0.0f);
		REQUIRE(NormalizePrecise(q) == test_util::Approx(Quat(2.f / 15, 5.f / 15, 14.f / 15, 0), 1e-6f));
	}
	SECTION("Overflow") {
		const Quat q(2e+20f, 5e+20f, 14e+20f, 0.0f);
		REQUIRE(NormalizePrecise(q) == test_util::Approx(Quat(2.f / 15, 5.f / 15, 14.f / 15, 0), 1e-6f));
	}
	SECTION("Zero") {
		const Quat q(0.0f, 0.0f, 0.0f, 0.0f);
		REQUIRE(NormalizePrecise(q) == test_util::Approx(Quat(1, 0, 0, 0), 1e-6f));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Abs", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	SECTION("Underflow") {
		const Quat q(2e-30f, 5e-30f, 14e-30f, 0.0f);
		REQUIRE(Abs(q) == Catch::Approx(15e-30));
	}
	SECTION("Overflow") {
		const Quat q(2e+20f, 5e+20f, 14e+20f, 0.0f);
		REQUIRE(Abs(q) == Catch::Approx(15e+20));
	}
	SECTION("Zero") {
		const Quat q(0.0f, 0.0f, 0.0f, 0.0f);
		REQUIRE(Abs(q) == Catch::Approx(0));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Conj", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat q(2, 5, 14, 1);
	const auto result = Conj(q);
	REQUIRE(result.s == 2);
	REQUIRE(result.i == -5);
	REQUIRE(result.j == -14);
	REQUIRE(result.k == -1);
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Exp", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	SECTION("General") {
		const Quat q(2, 5, 3, 1);
		const auto result = Exp(q);
		const auto expected = ExpSeries(q);
		REQUIRE(result == test_util::Approx(expected, 1e-3f)); // Series method is rather inaccurate.
	}
	SECTION("Vector") {
		const Quat q(0, 5, 3, 1);
		const auto result = Exp(q);
		const auto expected = ExpSeries(q);
		REQUIRE(result == test_util::Approx(expected, 1e-2f)); // Series method is rather inaccurate.
	}
	SECTION("Scalar") {
		const Quat q(2);
		const auto result = Exp(q);
		const auto expected = Quat(std::exp(scalar_type_t<Quat>(2)));
		REQUIRE(result == test_util::Approx(expected));
	}
	SECTION("Zero") {
		const Quat q(0);
		const auto result = Exp(q);
		const auto expected = Quat(1);
		REQUIRE(result == test_util::Approx(expected));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Log", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	SECTION("General") {
		const Quat q(2, 5, 3, 1);
		const auto result = Exp(Log(q));
		REQUIRE(result == test_util::Approx(q));
	}
	SECTION("Vector") {
		const Quat q(0, 5, 3, 1);
		const auto result = Exp(Log(q));
		REQUIRE(result == test_util::Approx(q));
	}
	SECTION("Scalar") {
		const Quat q(2);
		const auto result = Exp(Log(q));
		REQUIRE(result == test_util::Approx(q));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Pow", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	SECTION("General") {
		const Quat q(2, 5, 3, 1);
		const auto result = Pow(q, scalar_type_t<Quat>(2));
		REQUIRE(result == test_util::Approx(q * q));
	}
	SECTION("Vector") {
		const Quat q(0, 5, 3, 1);
		const auto result = Pow(q, scalar_type_t<Quat>(2));
		REQUIRE(result == test_util::Approx(q * q));
	}
	SECTION("Scalar") {
		const Quat q(2);
		const auto result = Pow(q, scalar_type_t<Quat>(2));
		REQUIRE(result == test_util::Approx(q * q));
	}
	SECTION("Zero") {
		const Quat q(2);
		const auto result = Pow(q, scalar_type_t<Quat>(2));
		REQUIRE(result == test_util::Approx(q * q));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Inverse", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const auto unit = Quat(1);

	const Quat q(2, 5, 3, 1);
	const auto inverse = Inverse(q);
	REQUIRE(q * inverse == test_util::Approx(unit));
}