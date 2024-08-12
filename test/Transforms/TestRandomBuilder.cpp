// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Transforms/RandomBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <random>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


TEMPLATE_LIST_TEST_CASE("Transform: Random - Matrix", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using M23 = typename TestType::template Matrix<2, 3>;
	using Scalar = scalar_type_t<M23>;
	using Real = remove_complex_t<Scalar>;

	std::mt19937_64 rne;
	std::uniform_real_distribution<float> rng(5.0f, 6.0f);

	const M23 value = Random(rng, rne);
	REQUIRE(std::abs(value(0, 0) - Real(5)) <= Real(1));
	REQUIRE(std::abs(value(0, 1) - Real(5)) <= Real(1));
	REQUIRE(std::abs(value(0, 2) - Real(5)) <= Real(1));
	REQUIRE(std::abs(value(1, 0) - Real(5)) <= Real(1));
	REQUIRE(std::abs(value(1, 1) - Real(5)) <= Real(1));
	REQUIRE(std::abs(value(1, 2) - Real(5)) <= Real(1));
}


TEMPLATE_LIST_TEST_CASE("Transform: Random - Quaternion", "[Transforms]",
						decltype(VectorCaseList<ScalarsFloatingAndComplex, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;
	using Scalar = scalar_type_t<Vec>;
	using Real = remove_complex_t<Scalar>;

	std::mt19937_64 rne;
	std::uniform_real_distribution<float> rng(5.0f, 6.0f);

	const Vec value = Random(rng, rne);
	REQUIRE(std::abs(value[0] - Real(5)) <= Real(1));
	REQUIRE(std::abs(value[1] - Real(5)) <= Real(1));
	REQUIRE(std::abs(value[2] - Real(5)) <= Real(1));
}


TEMPLATE_LIST_TEST_CASE("Transform: Random - Quaternion", "[Transforms]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;
	using Scalar = scalar_type_t<Quat>;

	std::mt19937_64 rne;
	std::uniform_real_distribution<float> rng(5.0f, 6.0f);

	const Quat value = Random(rng, rne);
	REQUIRE(std::abs(value.s - Scalar(5)) <= Scalar(1));
	REQUIRE(std::abs(value.i - Scalar(5)) <= Scalar(1));
	REQUIRE(std::abs(value.j - Scalar(5)) <= Scalar(1));
	REQUIRE(std::abs(value.k - Scalar(5)) <= Scalar(1));
}