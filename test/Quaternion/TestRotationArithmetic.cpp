// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"
#include "../Rotation.hpp"

#include <Mathter/Quaternion/RotationArithmetic.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif


using namespace mathter;
using namespace test_util;



//------------------------------------------------------------------------------
// Quat x vec
//------------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("Quaternion - Multiplication (quat x vector)", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												VectorCaseList<ScalarsFloating, PackingsAll>>{})) {
	using Quat = typename TestType::Lhs::Quat;
	using Vec3 = typename TestType::Rhs::template Vector<3>;
	using Vec4 = typename TestType::Rhs::template Vector<4>;
	using Scalar = scalar_type_t<Quat>;
	using ResultScalar = common_arithmetic_type_t<scalar_type_t<Quat>, scalar_type_t<Vec3>>;
	using ResultVec3 = Vector<ResultScalar, 3, false>;

	const auto theta = Scalar(-0.78);
	const auto axis = Normalize(Vector<Scalar, 3, is_packed_v<Quat>>(1, 2, 3));
	const auto scalar = std::cos(theta / 2);
	const auto vector = std::sin(theta / 2) * axis;

	const Quat q(scalar, vector);
	const Vec3 v(4, 1, 2);
	const Vec4 vh(v, 1);
	const auto expected = Rotate(ResultVec3(v), ResultVec3(axis), ResultScalar(theta));

	SECTION("operator()") {
		const auto result = q(v);
		static_assert(is_vector_v<std::decay_t<decltype(result)>>);
		static_assert(dimension_v<std::decay_t<decltype(result)>> == 3);
		REQUIRE((result == test_util::Approx(expected, 1e-6f)));
	}

	// These are testing deprecated functionality.
	SECTION("Quat x Vec3") {
		const auto result = q * v;
		static_assert(is_vector_v<std::decay_t<decltype(result)>>);
		static_assert(dimension_v<std::decay_t<decltype(result)>> == 3);
		REQUIRE((result == test_util::Approx(expected, 1e-6f)));
	}
	SECTION("Vec3 x Quat") {
		const auto result = v * q;
		static_assert(is_vector_v<std::decay_t<decltype(result)>>);
		static_assert(dimension_v<std::decay_t<decltype(result)>> == 3);
		REQUIRE((result == test_util::Approx(expected, 1e-6f)));
	}
	SECTION("Vec3 x Quat") {
		auto copy = v;
		REQUIRE(&(copy *= q) == &copy);
		REQUIRE((copy == test_util::Approx(expected, 1e-6f)));
	}
}