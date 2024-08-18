// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)
#if defined(__clang__)
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif


#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Quaternion/Quaternion.hpp>
#include <Mathter/Transforms/Rotation3DBuilder.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Quaternion - Construct default", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	// The elements are initialized by the vector default ctor, so no need to test that.
	// Still important to test that this compiles.
	using Quat = typename TestType::Quat;
	Quat q;
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Construct from vec4", "[Quaternion]",
						decltype(BinaryCaseList<
								 QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
								 VectorCaseList<ScalarsFloating, PackingsAll>>{})) {
	// This test is a prerequisite for testing element accessors. Do not just delete this.
	using Quat = typename TestType::Lhs::Quat;
	using Vec = typename TestType::Rhs::template Vector<4>;

	const Vec v = { 1, 2, 3, 4 };
	const Quat q(v);

	REQUIRE(q.elements.array[0] == 1);
	REQUIRE(q.elements.array[1] == 2);
	REQUIRE(q.elements.array[2] == 3);
	REQUIRE(q.elements.array[3] == 4);
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Element access", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	// This test relies on the vec4 converting constructor being tested properly.
	using Quat = typename TestType::Quat;
	using Vec = Vector<scalar_type_t<Quat>, 4, is_packed_v<Quat>>;

	const auto v = layout_v<Quat> == eQuaternionLayout::SCALAR_FIRST ? Vec{ 1, 2, 3, 4 } : Vec{ 2, 3, 4, 1 };
	const Quat q(v);

	REQUIRE(q.s == 1);
	REQUIRE(q.i == 2);
	REQUIRE(q.j == 3);
	REQUIRE(q.k == 4);

	REQUIRE(q.w == 1);
	REQUIRE(q.x == 2);
	REQUIRE(q.y == 3);
	REQUIRE(q.z == 4);

	REQUIRE(q.scalar == 1);
	REQUIRE(q.vector[0] == 2);
	REQUIRE(q.vector[1] == 3);
	REQUIRE(q.vector[2] == 4);

	REQUIRE(q.canonical[0] == 1);
	REQUIRE(q.canonical[1] == 2);
	REQUIRE(q.canonical[2] == 3);
	REQUIRE(q.canonical[3] == 4);
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Construct from 4 coefficients", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat q(1, 2, 3, 4);

	REQUIRE(q.scalar == 1);
	REQUIRE(q.vector[0] == 2);
	REQUIRE(q.vector[1] == 3);
	REQUIRE(q.vector[2] == 4);
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Construct from scalar and vector", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat q(1, Vector(2, 3, 4));

	REQUIRE(q.scalar == 1);
	REQUIRE(q.vector[0] == 2);
	REQUIRE(q.vector[1] == 3);
	REQUIRE(q.vector[2] == 4);
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Construct from rotation matrix", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>>{})) {
	using Quat = typename TestType::Lhs::Quat;
	using Vec = Vector<scalar_type_t<Quat>, 3, is_packed_v<Quat>>;

	const auto builder = RotationAxisAngle(Normalize(Vec(1, 2, 3)), scalar_type_t<Quat>(0.7854626));

	SECTION("3x3") {
		using Mat = typename TestType::Rhs::template Matrix<3, 3>;
		const Mat m = builder;
		const Quat result(m);
		const Quat expected = builder;
		REQUIRE(result == test_util::Approx(expected, 2e-6f));
	}
	SECTION("4x3") {
		using Mat = typename TestType::Rhs::template Matrix<4, 3>;
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			const Mat m = builder;
			const Quat result(m);
			const Quat expected = builder;
			REQUIRE(result == test_util::Approx(expected, 2e-6f));
		}
	}
	SECTION("3x4") {
		using Mat = typename TestType::Rhs::template Matrix<3, 4>;
		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			const Mat m = builder;
			const Quat result(m);
			const Quat expected = builder;
			REQUIRE(result == test_util::Approx(expected, 2e-6f));
		}
	}
	SECTION("4x4") {
		using Mat = typename TestType::Rhs::template Matrix<4, 4>;
		const Mat m = builder;
		const Quat result(m);
		const Quat expected = builder;
		REQUIRE(result == test_util::Approx(expected, 2e-6f));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Convert to rotation matrix", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												MatrixCaseList<ScalarsFloating, OrdersAll, LayoutsAll, PackingsAll>>{})) {
	using Quat = typename TestType::Lhs::Quat;
	using Vec = Vector<scalar_type_t<Quat>, 3, is_packed_v<Quat>>;

	const auto builder = RotationAxisAngle(Normalize(Vec(1, 2, 3)), scalar_type_t<Quat>(0.7854626));

	SECTION("3x3") {
		using Mat = typename TestType::Rhs::template Matrix<3, 3>;
		const Quat q = builder;
		const Mat result(q);
		const Mat expected = builder;
		REQUIRE(result == test_util::Approx(expected, 2e-6f));
	}
	SECTION("4x3") {
		using Mat = typename TestType::Rhs::template Matrix<4, 3>;
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			const Quat q = builder;
			const Mat result(q);
			const Mat expected = builder;
			REQUIRE(result == test_util::Approx(expected, 2e-6f));
		}
	}
	SECTION("3x4") {
		using Mat = typename TestType::Rhs::template Matrix<3, 4>;
		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			const Quat q = builder;
			const Mat result(q);
			const Mat expected = builder;
			REQUIRE(result == test_util::Approx(expected, 2e-6f));
		}
	}
	SECTION("4x4") {
		using Mat = typename TestType::Rhs::template Matrix<4, 4>;
		const Quat q = builder;
		const Mat result(q);
		const Mat expected = builder;
		REQUIRE(result == test_util::Approx(expected, 2e-6f));
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Convert to vector", "[Quaternion]",
						decltype(BinaryCaseList<QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>,
												VectorCaseList<ScalarsFloating, PackingsAll>>{})) {
	using Quat = typename TestType::Lhs::Quat;
	using Vec = typename TestType::Rhs::template Vector<4>;

	const Quat q(1, 2, 3, 4);
	const Vec v(q);

	if constexpr (layout_v<Quat> == eQuaternionLayout::SCALAR_FIRST) {
		REQUIRE(v[0] == 1);
		REQUIRE(v[1] == 2);
		REQUIRE(v[2] == 3);
		REQUIRE(v[3] == 4);
	}
	else {
		REQUIRE(v[0] == 2);
		REQUIRE(v[1] == 3);
		REQUIRE(v[2] == 4);
		REQUIRE(v[3] == 1);
	}
}


TEMPLATE_LIST_TEST_CASE("Quaternion - Scalar part & vector part", "[Quaternion]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;

	const Quat q(1, 2, 3, 4);

	REQUIRE(q.ScalarPart() == 1);
	REQUIRE(q.VectorPart()[0] == 2);
	REQUIRE(q.VectorPart()[1] == 3);
	REQUIRE(q.VectorPart()[2] == 4);
}