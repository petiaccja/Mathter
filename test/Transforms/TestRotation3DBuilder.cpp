// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2024 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../ApplyTransform.hpp"
#include "../Approx.hpp"
#include "../Cases.hpp"
#include "../Rotation.hpp"

#include <Mathter/Quaternion/RotationArithmetic.hpp>
#include <Mathter/Transforms/Rotation3DBuilder.hpp>

#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;
using namespace std::complex_literals;


template <class Mat, class Builder, class Vec, class Scalar>
void TestCaseRotationMatrix(const Builder& builder, const Vec& testPoint, const Vec& axis, const Scalar& angle) {
	const auto expected = test_util::Rotate(testPoint, axis, angle);

	const Mat m = builder;
	const auto result = ApplyTransform(m, testPoint);
	REQUIRE(result == test_util::Approx(expected, 2e-6f));
}


template <class Quat, class Builder, class Vec, class Scalar>
void TestCaseRotationQuat(const Builder& builder, const Vec& testPoint, const Vec& axis, const Scalar& angle) {
	const auto expected = test_util::Rotate(testPoint, axis, angle);

	const Quat q = builder;
	const auto result = q * testPoint;
	REQUIRE(result == test_util::Approx(expected, 2e-6f));
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- Single axis matrix", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Scalar = scalar_type_t<ExampleMat>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;

	const Vec testPoint = { 1, 2, 3 };
	const Scalar angle = 0.64350110879328;


	SECTION("X axis") {
		const Vec axis = { 1, 0, 0 };
		const auto builder = RotationX(angle);

		SECTION("3x3") {
			using Mat = typename TestType::template Matrix<3, 3>;
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
		SECTION("4x4") {
			using Mat = typename TestType::template Matrix<4, 4>;
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
		SECTION("3x4") {
			using Mat = typename TestType::template Matrix<3, 4>;
			if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
				TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
			}
		}
		SECTION("4x3") {
			using Mat = typename TestType::template Matrix<4, 3>;
			if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
				TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
			}
		}
	}

	SECTION("Y axis") {
		const Vec axis = { 0, 1, 0 };
		const auto builder = RotationY(angle);

		SECTION("3x3") {
			using Mat = typename TestType::template Matrix<3, 3>;
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
		SECTION("4x4") {
			using Mat = typename TestType::template Matrix<4, 4>;
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
		SECTION("3x4") {
			using Mat = typename TestType::template Matrix<3, 4>;
			if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
				TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
			}
		}
		SECTION("4x3") {
			using Mat = typename TestType::template Matrix<4, 3>;
			if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
				TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
			}
		}
	}

	SECTION("Z axis") {
		const Vec axis = { 0, 0, 1 };
		const auto builder = RotationZ(angle);

		SECTION("3x3") {
			using Mat = typename TestType::template Matrix<3, 3>;
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
		SECTION("4x4") {
			using Mat = typename TestType::template Matrix<4, 4>;
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
		SECTION("3x4") {
			using Mat = typename TestType::template Matrix<3, 4>;
			if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
				TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
			}
		}
		SECTION("4x3") {
			using Mat = typename TestType::template Matrix<4, 3>;
			if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
				TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
			}
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- Single axis quaternion", "[Transforms]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;
	using Scalar = scalar_type_t<Quat>;
	using Vec = Vector<Scalar, 3, is_packed_v<Quat>>;

	const Vec testPoint = { 1, 2, 3 };
	const Scalar angle = 0.64350110879328;

	SECTION("X axis") {
		const Vec axis = { 1, 0, 0 };
		const auto builder = RotationX(angle);
		TestCaseRotationQuat<Quat>(builder, testPoint, axis, angle);
	}

	SECTION("Y axis") {
		const Vec axis = { 0, 1, 0 };
		const auto builder = RotationY(angle);
		TestCaseRotationQuat<Quat>(builder, testPoint, axis, angle);
	}

	SECTION("Z axis") {
		const Vec axis = { 0, 0, 1 };
		const auto builder = RotationZ(angle);
		TestCaseRotationQuat<Quat>(builder, testPoint, axis, angle);
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- RPY matrix", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Scalar = remove_complex_t<scalar_type_t<ExampleMat>>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;
	using Quat = Quaternion<Scalar>;

	const std::array<Scalar, 3> angles = { 0.5f, 0.4f, -0.95f };

	const auto builder = RotationRPY(angles[0], angles[1], angles[2]);

	// This test primarily verifies if the order of axes is correct.
	const Quat r0 = RotationX(angles[0]);
	const Quat r1 = RotationY(angles[1]);
	const Quat r2 = RotationZ(angles[2]);

	const auto p0 = Vec(1, 2, 3);
	const auto p1 = r0 * p0;
	const auto p2 = r1 * p1;
	const auto p3 = r2 * p2;

	SECTION("3x3") {
		using Mat = typename TestType::template Matrix<3, 3>;
		const Mat r = builder;
		const auto result = ApplyTransform(r, p0);
		REQUIRE(result == test_util::Approx(p3));
	}
	SECTION("4x4") {
		using Mat = typename TestType::template Matrix<4, 4>;
		const Mat r = builder;
		const auto result = ApplyTransform(r, p0);
		REQUIRE(result == test_util::Approx(p3));
	}
	SECTION("3x4") {
		using Mat = typename TestType::template Matrix<3, 4>;
		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			const Mat r = builder;
			const auto result = ApplyTransform(r, p0);
			REQUIRE(result == test_util::Approx(p3));
		}
	}
	SECTION("4x3") {
		using Mat = typename TestType::template Matrix<4, 3>;
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			const Mat r = builder;
			const auto result = ApplyTransform(r, p0);
			REQUIRE(result == test_util::Approx(p3));
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- Euler angles matrix", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Scalar = remove_complex_t<scalar_type_t<ExampleMat>>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;
	using Quat = Quaternion<Scalar>;

	const std::array<Scalar, 3> angles = { 0.5f, 0.4f, -0.95f };

	const auto builder = RotationEuler(angles[0], angles[1], angles[2]);

	// This test primarily verifies if the order of axes is correct.
	const Quat r0 = RotationZ(angles[0]);
	const Quat r1 = RotationX(angles[1]);
	const Quat r2 = RotationZ(angles[2]);

	const auto p0 = Vec(1, 2, 3);
	const auto p1 = r0 * p0;
	const auto p2 = r1 * p1;
	const auto p3 = r2 * p2;

	SECTION("3x3") {
		using Mat = typename TestType::template Matrix<3, 3>;
		const Mat r = builder;
		const auto result = ApplyTransform(r, p0);
		REQUIRE(result == test_util::Approx(p3));
	}
	SECTION("4x4") {
		using Mat = typename TestType::template Matrix<4, 4>;
		const Mat r = builder;
		const auto result = ApplyTransform(r, p0);
		REQUIRE(result == test_util::Approx(p3));
	}
	SECTION("3x4") {
		using Mat = typename TestType::template Matrix<3, 4>;
		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			const Mat r = builder;
			const auto result = ApplyTransform(r, p0);
			REQUIRE(result == test_util::Approx(p3));
		}
	}
	SECTION("4x3") {
		using Mat = typename TestType::template Matrix<4, 3>;
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			const Mat r = builder;
			const auto result = ApplyTransform(r, p0);
			REQUIRE(result == test_util::Approx(p3));
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- RPY quaternion", "[Transforms]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;
	using Scalar = scalar_type_t<Quat>;
	using Vec = Vector<Scalar, 3, is_packed_v<Quat>>;

	const std::array<Scalar, 3> angles = { 0.5f, 0.4f, -0.95f };

	// This test primarily verifies if the order of axes is correct.
	const Quat r0 = RotationX(angles[0]);
	const Quat r1 = RotationY(angles[1]);
	const Quat r2 = RotationZ(angles[2]);
	const Quat r = RotationRPY(angles[0], angles[1], angles[2]);

	const auto p0 = Vec(1, 2, 3);
	const auto p1 = r0 * p0;
	const auto p2 = r1 * p1;
	const auto p3 = r2 * p2;

	const auto result = r * p0;
	REQUIRE(result == test_util::Approx(p3));
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- Euler angles quaternion", "[Transforms]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;
	using Scalar = scalar_type_t<Quat>;
	using Vec = Vector<Scalar, 3, is_packed_v<Quat>>;

	const std::array<Scalar, 3> angles = { 0.5f, 0.4f, -0.95f };

	// This test primarily verifies if the order of axes is correct.
	const Quat r0 = RotationZ(angles[0]);
	const Quat r1 = RotationX(angles[1]);
	const Quat r2 = RotationZ(angles[2]);
	const Quat r = RotationEuler(angles[0], angles[1], angles[2]);

	const auto p0 = Vec(1, 2, 3);
	const auto p1 = r0 * p0;
	const auto p2 = r1 * p1;
	const auto p3 = r2 * p2;

	const auto result = r * p0;
	REQUIRE(result == test_util::Approx(p3));
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- Axis-angle matrix", "[Transforms]",
						decltype(MatrixCaseList<ScalarsFloatingAndComplex, OrdersAll, LayoutsAll, PackingsAll>{})) {
	using ExampleMat = typename TestType::template Matrix<1, 1>;
	using Scalar = remove_complex_t<scalar_type_t<ExampleMat>>;
	using Vec = Vector<Scalar, 3, is_packed_v<ExampleMat>>;
	using Quat = Quaternion<Scalar>;

	const Scalar angle = 0.577215664901532;
	const Vec axis = Normalize(Vec(1, 2, 3));
	const Vec testPoint = { 4, 2, 3 };
	const auto builder = RotationAxisAngle(axis, angle);

	SECTION("3x3") {
		using Mat = typename TestType::template Matrix<3, 3>;
		TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
	}
	SECTION("4x4") {
		using Mat = typename TestType::template Matrix<4, 4>;
		TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
	}
	SECTION("3x4") {
		using Mat = typename TestType::template Matrix<3, 4>;
		if constexpr (order_v<Mat> == eMatrixOrder::PRECEDE_VECTOR) {
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
	}
	SECTION("4x3") {
		using Mat = typename TestType::template Matrix<4, 3>;
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			TestCaseRotationMatrix<Mat>(builder, testPoint, axis, angle);
		}
	}
}


TEMPLATE_LIST_TEST_CASE("Transform: Rotation 3D -- Axis-angle quaternion", "[Transforms]",
						decltype(QuaternionCaseList<ScalarsFloating, QuatLayoutsAll, PackingsAll>{})) {
	using Quat = typename TestType::Quat;
	using Scalar = remove_complex_t<scalar_type_t<Quat>>;
	using Vec = Vector<Scalar, 3, is_packed_v<Quat>>;

	const Scalar angle = 0.577215664901532;
	const Vec axis = Normalize(Vec(1, 2, 3));
	const Vec testPoint = { 4, 2, 3 };
	const auto builder = RotationAxisAngle(axis, angle);

	TestCaseRotationQuat<Quat>(builder, testPoint, axis, angle);
}