#include "../Approx.hpp"
#include "../Rotation.hpp"

#include <catch2/catch_test_macros.hpp>


using namespace mathter;


TEST_CASE("TestUtils - Rotation", "[TestUtils]") {
	using Vec = Vector<float, 3>;

	SECTION("Aligned X") {
		const Vec axis = { 1, 0, 0 }; // X axis
		const float angle = 3.141592653f / 4.0f; // 45 degrees
		const Vec point = { 2, 0, 1 };
		const Vec expected = { 2, -0.70710678118655f, 0.70710678118655f };

		const auto result = test_util::Rotate(point, axis, angle);
		REQUIRE(result == test_util::Approx(expected));
	}
	SECTION("Aligned Y") {
		const Vec axis = { 0, 1, 0 }; // Y axis
		const float angle = 3.141592653f / 4.0f; // 45 degrees
		const Vec point = { 0, 2, 1 };
		const Vec expected = { 0.70710678118655f, 2, 0.70710678118655f };

		const auto result = test_util::Rotate(point, axis, angle);
		REQUIRE(result == test_util::Approx(expected));
	}
	SECTION("Aligned Z") {
		const Vec axis = { 0, 0, 1 }; // Y axis
		const float angle = 3.141592653f / 4.0f; // 45 degrees
		const Vec point = { 0, 1, 2 };
		const Vec expected = { -0.70710678118655f, 0.70710678118655f, 2 };

		const auto result = test_util::Rotate(point, axis, angle);
		REQUIRE(result == test_util::Approx(expected));
	}
}
