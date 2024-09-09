#include <Mathter/IoStream.hpp>
#include <Mathter/Matrix/Comparison.hpp>
#include <Mathter/Quaternion/Comparison.hpp>
#include <Mathter/Transforms/Rotation3DBuilder.hpp>

#include <catch2/catch_test_macros.hpp>
#include <sstream>


using namespace mathter;


TEST_CASE("IoStream: print vector", "[IoStream]") {
	std::stringstream ss;

	SECTION("Integer") {
		const auto v = Vector(1, 2, 3);
		ss << v;
		REQUIRE(ss.str() == "[1, 2, 3]");
	}
	SECTION("Float") {
		const auto v = Vector(1.1f, 2.2f, 3.3f);
		ss << v;
		REQUIRE(ss.str() == "[1.1, 2.2, 3.3]");
	}
	SECTION("Complex") {
		using namespace std::complex_literals;
		const auto v = Vector(1.1f + 0.1if, 2.2f + 0.2if, 3.3f + 0.3if);
		ss << v;
		REQUIRE(ss.str() == "[1.1 + 0.1j, 2.2 + 0.2j, 3.3 + 0.3j]");
	}
}


TEST_CASE("IoStream: print matrix", "[IoStream]") {
	std::stringstream ss;

	const std::string_view expected = "[[1, 2, 3], [4, 5, 6]]";

	SECTION("Row major") {
		const Matrix<float, 2, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false> m = {
			1, 2, 3,
			4, 5, 6
		};
		ss << m;
		REQUIRE(ss.str() == expected);
	}
	SECTION("Column major") {
		const Matrix<float, 2, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, false> m = {
			1, 2, 3,
			4, 5, 6
		};
		ss << m;
		REQUIRE(ss.str() == expected);
	}
}


TEST_CASE("IoStream: print quaternion", "[IoStream]") {
	std::stringstream ss;

	const Quaternion<float> q = RotationAxisAngle(Normalize(Vector(3.0f, 4.0f, 0.0f)), Deg2Rad(72.0f));

	ss << q;
	REQUIRE(ss.str() == "[72 deg @ [0.6, 0.8, 0]]");
}


TEST_CASE("IoStream: parse vector", "[IoStream]") {
	std::stringstream ss;

	SECTION("Good") {
		SECTION("Integer") {
			Vector<int, 3> v;
			ss.str("[1, 2, 3]");
			ss >> v;
			REQUIRE(v == Vector(1, 2, 3));
		}
		SECTION("Float") {
			Vector<float, 3> v;
			ss.str("[1.1, 2.2, 3.3]");
			ss >> v;
			REQUIRE(v == Vector(1.1f, 2.2f, 3.3f));
		}
		SECTION("Complex") {
			using namespace std::complex_literals;
			Vector<std::complex<float>, 3> v;
			ss.str("[1.1 + 0.1j, 2.2 + 0.2j, 3.3 + 0.3j]");
			ss >> v;
			REQUIRE(v == Vector(1.1f + 0.1if, 2.2f + 0.2if, 3.3f + 0.3if));
		}
	}
	SECTION("Bad") {
		SECTION("Incorrect dimension") {
			Vector<float, 3> v;
			ss.str("[1.1, 2.2]");
			REQUIRE_THROWS(ss >> v);
		}
		SECTION("Incorrect array") {
			// TODO: these tests are not super important
			Vector<float, 3> v;
			ss.str("[1.1, 2.2, ]");
			// REQUIRE_THROWS(ss >> v);
		}
		SECTION("Incorrect array lead") {
			Vector<float, 3> v;
			ss.str("{1.1, 2.2, 3.3]");
			REQUIRE_THROWS(ss >> v);
		}
	}
}


TEST_CASE("IoStream: parse matrix", "[IoStream]") {
	using Mat = Matrix<float, 2, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;


	const auto expected = Mat{
		1, 2, 3,
		4, 5, 6
	};

	std::stringstream ss;
	ss.str("[[1, 2, 3], [4, 5, 6]]");
	Mat m;
	ss >> m;
	REQUIRE(m == expected);
}


TEST_CASE("IoStream: parse quaternion", "[IoStream]") {
	using Quat = Quaternion<float>;

	const Quat expected = RotationAxisAngle(Normalize(Vector(3.0f, 4.0f, 0.0f)), Deg2Rad(72.0f));

	std::stringstream ss;
	ss.str("[72 deg @ [0.6, 0.8, 0]]");
	Quat q;
	ss >> q;
	REQUIRE(q == expected);
}