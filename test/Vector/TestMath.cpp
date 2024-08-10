// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Approx.hpp"
#include "../Cases.hpp"

#include <Mathter/Vector/Comparison.hpp>
#include <Mathter/Vector/Math.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>


using namespace mathter;
using namespace test_util;


TEMPLATE_LIST_TEST_CASE("Vector - Min / Max (elementwise binary)", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatAndInt32, PackingsAll>,
												VectorCaseList<ScalarsFloatAndInt32, PackingsAll>>{})) {
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;

	const VecLhs lhs = { -1, 2, 3 };
	const VecRhs rhs = { -2, 9, 0 };

	SECTION("Min") {
		const auto result = Min(lhs, rhs);
		REQUIRE(result[0] == -2);
		REQUIRE(result[1] == 2);
		REQUIRE(result[2] == 0);
	}
	SECTION("Max") {
		const auto result = Max(lhs, rhs);
		REQUIRE(result[0] == -1);
		REQUIRE(result[1] == 9);
		REQUIRE(result[2] == 3);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Min / Max (reduction)", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2, 3 };

	SECTION("Min") {
		REQUIRE(Min(value) == -1);
	}
	SECTION("Max") {
		REQUIRE(Max(value) == 3);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Abs (real)", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2, 3 };
	const auto result = Abs(value);
	REQUIRE(result[0] == 1);
	REQUIRE(result[1] == 2);
	REQUIRE(result[2] == 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - Abs (complex)", "[Vector]",
						decltype(VectorCaseList<ScalarsComplex, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2if, 3 };
	const auto result = Abs(value);
	using C = scalar_type_t<Vec>;
	using S = remove_complex_t<C>;
	static_assert(std::is_same_v<S, scalar_type_t<std::decay_t<decltype(result)>>>);
	REQUIRE(result[0] == 1);
	REQUIRE(result[1] == 2);
	REQUIRE(result[2] == 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - Real (real)", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2, 3 };
	const auto result = Real(value);
	REQUIRE(result[0] == -1);
	REQUIRE(result[1] == 2);
	REQUIRE(result[2] == 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - Real (complex)", "[Vector]",
						decltype(VectorCaseList<ScalarsComplex, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2if, 3 };
	const auto result = Real(value);
	using C = scalar_type_t<Vec>;
	using S = remove_complex_t<C>;
	static_assert(std::is_same_v<S, scalar_type_t<std::decay_t<decltype(result)>>>);
	REQUIRE(result[0] == -1);
	REQUIRE(result[1] == 0);
	REQUIRE(result[2] == 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - Imag (real)", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2, 3 };
	const auto result = Imag(value);
	REQUIRE(result[0] == 0);
	REQUIRE(result[1] == 0);
	REQUIRE(result[2] == 0);
}


TEMPLATE_LIST_TEST_CASE("Vector - Imag (complex)", "[Vector]",
						decltype(VectorCaseList<ScalarsComplex, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2if, 3 };
	const auto result = Imag(value);
	using C = scalar_type_t<Vec>;
	using S = remove_complex_t<C>;
	static_assert(std::is_same_v<S, scalar_type_t<std::decay_t<decltype(result)>>>);
	REQUIRE(result[0] == 0);
	REQUIRE(result[1] == 2);
	REQUIRE(result[2] == 0);
}


TEMPLATE_LIST_TEST_CASE("Vector - Conj (real)", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2, 3 };
	const auto result = Conj(value);
	REQUIRE(result[0] == -1);
	REQUIRE(result[1] == 2);
	REQUIRE(result[2] == 3);
}


TEMPLATE_LIST_TEST_CASE("Vector - Conj (complex)", "[Vector]",
						decltype(VectorCaseList<ScalarsComplex, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2if, 3 };
	const auto result = Conj(value);
	REQUIRE(result[0] == scalar_type_t<Vec>(-1));
	REQUIRE(result[1] == scalar_type_t<Vec>(-2if));
	REQUIRE(result[2] == scalar_type_t<Vec>(3));
}


TEMPLATE_LIST_TEST_CASE("Vector - Sum", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2, 3 };
	REQUIRE(Sum(value) == scalar_type_t<Vec>(4));
}


TEST_CASE("Vector - SumCompensated", "[Vector]") {
	const auto value = Vector(1.0f, 1e-12f, 1e-12f, 1e-12f, 1e-12f, 1.0f, 1e-12f);
	const auto [sum, compensation] = SumCompensated(value);
	REQUIRE(sum == 2.0f);
	REQUIRE(compensation == 5e-12f);
}


TEST_CASE("Vector - DotCompensated", "[Vector]") {
	const auto value = Vector(1.0f, 1e-6f, 1e-6f, 1e-6f, 1e-6f, 1.0f, 1e-6f);
	const auto [sum, compensation] = DotCompensated(value, value);
	REQUIRE(sum == 2.0f);
	REQUIRE(compensation == 5e-12f);
}


TEMPLATE_LIST_TEST_CASE("Vector - Dot (real)", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatAndInt32, PackingsAll>,
												VectorCaseList<ScalarsFloatAndInt32, PackingsAll>>{})) {
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;

	const VecLhs lhs = { -1, 2, 3 };
	const VecRhs rhs = { 4, 1, 2 };
	REQUIRE(Dot(lhs, rhs) == 4);
}


TEMPLATE_LIST_TEST_CASE("Vector - Dot (complex)", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsComplex32, PackingsAll>,
												VectorCaseList<ScalarsComplex32, PackingsAll>>{})) {
	using namespace std::complex_literals;
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;

	const VecLhs lhs = { -1, 2if, 3 };
	const VecRhs rhs = { 4, 1if, 2 };
	REQUIRE(std::complex<double>(Dot(lhs, rhs)) == std::complex<double>(4));
}


TEMPLATE_LIST_TEST_CASE("Vector - Length (real)", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatAndInt32, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2, 3 };
	REQUIRE(Length(value) == Catch::Approx(3.7416573867739));
}


TEMPLATE_LIST_TEST_CASE("Vector - Length (complex)", "[Vector]",
						decltype(VectorCaseList<ScalarsComplex, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { -1, 2if, 3 };
	static_assert(!is_complex_v<std::decay_t<decltype(Length(value))>>);
	REQUIRE(Length(value) == Catch::Approx(3.7416573867739));
}


TEMPLATE_LIST_TEST_CASE("Vector - LengthPrecise (real)", "[Vector]",
						decltype(VectorCaseList<ScalarsFloat32, PackingsAll>{})) {
	using Vec = typename TestType::template Vector<3>;

	SECTION("Underflow") {
		const Vec value = { -1e-20f, 2e-20f, 3e-20f };
		REQUIRE(LengthPrecise(value) == Catch::Approx(3.7416573867739e-20f));
	}
	SECTION("Overflow") {
		const Vec value = { -1e+20f, 2e+20f, 3e+20f };
		REQUIRE(LengthPrecise(value) == Catch::Approx(3.7416573867739e+20f));
	}
	SECTION("Zero") {
		const Vec value = { 0, 0, 0 };
		REQUIRE(LengthPrecise(value) == 0);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - LengthPrecise (complex)", "[Vector]",
						decltype(VectorCaseList<ScalarsComplex32, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	SECTION("Underflow") {
		const Vec value = { -1e-20f, 2e-20if, 3e-20f };
		static_assert(!is_complex_v<std::decay_t<decltype(Length(value))>>);
		REQUIRE(LengthPrecise(value) == Catch::Approx(3.7416573867739e-20f));
	}
	SECTION("Overflow") {
		const Vec value = { -1e+20f, 2e+20if, 3e+20f };
		REQUIRE(LengthPrecise(value) == Catch::Approx(3.7416573867739e+20f));
	}
	SECTION("Zero") {
		const Vec value = { 0, 0, 0 };
		REQUIRE(LengthPrecise(value) == 0);
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Distance", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloatAndInt32, PackingsAll>,
												VectorCaseList<ScalarsFloatAndInt32, PackingsAll>>{})) {
	using namespace std::complex_literals;
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;

	const VecLhs lhs = { 5, 3, 9 };
	const VecRhs rhs = { 8, 3, 5 };

	REQUIRE(Distance(lhs, rhs) == Catch::Approx(5));
}


TEMPLATE_LIST_TEST_CASE("Vector - DistancePrecise", "[Vector]",
						decltype(BinaryCaseList<VectorCaseList<ScalarsFloat32, PackingsAll>,
												VectorCaseList<ScalarsFloat32, PackingsAll>>{})) {
	using namespace std::complex_literals;
	using VecLhs = typename TestType::Lhs::template Vector<3>;
	using VecRhs = typename TestType::Rhs::template Vector<3>;

	SECTION("Underflow") {
		const VecLhs lhs = { 5e-20f, 3e-20f, 9e-20f };
		const VecRhs rhs = { 8e-20f, 3e-20f, 5e-20f };
		REQUIRE(DistancePrecise(lhs, rhs) == Catch::Approx(5e-20));
	}
	SECTION("Overflow") {
		const VecLhs lhs = { 5e+20f, 3e+20f, 9e+20f };
		const VecRhs rhs = { 8e+20f, 3e+20f, 5e+20f };
		REQUIRE(DistancePrecise(lhs, rhs) == Catch::Approx(5e+20));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Normalize", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatingAndComplex, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	const Vec value = { 2, 10, 11 };
	const auto normalized = Normalize(value);
	REQUIRE(std::real(normalized[0]) == Catch::Approx(2.0 / 15.0));
	REQUIRE(std::real(normalized[1]) == Catch::Approx(10.0 / 15.0));
	REQUIRE(std::real(normalized[2]) == Catch::Approx(11.0 / 15.0));
	REQUIRE(std::imag(normalized[0]) == Catch::Approx(0));
	REQUIRE(std::imag(normalized[1]) == Catch::Approx(0));
	REQUIRE(std::imag(normalized[2]) == Catch::Approx(0));
}


TEMPLATE_LIST_TEST_CASE("Vector - NormalizePrecise", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatingAndComplex, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;

	SECTION("Underflow") {
		const Vec value = { 2e-20f, 10e-20f, 11e-20f };
		const auto normalized = NormalizePrecise(value);
		REQUIRE(std::real(normalized[0]) == Catch::Approx(2.0 / 15.0));
		REQUIRE(std::real(normalized[1]) == Catch::Approx(10.0 / 15.0));
		REQUIRE(std::real(normalized[2]) == Catch::Approx(11.0 / 15.0));
		REQUIRE(std::imag(normalized[0]) == Catch::Approx(0));
		REQUIRE(std::imag(normalized[1]) == Catch::Approx(0));
		REQUIRE(std::imag(normalized[2]) == Catch::Approx(0));
	}
	SECTION("Overflow") {
		const Vec value = { 2e+20f, 10e+20f, 11e+20f };
		const auto normalized = NormalizePrecise(value);
		REQUIRE(std::real(normalized[0]) == Catch::Approx(2.0 / 15.0));
		REQUIRE(std::real(normalized[1]) == Catch::Approx(10.0 / 15.0));
		REQUIRE(std::real(normalized[2]) == Catch::Approx(11.0 / 15.0));
		REQUIRE(std::imag(normalized[0]) == Catch::Approx(0));
		REQUIRE(std::imag(normalized[1]) == Catch::Approx(0));
		REQUIRE(std::imag(normalized[2]) == Catch::Approx(0));
	}
	SECTION("Null") {
		const Vec value = { 0, 0, 0 };
		const auto normalized = NormalizePrecise(value);
		REQUIRE(std::real(normalized[0]) == Catch::Approx(1));
		REQUIRE(std::real(normalized[1]) == Catch::Approx(0));
		REQUIRE(std::real(normalized[2]) == Catch::Approx(0));
		REQUIRE(std::imag(normalized[0]) == Catch::Approx(0));
		REQUIRE(std::imag(normalized[1]) == Catch::Approx(0));
		REQUIRE(std::imag(normalized[2]) == Catch::Approx(0));
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Fill", "[Vector]",
						decltype(VectorCaseList<ScalarsAll, PackingsAll>{})) {
	using namespace std::complex_literals;
	using Vec = typename TestType::template Vector<3>;
	using T = scalar_type_t<Vec>;

	Vec value(static_cast<T>(0));
	REQUIRE(std::all_of(std::begin(value), std::end(value), [](const auto& v) { return v == static_cast<T>(0); }));
	Fill(value, 5);
	REQUIRE(std::all_of(std::begin(value), std::end(value), [](const auto& v) { return v == static_cast<T>(5); }));
}


TEMPLATE_LIST_TEST_CASE("Vector - Cross", "[Vector]",
						decltype(VectorCaseList<ScalarsFloating, PackingsAll>{})) {
	const auto triangleArea = [](auto a, auto b, auto c) {
		// Heron's formula
		const auto s = (a + b + c) / 2;
		return std::sqrt(s * (s - a) * (s - b) * (s - c));
	};


	SECTION("2D") {
		using Vec = typename TestType::template Vector<2>;
		const Vec v(1, 2);
		const auto result = Cross(v);
		REQUIRE(Dot(v, result) < 1e-6f);
		REQUIRE(Length(result) == Catch::Approx(Length(v)));

		// This validates the orientation of the 2D specialized and checks generalized values.
		const auto ilist = std::initializer_list{ v };
		const auto generalized = impl::CrossND(ilist.begin(), ilist.end());
		REQUIRE(generalized[0] == Catch::Approx(result[0]));
		REQUIRE(generalized[1] == Catch::Approx(result[1]));
	}
	SECTION("3D") {
		using Vec = typename TestType::template Vector<3>;
		const Vec a(0.2f, 1.0f, 0.3f);
		const Vec b(-1.0f, 0.05f, 0.2f);
		const auto result = Cross(a, b);
		const auto area = 2 * triangleArea(Length(a), Length(b), Length(a - b));
		REQUIRE((Dot(a, result) < 1e-6f && Dot(b, result) < 1e-6f)); // Check orthogonality.
		REQUIRE(area == Catch::Approx(Length(result))); // Check magnitude.
		REQUIRE(result.z > 0.8f); // Check orientation.

		// This validates the orientation of the 2D specialized and checks generalized values.
		const auto ilist = std::initializer_list{ a, b };
		const auto generalized = impl::CrossND(ilist.begin(), ilist.end());
		REQUIRE(generalized[0] == Catch::Approx(result[0]));
		REQUIRE(generalized[1] == Catch::Approx(result[1]));
		REQUIRE(generalized[2] == Catch::Approx(result[2]));
	}
	SECTION("4D") {
		using Vec = typename TestType::template Vector<4>;

		const Vec a(0.2f, 1.0f, 0.3f, 0.88f);
		const Vec b(-1.0f, 0.05f, 0.2f, 0.43f);
		const Vec c(0.22f, -0.15f, 0.37f, -1.0f);
		const auto result = Cross(a, b, c);

		REQUIRE((Dot(a, result) < 1e-6f && Dot(b, result) < 1e-6f && Dot(c, result) < 1e-6f)); // Check orthogonality.
	}
}


TEMPLATE_LIST_TEST_CASE("Vector - Gram-Schmidt", "[Vector]",
						decltype(VectorCaseList<ScalarsFloatingAndComplex, PackingsAll>{})) {

	SECTION("2D") {
		// This is the example from the Wikipedia article.
		using Vec = typename TestType::template Vector<2>;

		const std::array vectors = {
			Vec(3, 1),
			Vec(2, 2),
		};
		const std::array expected = {
			Vec(3, 1),
			Vec(-2.0 / 5, 6.0 / 5),
		};

		std::array<Vec, 2> result;
		GramSchmidtOrthogonalize(vectors.begin(), vectors.end(), result.begin());
		REQUIRE(result[0] == test_util::Approx(expected[0]));
		REQUIRE(result[1] == test_util::Approx(expected[1]));
		REQUIRE(std::abs(Dot(result[0], result[1])) < 1e-6f);
	}
	SECTION("4D") {
		// This is the example from the Wikipedia article.
		using Vec = typename TestType::template Vector<4>;

		const std::array vectors = {
			Vec(3, 1, 5, 2),
			Vec(2, 2, -1, 3),
			Vec(4, -7, -1, -2),
			Vec(1, 9, -1, -7),
		};

		std::array<Vec, 4> result;
		GramSchmidtOrthogonalize(vectors.begin(), vectors.end(), result.begin());

		for (auto it = result.begin(); it != result.end(); ++it) {
			REQUIRE(LengthPrecise(*it) > 0.1f); // No null vectors created.

			for (auto other = result.begin(); other != it; ++other) {
				REQUIRE(std::abs(Dot(NormalizePrecise(*it), NormalizePrecise(*other))) < 1e-6f);
			}
		}
	}
}