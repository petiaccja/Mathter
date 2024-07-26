// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include "../Cases.hpp"

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


// TODO: Tests for cross product!
// Can't do until matrices are also fixed.