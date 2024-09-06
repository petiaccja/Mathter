#pragma once

#include <Mathter/Matrix.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>


template <class Mat, class Tol>
void VerifyUnitary(const Mat& m, Tol tolerance) {
	using namespace mathter;

	if constexpr (row_count_v<Mat> > column_count_v<Mat>) {
		const auto iden = ConjTranspose(m) * m;
		REQUIRE(iden == test_util::Approx(decltype(iden)(Identity()), tolerance));
	}
	else if constexpr (row_count_v<Mat> < column_count_v<Mat>) {
		const auto iden = m * ConjTranspose(m);
		REQUIRE(iden == test_util::Approx(decltype(iden)(Identity()), tolerance));
	}
	else {
		const auto iden1 = ConjTranspose(m) * m;
		const auto iden2 = m * ConjTranspose(m);
		const auto det = std::abs(Determinant(m));
		REQUIRE_THAT(det, Catch::Matchers::WithinAbs(1.0, tolerance));
		REQUIRE(iden1 == test_util::Approx(decltype(iden1)(Identity()), tolerance));
		REQUIRE(iden2 == test_util::Approx(decltype(iden2)(Identity()), tolerance));
	}
}


template <class Mat, class MatPI, class Tol>
void VerifyPseudoinverse(const Mat& m, const MatPI& pi, Tol tolerance) {
	using namespace mathter;

	static_assert(row_count_v<std::decay_t<decltype(pi)>> == column_count_v<std::decay_t<decltype(m)>>);
	static_assert(column_count_v<std::decay_t<decltype(pi)>> == row_count_v<std::decay_t<decltype(m)>>);
	// The Moore-Penrose conditions
	REQUIRE(m * pi * m == test_util::Approx(m, tolerance));
	REQUIRE(pi * m * pi == test_util::Approx(pi, tolerance));
	REQUIRE(ConjTranspose(m * pi) == test_util::Approx(m * pi, tolerance));
	REQUIRE(ConjTranspose(pi * m) == test_util::Approx(pi * m, tolerance));
}