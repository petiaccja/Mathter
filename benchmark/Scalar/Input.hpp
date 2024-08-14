#pragma once

#include <Mathter/Common/OptimizationUtil.hpp>
#include <Mathter/Common/TypeTraits.hpp>

#include <array>
#include <complex>
#include <random>


using namespace mathter;


template <class Scalar, int Count>
std::array<Scalar, Count> MakeRandomInput(remove_complex_t<Scalar> min = remove_complex_t<Scalar>(0.9f),
										  remove_complex_t<Scalar> max = remove_complex_t<Scalar>(0.98f)) {
	using Real = remove_complex_t<Scalar>;
	using namespace std::complex_literals;

	std::mt19937_64 rne;
	std::uniform_real_distribution<Real> rng(min, max);

	std::array<Scalar, Count> r;
	for (auto& v : r) {
		if constexpr (is_complex_v<Scalar>) {
			const Scalar real = rng(rne);
			const Scalar imag = rng(rne);
			v = real + Scalar(1if) * imag;
		}
		else {
			v = rng(rne);
		}
	}
	return r;
};


template <class Scalar, int Count>
std::array<Scalar, Count> MakeConstantInput(Scalar value) {
	std::array<Scalar, Count> r;
	for (auto& v : r) {
		v = value;
	}
	return r;
};
