#pragma once

#include <Mathter/Common/OptimizationUtil.hpp>
#include <Mathter/Common/TypeTraits.hpp>
#include <Mathter/Transforms.hpp>

#include <array>
#include <complex>
#include <random>


template <class Obj, int Count>
std::array<Obj, Count> MakeRandomInput(
	mathter::remove_complex_t<mathter::scalar_type_t<Obj>> min = 0,
	mathter::remove_complex_t<mathter::scalar_type_t<Obj>> max = 1) {
	using Scalar = mathter::scalar_type_t<Obj>;
	using Real = mathter::remove_complex_t<Scalar>;
	using namespace std::complex_literals;

	std::mt19937_64 rne;
	std::uniform_real_distribution<Real> rng(Real(-0.5), Real(0.5));

	std::array<Obj, Count> r;
	for (auto& v : r) {
		if constexpr (mathter::is_scalar_v<Obj>) {
			if constexpr (mathter::is_complex_v<Scalar>) {
				const auto real = rng(rne);
				const auto imag = rng(rne);
				v = real + Scalar(1if) * imag;
			}
			else {
				v = rng(rne);
			}
		}
		else {
			if constexpr (mathter::is_complex_v<Scalar>) {
				const Obj real = mathter::Random(rng, rne);
				const Obj imag = mathter::Random(rng, rne);
				v = real + Scalar(1if) * imag;
			}
			else {
				v = mathter::Random(rng, rne);
			}
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
