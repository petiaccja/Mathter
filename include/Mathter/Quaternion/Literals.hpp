// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Quaternion.hpp"

namespace mathter {

namespace quat_literals {

	inline auto operator"" _i(unsigned long long int arg) {
		return Quaternion<double>(0.0, double(arg), 0.0, 0.0);
	}
	inline auto operator"" _j(unsigned long long int arg) {
		return Quaternion<double>(0.0, 0.0, double(arg), 0.0);
	}
	inline auto operator"" _k(unsigned long long int arg) {
		return Quaternion<double>(0.0, 0.0, 0.0, double(arg));
	}

	inline auto operator"" _i(long double arg) {
		return Quaternion<double>(0.0, double(arg), 0.0, 0.0);
	}
	inline auto operator"" _j(long double arg) {
		return Quaternion<double>(0.0, 0.0, double(arg), 0.0);
	}
	inline auto operator"" _k(long double arg) {
		return Quaternion<double>(0.0, 0.0, 0.0, double(arg));
	}


	inline auto operator"" _if(unsigned long long int arg) {
		return Quaternion<float>(0.0f, float(arg), 0.0f, 0.0f);
	}
	inline auto operator"" _jf(unsigned long long int arg) {
		return Quaternion<float>(0.0f, 0.0f, float(arg), 0.0f);
	}
	inline auto operator"" _kf(unsigned long long int arg) {
		return Quaternion<float>(0.0f, 0.0f, 0.0f, float(arg));
	}

	inline auto operator"" _if(long double arg) {
		return Quaternion<float>(0.0f, float(arg), 0.0f, 0.0f);
	}
	inline auto operator"" _jf(long double arg) {
		return Quaternion<float>(0.0f, 0.0f, float(arg), 0.0f);
	}
	inline auto operator"" _kf(long double arg) {
		return Quaternion<float>(0.0f, 0.0f, 0.0f, float(arg));
	}


	inline auto operator"" _il(unsigned long long int arg) {
		return Quaternion<long double>(0.0l, static_cast<long double>(arg), 0.0l, 0.0l);
	}
	inline auto operator"" _jl(unsigned long long int arg) {
		return Quaternion<long double>(0.0l, 0.0l, static_cast<long double>(arg), 0.0l);
	}
	inline auto operator"" _kl(unsigned long long int arg) {
		return Quaternion<long double>(0.0l, 0.0l, 0.0l, static_cast<long double>(arg));
	}

	inline auto operator"" _il(long double arg) {
		return Quaternion<long double>(0.0l, arg, 0.0l, 0.0l);
	}
	inline auto operator"" _jl(long double arg) {
		return Quaternion<long double>(0.0l, 0.0l, arg, 0.0l);
	}
	inline auto operator"" _kl(long double arg) {
		return Quaternion<long double>(0.0l, 0.0l, 0.0l, arg);
	}

} // namespace quat_literals

} // namespace mathter