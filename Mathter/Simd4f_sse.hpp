#pragma once

#include <emmintrin.h>


union alignas(16) Simd4f {
	__m128 reg;
	float v[4];


	static inline Simd4f mul(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.reg = _mm_mul_ps(lhs.reg, rhs.reg);
		return res;
	}

	static inline Simd4f div(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.reg = _mm_div_ps(lhs.reg, rhs.reg);
		return res;
	}

	static inline Simd4f add(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.reg = _mm_add_ps(lhs.reg, rhs.reg);
		return res;
	}

	static inline Simd4f sub(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.reg = _mm_sub_ps(lhs.reg, rhs.reg);
		return res;
	}

	static inline Simd4f mul(const Simd4f& lhs, float rhs) {
		Simd4f res;
		__m128 tmp = _mm_set1_ps(rhs);
		res.reg = _mm_mul_ps(lhs.reg, tmp);
		return res;
	}

	static inline Simd4f div(const Simd4f& lhs, float rhs) {
		Simd4f res;
		__m128 tmp = _mm_set1_ps(rhs);
		res.reg = _mm_div_ps(lhs.reg, tmp);
		return res;
	}

	static inline Simd4f add(const Simd4f& lhs, float rhs) {
		Simd4f res;
		__m128 tmp = _mm_set1_ps(rhs);
		res.reg = _mm_add_ps(lhs.reg, tmp);
		return res;
	}

	static inline Simd4f sub(const Simd4f& lhs, float rhs) {
		Simd4f res;
		__m128 tmp = _mm_set1_ps(rhs);
		res.reg = _mm_sub_ps(lhs.reg, tmp);
		return res;
	}

	static inline Simd4f mad(const Simd4f& a, const Simd4f& b, const Simd4f& c) {
		return add(mul(a, b), c);
	}

	static inline Simd4f spread(float value) {
		Simd4f res;
		res.reg = _mm_set1_ps(value);
		return res;
	}

	static inline Simd4f set(float x, float y, float z, float w) {
		Simd4f res;
		res.reg = _mm_setr_ps(x, y, z, w);
		return res;
	}

	static inline float dot(const Simd4f& lhs, const Simd4f& rhs) {
		float sum;
		Simd4f m = mul(lhs, rhs);
		sum = m.v[0] + m.v[1] + m.v[2] + m.v[3];
		return sum;
	}
};