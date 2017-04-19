#pragma once



// Deprecated
union Simd4f {
	float v[4];

	static inline Simd4f mul(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] * rhs.v[0];
		res.v[1] = lhs.v[1] * rhs.v[1];
		res.v[2] = lhs.v[2] * rhs.v[2];
		res.v[3] = lhs.v[3] * rhs.v[3];
		return res;
	}

	static inline Simd4f div(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] / rhs.v[0];
		res.v[1] = lhs.v[1] / rhs.v[1];
		res.v[2] = lhs.v[2] / rhs.v[2];
		res.v[3] = lhs.v[3] / rhs.v[3];
		return res;
	}

	static inline Simd4f add(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] + rhs.v[0];
		res.v[1] = lhs.v[1] + rhs.v[1];
		res.v[2] = lhs.v[2] + rhs.v[2];
		res.v[3] = lhs.v[3] + rhs.v[3];
		return res;
	}

	static inline Simd4f sub(const Simd4f& lhs, const Simd4f& rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] - rhs.v[0];
		res.v[1] = lhs.v[1] - rhs.v[1];
		res.v[2] = lhs.v[2] - rhs.v[2];
		res.v[3] = lhs.v[3] - rhs.v[3];
		return res;
	}

	static inline Simd4f mul(const Simd4f& lhs, float rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] * rhs;
		res.v[1] = lhs.v[1] * rhs;
		res.v[2] = lhs.v[2] * rhs;
		res.v[3] = lhs.v[3] * rhs;
		return res;
	}

	static inline Simd4f div(const Simd4f& lhs, float rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] / rhs;
		res.v[1] = lhs.v[1] / rhs;
		res.v[2] = lhs.v[2] / rhs;
		res.v[3] = lhs.v[3] / rhs;
		return res;
	}

	static inline Simd4f add(const Simd4f& lhs, float rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] + rhs;
		res.v[1] = lhs.v[1] + rhs;
		res.v[2] = lhs.v[2] + rhs;
		res.v[3] = lhs.v[3] + rhs;
		return res;
	}

	static inline Simd4f sub(const Simd4f& lhs, float rhs) {
		Simd4f res;
		res.v[0] = lhs.v[0] - rhs;
		res.v[1] = lhs.v[1] - rhs;
		res.v[2] = lhs.v[2] - rhs;
		res.v[3] = lhs.v[3] - rhs;
		return res;
	}

	static inline Simd4f mad(const Simd4f& a, const Simd4f& b, const Simd4f& c) {
		return add(mul(a, b), c);
	}

	static inline Simd4f spread(float value) {
		Simd4f res;
		res.v[0] = value;
		res.v[1] = value;
		res.v[2] = value;
		res.v[3] = value;
		return res;
	}

	static inline Simd4f set(float x, float y, float z, float w) {
		Simd4f res;
		res.v[0] = x;
		res.v[1] = y;
		res.v[2] = z;
		res.v[3] = w;
		return res;
	}

	static inline float dot(const Simd4f& lhs, const Simd4f& rhs) {
		float sum = lhs.v[0] * rhs.v[0];
		sum += lhs.v[1] * rhs.v[1];
		sum += lhs.v[2] * rhs.v[2];
		sum += lhs.v[3] * rhs.v[3];
		return sum;
	}
};