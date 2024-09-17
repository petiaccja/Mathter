// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Vector/Math.hpp"
#include "Quaternion.hpp"


namespace mathter {

/// <summary> Returns the angle of the rotation represented by quaternion. </summary>
/// <remarks> Only valid for unit quaternions. </remarks>
template <class T, eQuaternionLayout Layout, bool Packed>
T Angle(const Quaternion<T, Layout, Packed>& quat) {
	return T(2) * std::atan2(Length(Vector(quat.vector)), T(quat.scalar));
}


/// <summary> Returns the axis of rotation represented by quaternion. </summary>
/// <remarks> Only valid for unit quaternions. Returns (1,0,0) for near 180 degree rotations. </remarks>
template <class T, eQuaternionLayout Layout, bool Packed>
Vector<T, 3, Packed> Axis(const Quaternion<T, Layout, Packed>& quat) {
	return NormalizePrecise(Vector(quat.vector));
}


/// <summary> Returns the square of the absolute value. </summary>
/// <remarks> Just like complex numbers, it's the square of the length of the vector formed by the coefficients.
///			This is much faster than <see cref="Length">. </remarks>
template <class T, eQuaternionLayout Layout, bool Packed>
T LengthSquared(const Quaternion<T, Layout, Packed>& q) {
	return LengthSquared(Vector(q));
}


/// <summary> Returns the absolute value of the quaternion. </summary>
/// <remarks> Just like complex numbers, it's the length of the vector formed by the coefficients. </remarks>
template <class T, eQuaternionLayout Layout, bool Packed>
T Length(const Quaternion<T, Layout, Packed>& q) {
	return Length(Vector(q));
}


/// <summary> Returns the absolute value of the quaternion. </summary>
/// <remarks> Just like complex numbers, it's the length of the vector formed by the coefficients. </remarks>
template <class T, eQuaternionLayout Layout, bool Packed>
T LengthPrecise(const Quaternion<T, Layout, Packed>& q) {
	return LengthPrecise(Vector(q));
}


/// <summary> Returns the unit quaternion of the same direction. Does not change this object. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed> Normalize(const Quaternion<T, Layout, Packed>& q) {
	return Quaternion<T, Layout, Packed>{ Normalize(Vector(q)) };
}


/// <summary> Returns the unit quaternion of the same direction. Does not change this object. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed> NormalizePrecise(const Quaternion<T, Layout, Packed>& q) {
	return Quaternion<T, Layout, Packed>{ NormalizePrecise(Vector(q), Vector(Quaternion<T, Layout, Packed>(1, 0, 0, 0))) };
}


/// <summary> The euclidean length of the vector of the 4 elements of the quaternion. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
T Abs(const Quaternion<T, Layout, Packed>& q) {
	return LengthPrecise(Vector(q));
}


/// <summary> Negates the imaginary values of the quaternion. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed> Conj(const Quaternion<T, Layout, Packed>& q) {
	constexpr auto one = static_cast<T>(1);
	const Vector sign = { one, -one, -one, -one };
	return Quaternion<T, Layout, Packed>{ canonicalArg, q.canonical * sign };
}


/// <summary> Natural quaternion exponentiation, base e. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed> Exp(const Quaternion<T, Layout, Packed>& q) {
	const auto scalar = T(q.scalar);
	const auto vector = Vector(q.vector);
	const T vectorLength = LengthPrecise(vector);
	const T scalarExp = std::exp(scalar);
	const auto vectorNormalized = NormalizePrecise(vector);


	return scalarExp * Quaternion{ std::cos(vectorLength), vectorNormalized * std::sin(vectorLength) };
}


/// <summary> Natural quaternion logarithm, base e. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed> Log(const Quaternion<T, Layout, Packed>& q) {
	const auto scalar = T(q.scalar);
	const auto vector = Vector(q.vector);

	const auto vectorLength = LengthPrecise(vector);
	const auto quatLength = std::hypot(scalar, vectorLength);
	const auto vectorNormalized = NormalizePrecise(vector);

	return { std::log(quatLength), vectorNormalized * std::acos(scalar / quatLength) };
}


/// <summary> Raises <paramref name="q"/> to the power of <paramref name="a"/>. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed> Pow(const Quaternion<T, Layout, Packed>& q, T a) {
	return Exp(a * Log(q));
}


/// <summary> Returns the quaternion of opposite rotation. </summary>
template <class T, eQuaternionLayout Layout, bool Packed>
Quaternion<T, Layout, Packed> Inverse(const Quaternion<T, Layout, Packed>& q) {
	const auto length = LengthSquared(q);
	assert(length != T(0) && "Zero Quaternion is not invertable.");
	return Conj(q) / LengthSquared(q);
}

} // namespace mathter