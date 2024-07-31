#pragma once

#include <Mathter/Vector/Math.hpp>
#include <Mathter/Vector/Vector.hpp>

#include <cmath>


namespace test_util {

template <class T, bool Packed>
mathter::Vector<T, 3, Packed> Rotate(const mathter::Vector<T, 3, Packed>& vector, const mathter::Vector<T, 3, Packed>& axis, T angle) {
	const auto basisS = Cross(axis, vector);
	const auto basisC = Cross(basisS, axis);
	const auto basisAxis = axis * Dot(axis, vector);
	return basisAxis + basisS * std::sin(angle) + basisC * std::cos(angle);
}

} // namespace test_util