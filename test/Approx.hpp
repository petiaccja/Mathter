#pragma once

#include <Mathter/Vector/Vector.hpp>

#include <utility>


namespace test_util {

template <class T>
struct Approx;


template <class T, int Dim, bool Packed>
struct Approx<mathter::Vector<T, Dim, Packed>> {
	Approx(const mathter::Vector<T, Dim, Packed>& object, T tolerance) : object(object), tolerance(tolerance) {}


	mathter::Vector<T, Dim, Packed> object;
	T tolerance;
};


template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator==(const Approx<mathter::Vector<T1, Dim, Packed1>>& lhs, const Approx<mathter::Vector<T1, Dim, Packed1>>& rhs) {
	return LengthPrecise(lhs.object - rhs.object) < std::min(lhs.tolerance, rhs.tolerance) * LengthPrecise(lhs);
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator==(mathter::Vector<T1, Dim, Packed1>& lhs, const Approx<mathter::Vector<T1, Dim, Packed1>>& rhs) {
	return LengthPrecise(lhs - rhs.object) < rhs.tolerance * LengthPrecise(lhs);
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator==(const Approx<mathter::Vector<T1, Dim, Packed1>>& lhs, mathter::Vector<T1, Dim, Packed1>& rhs) {
	return LengthPrecise(lhs.object - rhs) < lhs.tolerance * LengthPrecise(rhs);
}

} // namespace test_util