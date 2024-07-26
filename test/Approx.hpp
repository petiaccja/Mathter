#pragma once

#include <Mathter/Common/TypeTraits.hpp>
#include <Mathter/Matrix/Comparison.hpp>
#include <Mathter/Matrix/Matrix.hpp>
#include <Mathter/Vector/Comparison.hpp>
#include <Mathter/Vector/Vector.hpp>

#include <utility>


namespace test_util {

template <class T>
struct Approx;


template <class T, class... Args>
Approx(T&&, Args&&...) -> Approx<std::decay_t<T>>;


template <class T>
constexpr auto DefaultTolerance() {
	return static_cast<remove_complex_t<T>>(10) * std::numeric_limits<mathter::remove_complex_t<std::decay_t<T>>>::epsilon();
}


template <class T, int Dim, bool Packed>
struct Approx<mathter::Vector<T, Dim, Packed>> {
	Approx(const mathter::Vector<T, Dim, Packed>& object, mathter::remove_complex_t<T> tolerance = DefaultTolerance<T>())
		: object(object), tolerance(tolerance) {}


	mathter::Vector<T, Dim, Packed> object;
	mathter::remove_complex_t<T> tolerance;
};


template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator==(const Approx<mathter::Vector<T1, Dim, Packed1>>& lhs, const Approx<mathter::Vector<T2, Dim, Packed2>>& rhs) {
	return LengthPrecise(lhs.object - rhs.object) < std::min(lhs.tolerance, rhs.tolerance) * LengthPrecise(lhs);
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator==(mathter::Vector<T1, Dim, Packed1>& lhs, const Approx<mathter::Vector<T2, Dim, Packed2>>& rhs) {
	return LengthPrecise(lhs - rhs.object) < rhs.tolerance * LengthPrecise(lhs);
}


template <class T1, class T2, int Dim, bool Packed1, bool Packed2>
bool operator==(const Approx<mathter::Vector<T1, Dim, Packed1>>& lhs, mathter::Vector<T2, Dim, Packed2>& rhs) {
	return LengthPrecise(lhs.object - rhs) < lhs.tolerance * LengthPrecise(rhs);
}


template <class T, int Rows, int Columns, mathter::eMatrixOrder Order, mathter::eMatrixLayout Layout, bool Packed>
struct Approx<mathter::Matrix<T, Rows, Columns, Order, Layout, Packed>> {
	Approx(const mathter::Matrix<T, Rows, Columns, Order, Layout, Packed>& object, mathter::remove_complex_t<T> tolerance = DefaultTolerance<T>())
		: object(object), tolerance(tolerance) {}


	mathter::Matrix<T, Rows, Columns, Order, Layout, Packed> object;
	mathter::remove_complex_t<T> tolerance;
};


template <class T1, mathter::eMatrixOrder Order1, mathter::eMatrixLayout Layout1, bool Packed1,
		  class T2, mathter::eMatrixOrder Order2, mathter::eMatrixLayout Layout2, bool Packed2,
		  int Rows, int Columns>
bool operator==(const Approx<mathter::Matrix<T1, Rows, Columns, Order1, Layout1, Packed1>>& lhs, const Approx<mathter::Matrix<T2, Rows, Columns, Order2, Layout2, Packed2>>& rhs) {
	return NormPrecise(lhs.object - rhs.object) < std::min(lhs.tolerance, rhs.tolerance) * NormPrecise(lhs);
}


template <class T1, mathter::eMatrixOrder Order1, mathter::eMatrixLayout Layout1, bool Packed1,
		  class T2, mathter::eMatrixOrder Order2, mathter::eMatrixLayout Layout2, bool Packed2,
		  int Rows, int Columns>
bool operator==(const mathter::Matrix<T1, Rows, Columns, Order1, Layout1, Packed1>& lhs, const Approx<mathter::Matrix<T2, Rows, Columns, Order2, Layout2, Packed2>>& rhs) {
	return NormPrecise(lhs - rhs.object) < rhs.tolerance * NormPrecise(lhs);
}


template <class T1, mathter::eMatrixOrder Order1, mathter::eMatrixLayout Layout1, bool Packed1,
		  class T2, mathter::eMatrixOrder Order2, mathter::eMatrixLayout Layout2, bool Packed2,
		  int Rows, int Columns>
bool operator==(const Approx<mathter::Matrix<T1, Rows, Columns, Order1, Layout1, Packed1>>& lhs, const mathter::Matrix<T2, Rows, Columns, Order2, Layout2, Packed2>& rhs) {
	return NormPrecise(lhs.object - rhs) < lhs.tolerance * NormPrecise(rhs);
}

} // namespace test_util