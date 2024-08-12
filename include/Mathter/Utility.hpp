// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Common/TypeTraits.hpp"
#include "Common/Types.hpp"
#include "Vector/Math.hpp"
#include "Vector/Vector.hpp"

#include <algorithm>
#include <type_traits>


namespace mathter {

//------------------------------------------------------------------------------
// Public utility stuff
//------------------------------------------------------------------------------

// Common mathematical constants

/// <summary> Accurate mathematical constants. </summary>
template <class Scalar>
class Constants {
public:
	static constexpr auto Pi = Scalar(3.1415926535897932384626433832795028841971693993751);
	static constexpr auto PiHalf = Scalar(1.5707963267948966192313216916397514420985846996876);
	static constexpr auto PiFourth = Scalar(0.78539816339744830961566084581987572104929234984378);
	static constexpr auto E = Scalar(2.7182818284590452353602874713526624977572470937);
	static constexpr auto Sqrt2 = Scalar(1.4142135623730950488016887242096980785696718753769);
	static constexpr auto SqrtHalf = Scalar(0.70710678118654752440084436210484903928483593768847);
};


// Radians and degrees

/// <summary> Converts radians to degrees. </summary>
template <class Scalar>
auto Rad2Deg(Scalar rad) {
	using Real = remove_complex_t<Scalar>;
	using ComputeT = std::conditional_t<std::is_floating_point_v<Real>, Real, long double>;
	return rad / Constants<ComputeT>::Pi * ComputeT(180);
}

/// <summary> Converts degrees to radians. </summary>
template <class Scalar>
auto Deg2Rad(Scalar deg) {
	using Real = remove_complex_t<Scalar>;
	using ComputeT = std::conditional_t<std::is_floating_point_v<Real>, Real, long double>;
	return deg / ComputeT(180) * Constants<ComputeT>::Pi;
}


// Clamp and saturate

/// <summary> Limits arg to the range [lower, upper], making it either lower or upper if out of range. </summary>
template <class Scalar>
Scalar Clamp(Scalar arg, Scalar lower, Scalar upper) {
	return std::clamp(arg, lower, upper);
}

/// <summary> Clamps all elements of the vector according to the scalar clamp. </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Clamp(const Vector<T, Dim, Packed>& arg, T lower, T upper);

/// <summary> Clamps argument into range [0,1]. </summary>
template <class Scalar>
Scalar Saturate(Scalar arg) {
	return Clamp(arg, Scalar(0), Scalar(1));
}

/// <summary> Clamps all elements into range [0,1]. </summary>
template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Saturate(const Vector<T, Dim, Packed>& arg);


template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Clamp(const Vector<T, Dim, Packed>& arg, T lower, T upper) {
	using Vec = Vector<T, Dim, Packed>;
	return Min(Max(arg, Vec(lower)), Vec(upper));
}

template <class T, int Dim, bool Packed>
Vector<T, Dim, Packed> Saturate(const Vector<T, Dim, Packed>& arg) {
	return Clamp(arg, T(0), T(1));
}

} // namespace mathter