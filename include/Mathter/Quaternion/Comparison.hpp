// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "Quaternion.hpp"

namespace mathter {

template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
bool operator==(const Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return lhs.canonical == rhs.canonical;
}


template <class T1, eQuaternionLayout Layout1, bool Packed1,
		  class T2, eQuaternionLayout Layout2, bool Packed2>
bool operator!=(const Quaternion<T1, Layout1, Packed1>& lhs, const Quaternion<T2, Layout2, Packed2>& rhs) {
	return !(lhs == rhs);
}

} // namespace mathter