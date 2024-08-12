#pragma once

#include <Mathter/Common/TypeTraits.hpp>
#include <Mathter/Matrix/Matrix.hpp>


namespace test_util {

template <class Mat, class Vec, class = std::enable_if_t<mathter::is_matrix_v<Mat> && (mathter::is_vector_v<Vec> || mathter::is_matrix_v<Vec>)>>
auto ApplyTransform(const Mat& transform, const Vec& vector) {
	if constexpr (mathter::order_v<Mat> == mathter::eMatrixOrder::PRECEDE_VECTOR) {
		return transform * vector;
	}
	else {
		return vector * transform;
	}
}

} // namespace test_util