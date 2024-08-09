#pragma once

#include <Mathter/Matrix/Matrix.hpp>


namespace test_util {

template <class Mat>
Mat ZeroLowerTriangle(const Mat& m) {
	auto copy = m;
	for (size_t row = 0; row < row_count_v<Mat>; ++row) {
		for (size_t col = 0; col <= row && col < column_count_v<Mat>; ++col) {
			copy(row, col) = scalar_type_t<Mat>(0);
		}
	}
	return copy;
}


template <class Mat>
Mat ZeroUpperTriangle(const Mat& m) {
	auto copy = m;
	for (size_t row = 0; row < row_count_v<Mat>; ++row) {
		for (size_t col = row; col < column_count_v<Mat>; ++col) {
			copy(row, col) = scalar_type_t<Mat>(0);
		}
	}
	return copy;
}

} // namespace test_util