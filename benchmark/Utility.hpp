#pragma once

#include <Mathter/Matrix.hpp>
#include <Mathter/Vector.hpp>

#include <array>


using namespace mathter;


using Vec2f = Vector<float, 2, false>;
using Vec3f = Vector<float, 3, false>;
using Vec4f = Vector<float, 4, false>;

using Vec2fp = Vector<float, 2, true>;
using Vec3fp = Vector<float, 3, true>;
using Vec4fp = Vector<float, 4, true>;

using Vec2d = Vector<double, 2, false>;
using Vec3d = Vector<double, 3, false>;
using Vec4d = Vector<double, 4, false>;

using Vec2dp = Vector<double, 2, true>;
using Vec3dp = Vector<double, 3, true>;
using Vec4dp = Vector<double, 4, true>;


using Mat22f_fc = Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, false>;
using Mat33f_fc = Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, false>;
using Mat44f_fc = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, false>;
using Mat22f_fr = Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat33f_fr = Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat44f_fr = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;

using Mat22fp_fc = Matrix<float, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;
using Mat33fp_fc = Matrix<float, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;
using Mat44fp_fc = Matrix<float, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;

using Mat22d_fc = Matrix<double, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, false>;
using Mat33d_fc = Matrix<double, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, false>;
using Mat44d_fc = Matrix<double, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, false>;
using Mat22d_fr = Matrix<double, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat33d_fr = Matrix<double, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;
using Mat44d_fr = Matrix<double, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, false>;

using Mat22dp_fc = Matrix<double, 2, 2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;
using Mat33dp_fc = Matrix<double, 3, 3, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;
using Mat44dp_fc = Matrix<double, 4, 4, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, true>;


template <class VecT, int N>
auto MakeVectorArray() {
	std::array<VecT, N> items;
	for (auto& v : items) {
		Fill(v, 1);
	}
	return items;
}

template <class VecT, int N>
auto MakeMatrixArray() {
	std::array<VecT, N> items;
	for (auto& v : items) {
		v = Identity();
	}
	return items;
}

template <class... Args, size_t N>
auto TuplizeArrays(const std::array<Args, N>&... arrays) {
	std::array<std::tuple<Args...>, N> items;
	for (size_t i = 0; i < N; ++i) {
		items[i] = { arrays[i]... };
	}
	return items;
}


inline auto opAdd = [](const auto& v) { return std::get<0>(v) + std::get<1>(v); };
inline auto opMul = [](const auto& v) { return std::get<0>(v) * std::get<1>(v); };
inline auto opDiv = [](const auto& v) { return std::get<0>(v) / std::get<1>(v); };
inline auto feedBinary = [](const auto& result, const auto& init) { return std::tuple{ result, std::get<0>(init) }; };