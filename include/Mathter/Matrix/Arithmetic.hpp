// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma once

#include "../Common/LoopUtil.hpp"
#include "../Vector/Arithmetic.hpp"
#include "../Vector/Math.hpp"
#include "Matrix.hpp"

namespace mathter {

namespace impl {

	template <class T1, bool Packed1,
			  class T2, bool Packed2,
			  int Rows1, int Match, int Columns2, eMatrixOrder Order>
	auto Multiply(const Matrix<T1, Rows1, Match, Order, eMatrixLayout::ROW_MAJOR, Packed1>& lhs,
				  const Matrix<T2, Match, Columns2, Order, eMatrixLayout::ROW_MAJOR, Packed2>& rhs) {
		using T = common_arithmetic_type_t<T1, T2>;
		using M = Matrix<T, Rows1, Columns2, Order, eMatrixLayout::ROW_MAJOR, Packed1 && Packed2>;

		const auto calcPartialRow = [&lhs, &rhs](size_t rowIdx, size_t runIdx) {
			return lhs(rowIdx, runIdx) * rhs.Row(runIdx);
		};

		const auto calcRow = [&calcPartialRow](size_t rowIdx) {
			return ::mathter::LoopUnroll<Match>([rowIdx, &calcPartialRow](auto... runIdx) {
				return (... + calcPartialRow(rowIdx, runIdx));
			});
		};

		return ::mathter::LoopUnroll<Rows1>([&calcRow](auto... rowIdx) {
			return M(stripeArg, calcRow(rowIdx)...);
		});
	}


	template <class T1, bool Packed1,
			  class T2, bool Packed2,
			  int Rows1, int Match, int Columns2, eMatrixOrder Order>
	auto Multiply(const Matrix<T1, Rows1, Match, Order, eMatrixLayout::ROW_MAJOR, Packed1>& lhs,
				  const Matrix<T2, Match, Columns2, Order, eMatrixLayout::COLUMN_MAJOR, Packed2>& rhs) {
		using T = common_arithmetic_type_t<T1, T2>;
		using M = Matrix<T, Rows1, Columns2, Order, eMatrixLayout::ROW_MAJOR, Packed1 && Packed2>;

		const auto calcElement = [&lhs, &rhs](size_t rowIdx, size_t colIdx) {
			return Sum(lhs.Row(rowIdx) * rhs.Column(colIdx));
		};

		const auto calcRow = [&calcElement](size_t rowIdx) {
			return ::mathter::LoopUnroll<Columns2>([rowIdx, &calcElement](auto... runIdx) {
				return Vector(calcElement(rowIdx, runIdx)...);
			});
		};

		return ::mathter::LoopUnroll<Rows1>([&calcRow](auto... rowIdx) {
			return M(stripeArg, calcRow(rowIdx)...);
		});
	}


	template <class T1, bool Packed1,
			  class T2, eMatrixLayout Layout2, bool Packed2,
			  int Rows1, int Match, int Columns2, eMatrixOrder Order>
	auto Multiply(const Matrix<T1, Rows1, Match, Order, eMatrixLayout::COLUMN_MAJOR, Packed1>& lhs,
				  const Matrix<T2, Match, Columns2, Order, Layout2, Packed2>& rhs) {
		using T = common_arithmetic_type_t<T1, T2>;
		using M = Matrix<T, Rows1, Columns2, Order, eMatrixLayout::COLUMN_MAJOR, Packed1 && Packed2>;

		const auto calcPartialColumn = [&lhs, &rhs](size_t colIdx, size_t runIdx) {
			return lhs.Column(runIdx) * rhs(runIdx, colIdx);
		};

		const auto calcColumn = [&calcPartialColumn](size_t colIdx) {
			return ::mathter::LoopUnroll<Match>([colIdx, &calcPartialColumn](auto... runIdx) {
				return (... + calcPartialColumn(colIdx, runIdx));
			});
		};

		return ::mathter::LoopUnroll<Columns2>([&calcColumn](auto... colIdx) {
			return M(stripeArg, calcColumn(colIdx)...);
		});
	}


	template <class T1, eMatrixLayout Layout1, bool Packed1,
			  class T2, eMatrixLayout Layout2, bool Packed2,
			  int Rows, int Columns, eMatrixOrder Order,
			  class Func>
	auto Elementwise(const Matrix<T1, Rows, Columns, Order, Layout1, Packed1>& lhs,
					 const Matrix<T2, Rows, Columns, Order, Layout2, Packed2>& rhs,
					 Func&& func) {
		using T = std::invoke_result_t<Func, T1, T2>;
		using M = Matrix<T, Rows, Columns, Order, Layout1, Packed1 && Packed2>;

		const auto calcStripe = [&lhs, &rhs, &func](size_t stripeIdx) {
			if constexpr (Layout1 == eMatrixLayout::ROW_MAJOR) {
				return func(lhs.Row(stripeIdx), rhs.Row(stripeIdx));
			}
			else {
				return func(lhs.Column(stripeIdx), rhs.Column(stripeIdx));
			}
		};

		return ::mathter::LoopUnroll<M::stripeCount>([&calcStripe](auto... stripeIdx) {
			return M(stripeArg, calcStripe(stripeIdx)...);
		});
	}


	template <class T1, int Rows1, int Columns1, eMatrixOrder Order1, eMatrixLayout Layout1, bool Packed1,
			  class T2,
			  class Func,
			  class = std::enable_if_t<is_scalar_v<T2>>>
	auto Scalar(const Matrix<T1, Rows1, Columns1, Order1, Layout1, Packed1>& lhs,
				const T2& rhs,
				Func&& func) {
		using T = std::invoke_result_t<Func, T1, T2>;
		using M = Matrix<T, Rows1, Columns1, Order1, Layout1, Packed1>;

		return ::mathter::LoopUnroll<M::stripeCount>([&lhs, &rhs, &func](auto... stripeIdx) {
			return M(stripeArg, func(lhs.stripes[stripeIdx], rhs)...);
		});
	}


	template <class T1,
			  class T2, int Rows2, int Columns2, eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2,
			  class Func,
			  class = std::enable_if_t<is_scalar_v<T1>>>
	auto Scalar(const T1& lhs,
				const Matrix<T2, Rows2, Columns2, Order2, Layout2, Packed2>& rhs,
				Func&& func) {
		using T = std::invoke_result_t<Func, T1, T2>;
		using M = Matrix<T, Rows2, Columns2, Order2, Layout2, Packed2>;

		return ::mathter::LoopUnroll<M::stripeCount>([&lhs, &rhs, &func](auto... stripeIdx) {
			return M(stripeArg, func(lhs, rhs.stripes[stripeIdx])...);
		});
	}

	template <class T1, bool Packed1,
			  class T2, int Columns2, bool Packed2,
			  int Match>
	auto Multiply(const Vector<T1, Match, Packed1>& lhs,
				  const Matrix<T2, Match, Columns2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::ROW_MAJOR, Packed2>& rhs) {
		return ::mathter::LoopUnroll<Match>([lhs, rhs](auto... runIdx) {
			return (... + (lhs(runIdx) * rhs.Row(runIdx)));
		});
	}

	template <class T1, bool Packed1,
			  class T2, int Columns2, bool Packed2,
			  int Match>
	auto Multiply(const Vector<T1, Match, Packed1>& lhs,
				  const Matrix<T2, Match, Columns2, eMatrixOrder::FOLLOW_VECTOR, eMatrixLayout::COLUMN_MAJOR, Packed2>& rhs) {
		return ::mathter::LoopUnroll<Columns2>([lhs, rhs](auto... elementIdx) {
			return Vector(Sum(lhs * rhs.Column(elementIdx))...);
		});
	}

	template <class T1, int Rows1, bool Packed1,
			  class T2, bool Packed2,
			  int Match>
	auto Multiply(const Matrix<T1, Rows1, Match, eMatrixOrder::PRECEDE_VECTOR, eMatrixLayout::COLUMN_MAJOR, Packed1>& lhs,
				  const Vector<T2, Match, Packed2>& rhs) {
		return ::mathter::LoopUnroll<Match>([lhs, rhs](auto... runIdx) {
			return (... + (lhs.Column(runIdx) * rhs(runIdx)));
		});
	}

	template <class T1, int Rows1, bool Packed1,
			  class T2, bool Packed2,
			  int Match>
	auto Multiply(const Matrix<T1, Rows1, Match, eMatrixOrder::PRECEDE_VECTOR, eMatrixLayout::ROW_MAJOR, Packed1>& lhs,
				  const Vector<T2, Match, Packed2>& rhs) {
		return ::mathter::LoopUnroll<Rows1>([lhs, rhs](auto... elementIdx) {
			return Vector(Sum(lhs.Row(elementIdx) * rhs)...);
		});
	}

} // namespace impl


//------------------------------------------------------------------------------
// Matrix product
//------------------------------------------------------------------------------

template <class T1, eMatrixLayout Layout1, bool Packed1,
		  class T2, eMatrixLayout Layout2, bool Packed2,
		  int Rows1, int Match, int Columns2, eMatrixOrder Order>
auto operator*(const Matrix<T1, Rows1, Match, Order, Layout1, Packed1>& lhs,
			   const Matrix<T2, Match, Columns2, Order, Layout2, Packed2>& rhs) {
	return impl::Multiply(lhs, rhs);
}

template <class T1, eMatrixLayout Layout1, bool Packed1,
		  class T2, eMatrixLayout Layout2, bool Packed2,
		  int Match, eMatrixOrder Order>
auto& operator*=(Matrix<T1, Match, Match, Order, Layout1, Packed1>& lhs,
				 const Matrix<T2, Match, Match, Order, Layout2, Packed2>& rhs) {
	return lhs = lhs * rhs;
}


//------------------------------------------------------------------------------
// Matrix x matrix elementwise operations
//------------------------------------------------------------------------------

#define MATHTER_MATRIX_ELEMENTWISE(NAME, FUNC)                                 \
	template <class T1, eMatrixLayout Layout1, bool Packed1,                   \
			  class T2, eMatrixLayout Layout2, bool Packed2,                   \
			  int Rows, int Columns, eMatrixOrder Order>                       \
	auto NAME(const Matrix<T1, Rows, Columns, Order, Layout1, Packed1>& lhs,   \
			  const Matrix<T2, Rows, Columns, Order, Layout2, Packed2>& rhs) { \
		return impl::Elementwise(lhs, rhs, FUNC{});                            \
	}


#define MATHTER_MATRIX_ELEMENTWISE_ASSIGN(NAME, FUNC)                           \
	template <class T1, eMatrixLayout Layout1, bool Packed1,                    \
			  class T2, eMatrixLayout Layout2, bool Packed2,                    \
			  int Rows, int Columns, eMatrixOrder Order>                        \
	auto& NAME(Matrix<T1, Rows, Columns, Order, Layout1, Packed1>& lhs,         \
			   const Matrix<T2, Rows, Columns, Order, Layout2, Packed2>& rhs) { \
		return lhs = impl::Elementwise(lhs, rhs, FUNC{});                       \
	}


MATHTER_MATRIX_ELEMENTWISE(operator+, std::plus);
MATHTER_MATRIX_ELEMENTWISE(operator-, std::minus);
MATHTER_MATRIX_ELEMENTWISE(operator/, std::divides);
MATHTER_MATRIX_ELEMENTWISE(Hadamard, std::multiplies);

MATHTER_MATRIX_ELEMENTWISE_ASSIGN(operator+=, std::plus);
MATHTER_MATRIX_ELEMENTWISE_ASSIGN(operator-=, std::minus);
MATHTER_MATRIX_ELEMENTWISE_ASSIGN(operator/=, std::divides);


//------------------------------------------------------------------------------
// Matrix x scalar operations
//------------------------------------------------------------------------------

#define MATHTER_MATRIX_SCALAR(OP, FUNC)                                                                    \
	template <class T1, int Rows1, int Columns1, eMatrixOrder Order1, eMatrixLayout Layout1, bool Packed1, \
			  class T2,                                                                                    \
			  class = std::enable_if_t<is_scalar_v<T2>>>                                                   \
	auto operator OP(const Matrix<T1, Rows1, Columns1, Order1, Layout1, Packed1>& lhs,                     \
					 const T2& rhs) {                                                                      \
		return impl::Scalar(lhs, rhs, FUNC{});                                                             \
	}


#define MATHTER_MATRIX_SCALAR_ASSIGN(OP, FUNC)                                                             \
	template <class T1, int Rows1, int Columns1, eMatrixOrder Order1, eMatrixLayout Layout1, bool Packed1, \
			  class T2,                                                                                    \
			  class = std::enable_if_t<is_scalar_v<T2>>>                                                   \
	auto& operator OP(Matrix<T1, Rows1, Columns1, Order1, Layout1, Packed1>& lhs,                          \
					  const T2 & rhs) {                                                                    \
		return lhs = impl::Scalar(lhs, rhs, FUNC{});                                                       \
	}


#define MATHTER_MATRIX_SCALAR_REVERSE(OP, FUNC)                                                            \
	template <class T1,                                                                                    \
			  class T2, int Rows2, int Columns2, eMatrixOrder Order2, eMatrixLayout Layout2, bool Packed2, \
			  class = std::enable_if_t<is_scalar_v<T1>>>                                                   \
	auto operator OP(const T1& lhs,                                                                        \
					 const Matrix<T2, Rows2, Columns2, Order2, Layout2, Packed2>& rhs) {                   \
		return impl::Scalar(lhs, rhs, FUNC{});                                                             \
	}


MATHTER_MATRIX_SCALAR(*, std::multiplies);
MATHTER_MATRIX_SCALAR(/, std::divides);
MATHTER_MATRIX_SCALAR(+, std::plus);
MATHTER_MATRIX_SCALAR(-, std::minus);

MATHTER_MATRIX_SCALAR_ASSIGN(*=, std::multiplies);
MATHTER_MATRIX_SCALAR_ASSIGN(/=, std::divides);
MATHTER_MATRIX_SCALAR_ASSIGN(+=, std::plus);
MATHTER_MATRIX_SCALAR_ASSIGN(-=, std::minus);

MATHTER_MATRIX_SCALAR_REVERSE(*, std::multiplies);
MATHTER_MATRIX_SCALAR_REVERSE(/, std::divides);
MATHTER_MATRIX_SCALAR_REVERSE(+, std::plus);
MATHTER_MATRIX_SCALAR_REVERSE(-, std::minus);


//------------------------------------------------------------------------------
// Matrix x vector operations
//------------------------------------------------------------------------------

template <class T1, bool Packed1,
		  class T2, int Columns2, eMatrixLayout Layout2, bool Packed2,
		  int Match>
auto operator*(const Vector<T1, Match, Packed1>& lhs,
			   const Matrix<T2, Match, Columns2, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed2>& rhs) {
	return impl::Multiply(lhs, rhs);
}


template <class T1, int Rows1, eMatrixLayout Layout1, bool Packed1,
		  class T2, bool Packed2,
		  int Match>
auto operator*(const Matrix<T1, Rows1, Match, eMatrixOrder::PRECEDE_VECTOR, Layout1, Packed1>& lhs,
			   const Vector<T2, Match, Packed2>& rhs) {
	return impl::Multiply(lhs, rhs);
}


/// <summary> Applies a homogeneous coordinate transform to the vector. </summary>
/// <remarks>
/// The vector is augmented by appending a 1 to make it a homogeneous vector.
///	After multiplying the augmented vector, the perspective division is performed,
///	and the truncated non-homogeneous vector is returned.
/// </remarks>
template <class T1, bool Packed1,
		  class T2, eMatrixLayout Layout2, bool Packed2,
		  int Match>
auto operator*(const Vector<T1, Match - 1, Packed1>& lhs,
			   const Matrix<T2, Match, Match, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed2>& rhs) {
	const auto homogeneous = impl::Multiply(Vector(lhs, static_cast<T1>(1)), rhs);
	using VecHomo = std::decay_t<decltype(homogeneous)>;
	using VecOut = Vector<scalar_type_t<VecHomo>, dimension_v<VecHomo> - 1, is_packed_v<VecHomo>>;
	return VecOut(homogeneous);
}


/// <summary> Applies a homogeneous coordinate transform to the vector. </summary>
/// <remarks>
/// The vector is augmented by appending a 1 to make it a homogeneous vector.
///	After multiplying the augmented vector, the perspective division is performed,
///	and the truncated non-homogeneous vector is returned.
/// </remarks>
template <class T1, eMatrixLayout Layout1, bool Packed1,
		  class T2, bool Packed2,
		  int Match>
auto operator*(const Matrix<T1, Match, Match, eMatrixOrder::PRECEDE_VECTOR, Layout1, Packed1>& lhs,
			   const Vector<T2, Match - 1, Packed2>& rhs) {
	const auto homogeneous = impl::Multiply(lhs, Vector(rhs, static_cast<T1>(1)));
	using VecHomo = std::decay_t<decltype(homogeneous)>;
	using VecOut = Vector<scalar_type_t<VecHomo>, dimension_v<VecHomo> - 1, is_packed_v<VecHomo>>;
	return VecOut(homogeneous);
}


/// <summary> Applies an affine transform to the vector. </summary>
/// <remarks> The vector is augmented by appending a 1, and then multiplied by the matrix. </remarks>
template <class T1, bool Packed1,
		  class T2, eMatrixLayout Layout2, bool Packed2,
		  int Match>
auto operator*(const Vector<T1, Match - 1, Packed1>& lhs,
			   const Matrix<T2, Match, Match - 1, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed2>& rhs) {
	return impl::Multiply(Vector(lhs, static_cast<T1>(1)), rhs);
}


/// <summary> Applies an affine transform to the vector. </summary>
/// <remarks> The vector is augmented by appending a 1, and then multiplied by the matrix. </remarks>
template <class T1, eMatrixLayout Layout1, bool Packed1,
		  class T2, bool Packed2,
		  int Match>
auto operator*(const Matrix<T1, Match - 1, Match, eMatrixOrder::PRECEDE_VECTOR, Layout1, Packed1>& lhs,
			   const Vector<T2, Match - 1, Packed2>& rhs) {
	return impl::Multiply(lhs, Vector(rhs, static_cast<T1>(1)));
}


template <class T1, bool Packed1,
		  class T2, eMatrixLayout Layout2, bool Packed2,
		  int Match>
auto& operator*=(Vector<T1, Match, Packed1>& lhs,
				 const Matrix<T2, Match, Match, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed2>& rhs) {
	return lhs = lhs * rhs;
}


template <class T1, bool Packed1,
		  class T2, eMatrixLayout Layout2, bool Packed2,
		  int Match>
auto& operator*=(Vector<T1, Match - 1, Packed1>& lhs,
				 const Matrix<T2, Match, Match, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed2>& rhs) {
	return lhs = lhs * rhs;
}


template <class T1, bool Packed1,
		  class T2, eMatrixLayout Layout2, bool Packed2,
		  int Match>
auto& operator*=(Vector<T1, Match - 1, Packed1>& lhs,
				 const Matrix<T2, Match, Match - 1, eMatrixOrder::FOLLOW_VECTOR, Layout2, Packed2>& rhs) {
	return lhs = lhs * rhs;
}


//------------------------------------------------------------------------------
// Matrix plus & minus
//------------------------------------------------------------------------------

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto operator+(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) {
	return mat;
}

template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
auto operator-(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) {
	return ::mathter::LoopUnroll<std::decay_t<decltype(mat)>::stripeCount>([&mat](auto... stripeIdx) {
		return Matrix<T, Rows, Columns, Order, Layout, Packed>(stripeArg, -mat.stripes[stripeIdx]...);
	});
}

} // namespace mathter