#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "../Input.hpp"

#include <Mathter/Decompositions/DecomposeLU.hpp>
#include <Mathter/Decompositions/DecomposeQR.hpp>
#include <Mathter/Matrix.hpp>


using namespace mathter;

namespace {


struct NormOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		const auto result = Norm(mat);
		auto copy = mat;
		copy(0, 0) = -result;
		return copy;
	}
};


struct NormPreciseOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		const auto result = NormPrecise(mat);
		auto copy = mat;
		copy(0, 0) = -result;
		return copy;
	}
};


struct TraceOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		const auto result = Trace(mat);
		auto copy = mat;
		copy(0, 0) = -result;
		return copy;
	}
};


struct TransposeOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		return Transpose(mat);
	}
};


struct DeterminantOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		const auto result = Determinant(mat);
		auto copy = mat;
		copy(0, 0) = -result;
		return copy;
	}
};


struct DeterminantBareissOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		const auto result = mathter::impl::BareissAlgorithm(mat);
		auto copy = mat;
		copy(0, 0) = -result;
		return copy;
	}
};


struct InverseOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		return Inverse(mat);
	}
};


struct InverseLUPOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		return DecomposeLUP(mat).Inverse();
	}
};


struct InverseQROp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& mat) const {
		return DecomposeQR(mat).Inverse();
	}
};


static constexpr auto benchmarkCaseLayout_r = eMatrixLayout::ROW_MAJOR;
static constexpr auto benchmarkCaseLayout_c = eMatrixLayout::COLUMN_MAJOR;


#define MATRIX_UNARY_MATH_BENCHMARK_CASE(TYPE, ROWS, COLS, LAYOUT, PACKED, OP, OPTEXT)                                                   \
	BENCHMARK_CASE(#TYPE "." #ROWS #COLS #LAYOUT " " OPTEXT " (P=" #PACKED ")",                                                          \
				   "[Matrix][Math]",                                                                                                     \
				   50,                                                                                                                   \
				   64,                                                                                                                   \
				   GenericNAryFixture{ OP },                                                                                             \
				   MakeRandomInput<Matrix<TYPE, ROWS, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT, PACKED>, 1>()[0], \
				   MakeRandomInput<Matrix<TYPE, ROWS, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT, PACKED>, 4>(),    \
				   std::array<std::monostate, 64>());

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, NormOp{}, "Norm");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, NormOp{}, "Norm");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, NormOp{}, "Norm");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, NormPreciseOp{}, "NormPrecise");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, NormPreciseOp{}, "NormPrecise");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, NormPreciseOp{}, "NormPrecise");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, TraceOp{}, "Trace");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, TraceOp{}, "Trace");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, TraceOp{}, "Trace");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, TransposeOp{}, "Tranpose");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, TransposeOp{}, "Tranpose");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, TransposeOp{}, "Tranpose");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, DeterminantOp{}, "Determinant");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, DeterminantOp{}, "Determinant");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, DeterminantOp{}, "Determinant");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, DeterminantBareissOp{}, "Determinant:Bareiss");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, DeterminantBareissOp{}, "Determinant:Bareiss");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, DeterminantBareissOp{}, "Determinant:Bareiss");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, InverseOp{}, "Inverse");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, InverseOp{}, "Inverse");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, InverseOp{}, "Inverse");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, c, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, c, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, c, false, InverseLUPOp{}, "Invesrse:LUP");

MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, r, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, r, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, r, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 2, 2, c, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 3, 3, c, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(float, 4, 4, c, false, InverseQROp{}, "Inverse:QR");


MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, NormOp{}, "Norm");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, NormOp{}, "Norm");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, NormOp{}, "Norm");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, NormPreciseOp{}, "NormPrecise");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, NormPreciseOp{}, "NormPrecise");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, NormPreciseOp{}, "NormPrecise");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, TraceOp{}, "Trace");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, TraceOp{}, "Trace");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, TraceOp{}, "Trace");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, TransposeOp{}, "Tranpose");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, TransposeOp{}, "Tranpose");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, TransposeOp{}, "Tranpose");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, DeterminantOp{}, "Determinant");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, DeterminantOp{}, "Determinant");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, DeterminantOp{}, "Determinant");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, DeterminantBareissOp{}, "Determinant:Bareiss");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, DeterminantBareissOp{}, "Determinant:Bareiss");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, DeterminantBareissOp{}, "Determinant:Bareiss");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, InverseOp{}, "Inverse");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, InverseOp{}, "Inverse");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, InverseOp{}, "Inverse");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, c, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, c, false, InverseLUPOp{}, "Invesrse:LUP");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, c, false, InverseLUPOp{}, "Invesrse:LUP");

MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, r, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, r, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, r, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 2, 2, c, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 3, 3, c, false, InverseQROp{}, "Inverse:QR");
MATRIX_UNARY_MATH_BENCHMARK_CASE(double, 4, 4, c, false, InverseQROp{}, "Inverse:QR");

} // namespace