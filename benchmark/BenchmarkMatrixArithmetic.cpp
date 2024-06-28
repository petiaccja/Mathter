#include "Benchmark.hpp"
#include "Utility.hpp"

#include <Mathter/Matrix.hpp>


using namespace mathter;


static constexpr auto benchmarkCaseLayout_r = eMatrixLayout::ROW_MAJOR;
static constexpr auto benchmarkCaseLayout_c = eMatrixLayout::COLUMN_MAJOR;


#define MATRIX_MUL_BENCHMARK_CASE(TYPE, ROWS, MATCH, COLS, LAYOUT_L, LAYOUT_R, PACKED)                                                                   \
	BENCHMARK_CASE(#TYPE "." #ROWS #MATCH #LAYOUT_L " * " #TYPE "." #MATCH #COLS #LAYOUT_R " (P=" #PACKED ")",                                           \
				   "[Matrix]",                                                                                                                           \
				   50,                                                                                                                                   \
				   64,                                                                                                                                   \
				   opMul,                                                                                                                                \
				   feedBinary,                                                                                                                           \
				   (TuplizeArrays(MakeMatrixArray<Matrix<TYPE, ROWS, MATCH, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_L, PACKED>, 32>(), \
								  MakeMatrixArray<Matrix<TYPE, MATCH, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_R, PACKED>, 32>())));


#define MATRIX_EWISE_BENCHMARK_CASE(TYPE, ROWS, COLS, LAYOUT_L, LAYOUT_R, PACKED)                                                                       \
	BENCHMARK_CASE(#TYPE "." #ROWS #COLS #LAYOUT_L " + " #TYPE "." #ROWS #COLS #LAYOUT_R " (P=" #PACKED ")",                                            \
				   "[Matrix]",                                                                                                                          \
				   50,                                                                                                                                  \
				   64,                                                                                                                                  \
				   opAdd,                                                                                                                               \
				   feedBinary,                                                                                                                          \
				   (TuplizeArrays(MakeMatrixArray<Matrix<TYPE, ROWS, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_L, PACKED>, 32>(), \
								  MakeMatrixArray<Matrix<TYPE, ROWS, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_R, PACKED>, 32>())));


MATRIX_MUL_BENCHMARK_CASE(float, 2, 2, 2, r, r, false);
MATRIX_MUL_BENCHMARK_CASE(float, 3, 3, 3, r, r, false);
MATRIX_MUL_BENCHMARK_CASE(float, 4, 4, 4, r, r, false);
MATRIX_MUL_BENCHMARK_CASE(float, 2, 2, 2, r, r, true);
MATRIX_MUL_BENCHMARK_CASE(float, 3, 3, 3, r, r, true);
MATRIX_MUL_BENCHMARK_CASE(float, 4, 4, 4, r, r, true);

MATRIX_MUL_BENCHMARK_CASE(float, 2, 2, 2, c, c, false);
MATRIX_MUL_BENCHMARK_CASE(float, 3, 3, 3, c, c, false);
MATRIX_MUL_BENCHMARK_CASE(float, 4, 4, 4, c, c, false);
MATRIX_MUL_BENCHMARK_CASE(float, 2, 2, 2, c, c, true);
MATRIX_MUL_BENCHMARK_CASE(float, 3, 3, 3, c, c, true);
MATRIX_MUL_BENCHMARK_CASE(float, 4, 4, 4, c, c, true);

MATRIX_MUL_BENCHMARK_CASE(double, 2, 2, 2, r, r, false);
MATRIX_MUL_BENCHMARK_CASE(double, 3, 3, 3, r, r, false);
MATRIX_MUL_BENCHMARK_CASE(double, 4, 4, 4, r, r, false);
MATRIX_MUL_BENCHMARK_CASE(double, 2, 2, 2, r, r, true);
MATRIX_MUL_BENCHMARK_CASE(double, 3, 3, 3, r, r, true);
MATRIX_MUL_BENCHMARK_CASE(double, 4, 4, 4, r, r, true);

MATRIX_MUL_BENCHMARK_CASE(double, 2, 2, 2, c, c, false);
MATRIX_MUL_BENCHMARK_CASE(double, 3, 3, 3, c, c, false);
MATRIX_MUL_BENCHMARK_CASE(double, 4, 4, 4, c, c, false);
MATRIX_MUL_BENCHMARK_CASE(double, 2, 2, 2, c, c, true);
MATRIX_MUL_BENCHMARK_CASE(double, 3, 3, 3, c, c, true);
MATRIX_MUL_BENCHMARK_CASE(double, 4, 4, 4, c, c, true);

MATRIX_EWISE_BENCHMARK_CASE(float, 2, 2, c, c, false);
MATRIX_EWISE_BENCHMARK_CASE(float, 3, 3, c, c, false);
MATRIX_EWISE_BENCHMARK_CASE(float, 4, 4, c, c, false);
MATRIX_EWISE_BENCHMARK_CASE(float, 2, 2, c, c, true);
MATRIX_EWISE_BENCHMARK_CASE(float, 3, 3, c, c, true);
MATRIX_EWISE_BENCHMARK_CASE(float, 4, 4, c, c, true);

MATRIX_EWISE_BENCHMARK_CASE(double, 2, 2, c, c, false);
MATRIX_EWISE_BENCHMARK_CASE(double, 3, 3, c, c, false);
MATRIX_EWISE_BENCHMARK_CASE(double, 4, 4, c, c, false);
MATRIX_EWISE_BENCHMARK_CASE(double, 2, 2, c, c, true);
MATRIX_EWISE_BENCHMARK_CASE(double, 3, 3, c, c, true);
MATRIX_EWISE_BENCHMARK_CASE(double, 4, 4, c, c, true);