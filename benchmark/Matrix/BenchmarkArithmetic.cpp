#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "../Input.hpp"

#include <Mathter/Matrix.hpp>
#include <Mathter/Transforms/RandomBuilder.hpp>

#include <random>


using namespace mathter;

namespace {

static constexpr auto benchmarkCaseLayout_r = eMatrixLayout::ROW_MAJOR;
static constexpr auto benchmarkCaseLayout_c = eMatrixLayout::COLUMN_MAJOR;

static constexpr auto benchmarkCaseOrder_f = eMatrixOrder::FOLLOW_VECTOR;
static constexpr auto benchmarkCaseOrder_p = eMatrixOrder::PRECEDE_VECTOR;


struct MulVecOp {
	template <class Vec, class Mat>
	auto operator()(const Vec& v, const Mat& m) const {
		if constexpr (order_v<Mat> == eMatrixOrder::FOLLOW_VECTOR) {
			return v * m;
		}
		else {
			return m * v;
		}
	}
};


struct OuterProductOp {
	template <class MatLhs, class MatRhs>
	auto operator()(const MatLhs& lhs, const MatRhs& rhs) const {
		const auto result = lhs * rhs;
		const auto idx = size_t(std::real(result(0, 0)) < 0) / column_count_v<MatLhs>;
		return MatLhs(result.Column(idx));
	}
};


#define MATRIX_BINOP_BENCHMARK_CASE(TYPE, ROWS, MATCH, COLS, LAYOUT_L, LAYOUT_R, PACKED, OP, OPTEXT)                                        \
	BENCHMARK_CASE(#TYPE "." #ROWS #MATCH #LAYOUT_L " " OPTEXT " " #TYPE "." #MATCH #COLS #LAYOUT_R " (P=" #PACKED ")",                     \
				   "[Matrix][Arithmetic]",                                                                                                  \
				   50,                                                                                                                      \
				   64,                                                                                                                      \
				   GenericNAryFixture{ OP{} },                                                                                              \
				   MakeRandomInput<Matrix<TYPE, ROWS, MATCH, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_L, PACKED>, 1>()[0], \
				   MakeRandomInput<Matrix<TYPE, ROWS, MATCH, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_L, PACKED>, 4>(),    \
				   MakeRandomInput<Matrix<TYPE, MATCH, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_R, PACKED>, 64>());


#define MATRIX_BINOP_VEC_CASE(TYPE, ROWS, COLS, ORDER, LAYOUT, PACKED, DIM, OP, OPTEXT)                   \
	BENCHMARK_CASE(#TYPE "." #ROWS #COLS #ORDER #LAYOUT " " OPTEXT " " #TYPE "." #DIM " (P=" #PACKED ")", \
				   "[Matrix][Arithmetic]",                                                                \
				   50,                                                                                    \
				   64,                                                                                    \
				   GenericNAryFixture{ OP{} },                                                            \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 1>()[0],                                    \
				   MakeRandomInput<Vector<TYPE, DIM, PACKED>, 4>(),                                       \
				   MakeRandomInput<Matrix<TYPE, ROWS, COLS, benchmarkCaseOrder_##ORDER, benchmarkCaseLayout_##LAYOUT, PACKED>, 64>());


MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, r, r, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, r, r, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, r, r, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, r, r, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, r, r, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, r, r, true, std::multiplies<>, "*");

MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, c, c, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, c, c, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, c, c, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, c, c, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, c, c, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, c, c, true, std::multiplies<>, "*");

MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, r, r, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, r, r, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, r, r, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, r, r, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, r, r, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, r, r, true, std::multiplies<>, "*");

MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, c, c, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, c, c, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, c, c, false, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, c, c, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, c, c, true, std::multiplies<>, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, c, c, true, std::multiplies<>, "*");

MATRIX_BINOP_BENCHMARK_CASE(float, 2, 1, 2, r, r, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 1, 3, r, r, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 1, 4, r, r, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 2, 1, 2, r, r, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 1, 3, r, r, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 1, 4, r, r, true, OuterProductOp, "*");

MATRIX_BINOP_BENCHMARK_CASE(float, 2, 1, 2, c, c, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 1, 3, c, c, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 1, 4, c, c, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 2, 1, 2, c, c, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 1, 3, c, c, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 1, 4, c, c, true, OuterProductOp, "*");

MATRIX_BINOP_BENCHMARK_CASE(double, 3, 1, 3, r, r, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 1, 2, r, r, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 1, 4, r, r, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 1, 2, r, r, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 1, 3, r, r, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 1, 4, r, r, true, OuterProductOp, "*");

MATRIX_BINOP_BENCHMARK_CASE(double, 3, 1, 3, c, c, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 1, 2, c, c, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 1, 4, c, c, false, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 1, 2, c, c, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 1, 3, c, c, true, OuterProductOp, "*");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 1, 4, c, c, true, OuterProductOp, "*");

MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, c, c, false, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, c, c, false, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, c, c, false, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, c, c, true, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, c, c, true, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, c, c, true, std::plus<>, "+");

MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, c, c, false, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, c, c, false, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, c, c, false, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, c, c, true, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, c, c, true, std::plus<>, "+");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, c, c, true, std::plus<>, "+");

MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, c, c, false, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, c, c, false, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, c, c, false, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(float, 2, 2, 2, c, c, true, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(float, 3, 3, 3, c, c, true, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(float, 4, 4, 4, c, c, true, std::divides<>, "/");

MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, c, c, false, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, c, c, false, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, c, c, false, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(double, 2, 2, 2, c, c, true, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(double, 3, 3, 3, c, c, true, std::divides<>, "/");
MATRIX_BINOP_BENCHMARK_CASE(double, 4, 4, 4, c, c, true, std::divides<>, "/");


MATRIX_BINOP_VEC_CASE(float, 2, 2, f, r, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 3, 3, f, r, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, f, r, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, f, r, false, 3, MulVecOp, "*");

MATRIX_BINOP_VEC_CASE(float, 2, 2, f, c, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 3, 3, f, c, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, f, c, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, f, c, false, 3, MulVecOp, "*");

MATRIX_BINOP_VEC_CASE(float, 2, 2, p, r, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 3, 3, p, r, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, p, r, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, p, r, false, 3, MulVecOp, "*");

MATRIX_BINOP_VEC_CASE(float, 2, 2, p, c, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 3, 3, p, c, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, p, c, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(float, 4, 4, p, c, false, 3, MulVecOp, "*");

MATRIX_BINOP_VEC_CASE(double, 2, 2, f, r, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 3, 3, f, r, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, f, r, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, f, r, false, 3, MulVecOp, "*");

MATRIX_BINOP_VEC_CASE(double, 2, 2, f, c, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 3, 3, f, c, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, f, c, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, f, c, false, 3, MulVecOp, "*");

MATRIX_BINOP_VEC_CASE(double, 2, 2, p, r, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 3, 3, p, r, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, p, r, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, p, r, false, 3, MulVecOp, "*");

MATRIX_BINOP_VEC_CASE(double, 2, 2, p, c, false, 2, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 3, 3, p, c, false, 3, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, p, c, false, 4, MulVecOp, "*");
MATRIX_BINOP_VEC_CASE(double, 4, 4, p, c, false, 3, MulVecOp, "*");


} // namespace
