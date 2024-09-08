#include "../Benchmark.hpp"
#include "../Fixtures.hpp"
#include "../Input.hpp"

#include <Mathter/Decompositions/DecomposeLU.hpp>
#include <Mathter/Decompositions/DecomposeQR.hpp>
#include <Mathter/Decompositions/DecomposeSVD.hpp>
#include <Mathter/Matrix.hpp>


using namespace mathter;

namespace {


struct LUOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& salt, const Matrix<T, Rows, Columns, Order, Layout, Packed>& input) const {
		const auto salted = remove_complex_t<T>(0.00001) * salt + input;
		const auto [L, U] = DecomposeLU(salted);
		return L + U;
	}
};


struct LUPOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& salt, const Matrix<T, Rows, Columns, Order, Layout, Packed>& input) const {
		const auto salted = remove_complex_t<T>(0.00001) * salt + input;
		const auto [L, U, P] = DecomposeLUP(salted);
		return (L + U) * remove_complex_t<T>(P[0]);
	}
};


struct QROp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& salt, const Matrix<T, Rows, Columns, Order, Layout, Packed>& input) const {
		const auto salted = remove_complex_t<T>(0.00001) * salt + input;
		const auto [Q, R] = DecomposeQR(salted);
		return Q + R;
	}
};


struct LQOp {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& salt, const Matrix<T, Rows, Columns, Order, Layout, Packed>& input) const {
		const auto salted = remove_complex_t<T>(0.00001) * salt + input;
		const auto [L, Q] = DecomposeLQ(salted);
		return L + Q;
	}
};


struct SVD1Op {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& salt, const Matrix<T, Rows, Columns, Order, Layout, Packed>& input) const {
		const auto salted = remove_complex_t<T>(0.00001) * salt + input;
		const auto [U, S, V] = DecomposeSVD(salted, SVDAlgorithmOneSided);
		return (U + V) * S[0];
	}
};


struct SVD2Op {
	template <class T, int Rows, int Columns, eMatrixOrder Order, eMatrixLayout Layout, bool Packed>
	auto operator()(const Matrix<T, Rows, Columns, Order, Layout, Packed>& salt, const Matrix<T, Rows, Columns, Order, Layout, Packed>& input) const {
		const auto salted = remove_complex_t<T>(0.00001) * salt + input;
		const auto [U, S, V] = DecomposeSVD(salted, SVDAlgorithmTwoSided);
		return (U + V) * S[0];
	}
};


static constexpr auto benchmarkCaseLayout_r = eMatrixLayout::ROW_MAJOR;
static constexpr auto benchmarkCaseLayout_c = eMatrixLayout::COLUMN_MAJOR;


#define DECOMP_BENCHMARK_CASE(TYPE, ROWS, COLS, LAYOUT, PACKED, OP, OPTEXT)                                                              \
	BENCHMARK_CASE(#TYPE "." #ROWS #COLS #LAYOUT " " OPTEXT " (P=" #PACKED ")",                                                          \
				   "[Matrix][Math]",                                                                                                     \
				   50,                                                                                                                   \
				   12,                                                                                                                   \
				   GenericNAryFixture{ OP },                                                                                             \
				   MakeRandomInput<Matrix<TYPE, ROWS, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT, PACKED>, 1>()[0], \
				   MakeRandomInput<Matrix<TYPE, ROWS, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT, PACKED>, 4>(),    \
				   MakeRandomInput<Matrix<TYPE, ROWS, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT, PACKED>, 64>());

DECOMP_BENCHMARK_CASE(float, 2, 2, r, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, r, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, r, false, LUOp{}, "LU decomp");

DECOMP_BENCHMARK_CASE(float, 2, 2, c, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, c, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, c, false, LUOp{}, "LU decomp");

DECOMP_BENCHMARK_CASE(float, 2, 2, r, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, r, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, r, false, LUPOp{}, "LUP decomp");

DECOMP_BENCHMARK_CASE(float, 2, 2, c, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, c, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, c, false, LUPOp{}, "LUP decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, r, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, r, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, r, false, LUOp{}, "LU decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, c, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, c, false, LUOp{}, "LU decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, c, false, LUOp{}, "LU decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, r, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, r, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, r, false, LUPOp{}, "LUP decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, c, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, c, false, LUPOp{}, "LUP decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, c, false, LUPOp{}, "LUP decomp");


DECOMP_BENCHMARK_CASE(float, 2, 2, r, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, r, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, r, false, QROp{}, "QR decomp");

DECOMP_BENCHMARK_CASE(float, 2, 2, c, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, c, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, c, false, QROp{}, "QR decomp");

DECOMP_BENCHMARK_CASE(float, 2, 2, r, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, r, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, r, false, LQOp{}, "LQ decomp");

DECOMP_BENCHMARK_CASE(float, 2, 2, c, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(float, 3, 3, c, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(float, 4, 4, c, false, LQOp{}, "LQ decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, r, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, r, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, r, false, QROp{}, "QR decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, c, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, c, false, QROp{}, "QR decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, c, false, QROp{}, "QR decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, r, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, r, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, r, false, LQOp{}, "LQ decomp");

DECOMP_BENCHMARK_CASE(double, 2, 2, c, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(double, 3, 3, c, false, LQOp{}, "LQ decomp");
DECOMP_BENCHMARK_CASE(double, 4, 4, c, false, LQOp{}, "LQ decomp");


DECOMP_BENCHMARK_CASE(float, 2, 2, r, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(float, 3, 3, r, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(float, 4, 4, r, false, SVD1Op{}, "SVD (1-sided)");

DECOMP_BENCHMARK_CASE(float, 2, 2, c, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(float, 3, 3, c, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(float, 4, 4, c, false, SVD1Op{}, "SVD (1-sided)");

DECOMP_BENCHMARK_CASE(float, 2, 2, r, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(float, 3, 3, r, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(float, 4, 4, r, false, SVD2Op{}, "SVD (2-sided)");

DECOMP_BENCHMARK_CASE(float, 2, 2, c, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(float, 3, 3, c, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(float, 4, 4, c, false, SVD2Op{}, "SVD (2-sided)");

DECOMP_BENCHMARK_CASE(double, 2, 2, r, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(double, 3, 3, r, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(double, 4, 4, r, false, SVD1Op{}, "SVD (1-sided)");

DECOMP_BENCHMARK_CASE(double, 2, 2, c, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(double, 3, 3, c, false, SVD1Op{}, "SVD (1-sided)");
DECOMP_BENCHMARK_CASE(double, 4, 4, c, false, SVD1Op{}, "SVD (1-sided)");

DECOMP_BENCHMARK_CASE(double, 2, 2, r, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(double, 3, 3, r, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(double, 4, 4, r, false, SVD2Op{}, "SVD (2-sided)");

DECOMP_BENCHMARK_CASE(double, 2, 2, c, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(double, 3, 3, c, false, SVD2Op{}, "SVD (2-sided)");
DECOMP_BENCHMARK_CASE(double, 4, 4, c, false, SVD2Op{}, "SVD (2-sided)");



} // namespace