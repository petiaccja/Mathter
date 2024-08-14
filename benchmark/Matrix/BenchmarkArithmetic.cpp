#include "../Benchmark.hpp"
#include "../Fixtures.hpp"

#include <Mathter/Matrix.hpp>
#include <Mathter/Transforms/RandomBuilder.hpp>

#include <random>


using namespace mathter;

namespace {

template <class Mat, int Count>
std::array<Mat, Count> MakeInput() {
	using Scalar = scalar_type_t<Mat>;
	using Real = remove_complex_t<Scalar>;
	using namespace std::complex_literals;

	std::mt19937_64 rne;
	std::uniform_real_distribution<Real> rng(Real(-0.5), Real(0.5));

	std::array<Mat, Count> r;
	for (auto& v : r) {
		if constexpr (is_complex_v<Scalar>) {
			const Mat real = Random(rng, rne);
			const Mat imag = Random(rng, rne);
			v = real + Scalar(1if) * imag;
		}
		else {
			v = Random(rng, rne);
		}
	}
	return r;
};


static constexpr auto benchmarkCaseLayout_r = eMatrixLayout::ROW_MAJOR;
static constexpr auto benchmarkCaseLayout_c = eMatrixLayout::COLUMN_MAJOR;


#define MATRIX_BINOP_BENCHMARK_CASE(TYPE, ROWS, MATCH, COLS, LAYOUT_L, LAYOUT_R, PACKED, OP, OPTEXT)                                  \
	BENCHMARK_CASE(#TYPE "." #ROWS #MATCH #LAYOUT_L " " OPTEXT " " #TYPE "." #MATCH #COLS #LAYOUT_R " (P=" #PACKED ")",               \
				   "[Matrix][Arithmetic]",                                                                                            \
				   50,                                                                                                                \
				   64,                                                                                                                \
				   GenericBinaryFixture{ OP{} },                                                                                      \
				   MakeInput<Matrix<TYPE, ROWS, MATCH, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_L, PACKED>, 1>()[0], \
				   MakeInput<Matrix<TYPE, ROWS, MATCH, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_L, PACKED>, 4>(),    \
				   MakeInput<Matrix<TYPE, MATCH, COLS, eMatrixOrder::FOLLOW_VECTOR, benchmarkCaseLayout_##LAYOUT_R, PACKED>, 64>());


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


} // namespace
