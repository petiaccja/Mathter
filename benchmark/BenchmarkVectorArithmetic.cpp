#include "Benchmark.hpp"
#include "Utility.hpp"

#include <Mathter/Vector.hpp>


using namespace mathter;


#define VECTOR_BINOP_BENCHMARK_CASE(TYPE, DIM, PACKED, OP, OPTEXT)                  \
	BENCHMARK_CASE(#TYPE "." #DIM " " OPTEXT " " #TYPE "." #DIM " (P=" #PACKED ")", \
				   "[Vector]",                                                      \
				   50,                                                              \
				   64,                                                              \
				   OP,                                                              \
				   feedBinary,                                                      \
				   (TuplizeArrays(MakeVectorArray<Vector<TYPE, DIM, PACKED>, 32>(), MakeVectorArray<Vector<TYPE, DIM, PACKED>, 32>())));



VECTOR_BINOP_BENCHMARK_CASE(float, 2, false, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, false, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, false, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 2, true, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, true, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, true, opAdd, "+");

VECTOR_BINOP_BENCHMARK_CASE(float, 2, false, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, false, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, false, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 2, true, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, true, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, true, opMul, "*");

VECTOR_BINOP_BENCHMARK_CASE(float, 2, false, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, false, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, false, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 2, true, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 3, true, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(float, 4, true, opDiv, "/");

VECTOR_BINOP_BENCHMARK_CASE(double, 2, false, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, false, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, false, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 2, true, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, true, opAdd, "+");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, true, opAdd, "+");

VECTOR_BINOP_BENCHMARK_CASE(double, 2, false, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, false, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, false, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 2, true, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, true, opMul, "*");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, true, opMul, "*");

VECTOR_BINOP_BENCHMARK_CASE(double, 2, false, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, false, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, false, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 2, true, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 3, true, opDiv, "/");
VECTOR_BINOP_BENCHMARK_CASE(double, 4, true, opDiv, "/");