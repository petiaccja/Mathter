#include "Benchmark.hpp"
#include "Utility.hpp"

#include <Mathter/Vector.hpp>


using namespace mathter;


//------------------------------------------------------------------------------
// Single precision
//------------------------------------------------------------------------------

// Addition
BENCHMARK_CASE("Vector<float, 2, VECTORIZED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2f, 32>(), MakeVectorArray<Vec2f, 32>())));

BENCHMARK_CASE("Vector<float, 3, VECTORIZED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3f, 32>(), MakeVectorArray<Vec3f, 32>())));

BENCHMARK_CASE("Vector<float, 4, VECTORIZED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4f, 32>(), MakeVectorArray<Vec4f, 32>())));

BENCHMARK_CASE("Vector<float, 2, PACKED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2fp, 32>(), MakeVectorArray<Vec2fp, 32>())));

BENCHMARK_CASE("Vector<float, 3, PACKED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3fp, 32>(), MakeVectorArray<Vec3fp, 32>())));

BENCHMARK_CASE("Vector<float, 4, PACKED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4fp, 32>(), MakeVectorArray<Vec4fp, 32>())));


// Multiplication
BENCHMARK_CASE("Vector<float, 2, VECTORIZED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2f, 32>(), MakeVectorArray<Vec2f, 32>())));

BENCHMARK_CASE("Vector<float, 3, VECTORIZED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3f, 32>(), MakeVectorArray<Vec3f, 32>())));

BENCHMARK_CASE("Vector<float, 4, VECTORIZED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4f, 32>(), MakeVectorArray<Vec4f, 32>())));

BENCHMARK_CASE("Vector<float, 2, PACKED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2fp, 32>(), MakeVectorArray<Vec2fp, 32>())));

BENCHMARK_CASE("Vector<float, 3, PACKED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3fp, 32>(), MakeVectorArray<Vec3fp, 32>())));

BENCHMARK_CASE("Vector<float, 4, PACKED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4fp, 32>(), MakeVectorArray<Vec4fp, 32>())));


// Division
BENCHMARK_CASE("Vector<float, 2, VECTORIZED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2f, 32>(), MakeVectorArray<Vec2f, 32>())));

BENCHMARK_CASE("Vector<float, 3, VECTORIZED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3f, 32>(), MakeVectorArray<Vec3f, 32>())));

BENCHMARK_CASE("Vector<float, 4, VECTORIZED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4f, 32>(), MakeVectorArray<Vec4f, 32>())));

BENCHMARK_CASE("Vector<float, 2, PACKED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2fp, 32>(), MakeVectorArray<Vec2fp, 32>())));

BENCHMARK_CASE("Vector<float, 3, PACKED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3fp, 32>(), MakeVectorArray<Vec3fp, 32>())));

BENCHMARK_CASE("Vector<float, 4, PACKED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4fp, 32>(), MakeVectorArray<Vec4fp, 32>())));

//------------------------------------------------------------------------------
// Double precision
//------------------------------------------------------------------------------

// Addition
BENCHMARK_CASE("Vector<double, 2, VECTORIZED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2d, 32>(), MakeVectorArray<Vec2d, 32>())));

BENCHMARK_CASE("Vector<double, 3, VECTORIZED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3d, 32>(), MakeVectorArray<Vec3d, 32>())));

BENCHMARK_CASE("Vector<double, 4, VECTORIZED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4d, 32>(), MakeVectorArray<Vec4d, 32>())));

BENCHMARK_CASE("Vector<double, 2, PACKED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2dp, 32>(), MakeVectorArray<Vec2dp, 32>())));

BENCHMARK_CASE("Vector<double, 3, PACKED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3dp, 32>(), MakeVectorArray<Vec3dp, 32>())));

BENCHMARK_CASE("Vector<double, 4, PACKED> ADD", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4dp, 32>(), MakeVectorArray<Vec4dp, 32>())));


// Multiplication
BENCHMARK_CASE("Vector<double, 2, VECTORIZED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2d, 32>(), MakeVectorArray<Vec2d, 32>())));

BENCHMARK_CASE("Vector<double, 3, VECTORIZED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3d, 32>(), MakeVectorArray<Vec3d, 32>())));

BENCHMARK_CASE("Vector<double, 4, VECTORIZED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4d, 32>(), MakeVectorArray<Vec4d, 32>())));

BENCHMARK_CASE("Vector<double, 2, PACKED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2dp, 32>(), MakeVectorArray<Vec2dp, 32>())));

BENCHMARK_CASE("Vector<double, 3, PACKED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3dp, 32>(), MakeVectorArray<Vec3dp, 32>())));

BENCHMARK_CASE("Vector<double, 4, PACKED> MUL", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4dp, 32>(), MakeVectorArray<Vec4dp, 32>())));


// Division
BENCHMARK_CASE("Vector<double, 2, VECTORIZED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2d, 32>(), MakeVectorArray<Vec2d, 32>())));

BENCHMARK_CASE("Vector<double, 3, VECTORIZED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3d, 32>(), MakeVectorArray<Vec3d, 32>())));

BENCHMARK_CASE("Vector<double, 4, VECTORIZED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4d, 32>(), MakeVectorArray<Vec4d, 32>())));

BENCHMARK_CASE("Vector<double, 2, PACKED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2dp, 32>(), MakeVectorArray<Vec2dp, 32>())));

BENCHMARK_CASE("Vector<double, 3, PACKED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3dp, 32>(), MakeVectorArray<Vec3dp, 32>())));

BENCHMARK_CASE("Vector<double, 4, PACKED> DIV", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4dp, 32>(), MakeVectorArray<Vec4dp, 32>())));