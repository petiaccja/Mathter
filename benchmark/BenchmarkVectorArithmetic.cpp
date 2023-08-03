#include "Benchmark.hpp"
#include "Utility.hpp"

#include <Mathter/Vector.hpp>


using namespace mathter;


//------------------------------------------------------------------------------
// Single precision
//------------------------------------------------------------------------------

// Addition
BENCHMARK_CASE("Vec2 ADD float vec", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2f, 32>(), MakeVectorArray<Vec2f, 32>())));

BENCHMARK_CASE("Vec3 ADD float vec", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3f, 32>(), MakeVectorArray<Vec3f, 32>())));

BENCHMARK_CASE("Vec4 ADD float vec", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4f, 32>(), MakeVectorArray<Vec4f, 32>())));

BENCHMARK_CASE("Vec2 ADD float packed", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2fp, 32>(), MakeVectorArray<Vec2fp, 32>())));

BENCHMARK_CASE("Vec3 ADD float packed", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3fp, 32>(), MakeVectorArray<Vec3fp, 32>())));

BENCHMARK_CASE("Vec4 ADD float packed", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4fp, 32>(), MakeVectorArray<Vec4fp, 32>())));


// Multiplication
BENCHMARK_CASE("Vec2 MUL float vec", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2f, 32>(), MakeVectorArray<Vec2f, 32>())));

BENCHMARK_CASE("Vec3 MUL float vec", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3f, 32>(), MakeVectorArray<Vec3f, 32>())));

BENCHMARK_CASE("Vec4 MUL float vec", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4f, 32>(), MakeVectorArray<Vec4f, 32>())));

BENCHMARK_CASE("Vec2 MUL float packed", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2fp, 32>(), MakeVectorArray<Vec2fp, 32>())));

BENCHMARK_CASE("Vec3 MUL float packed", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3fp, 32>(), MakeVectorArray<Vec3fp, 32>())));

BENCHMARK_CASE("Vec4 MUL float packed", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4fp, 32>(), MakeVectorArray<Vec4fp, 32>())));


// Division
BENCHMARK_CASE("Vec2 DIV float vec", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2f, 32>(), MakeVectorArray<Vec2f, 32>())));

BENCHMARK_CASE("Vec3 DIV float vec", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3f, 32>(), MakeVectorArray<Vec3f, 32>())));

BENCHMARK_CASE("Vec4 DIV float vec", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4f, 32>(), MakeVectorArray<Vec4f, 32>())));

BENCHMARK_CASE("Vec2 DIV float packed", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2fp, 32>(), MakeVectorArray<Vec2fp, 32>())));

BENCHMARK_CASE("Vec3 DIV float packed", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3fp, 32>(), MakeVectorArray<Vec3fp, 32>())));

BENCHMARK_CASE("Vec4 DIV float packed", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4fp, 32>(), MakeVectorArray<Vec4fp, 32>())));

//------------------------------------------------------------------------------
// Double precision
//------------------------------------------------------------------------------

// Addition
BENCHMARK_CASE("Vec2 ADD double vec", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2d, 32>(), MakeVectorArray<Vec2d, 32>())));

BENCHMARK_CASE("Vec3 ADD double vec", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3d, 32>(), MakeVectorArray<Vec3d, 32>())));

BENCHMARK_CASE("Vec4 ADD double vec", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4d, 32>(), MakeVectorArray<Vec4d, 32>())));

BENCHMARK_CASE("Vec2 ADD double packed", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2dp, 32>(), MakeVectorArray<Vec2dp, 32>())));

BENCHMARK_CASE("Vec3 ADD double packed", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3dp, 32>(), MakeVectorArray<Vec3dp, 32>())));

BENCHMARK_CASE("Vec4 ADD double packed", "[Vector]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4dp, 32>(), MakeVectorArray<Vec4dp, 32>())));


// Multiplication
BENCHMARK_CASE("Vec2 MUL double vec", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2d, 32>(), MakeVectorArray<Vec2d, 32>())));

BENCHMARK_CASE("Vec3 MUL double vec", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3d, 32>(), MakeVectorArray<Vec3d, 32>())));

BENCHMARK_CASE("Vec4 MUL double vec", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4d, 32>(), MakeVectorArray<Vec4d, 32>())));

BENCHMARK_CASE("Vec2 MUL double packed", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2dp, 32>(), MakeVectorArray<Vec2dp, 32>())));

BENCHMARK_CASE("Vec3 MUL double packed", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3dp, 32>(), MakeVectorArray<Vec3dp, 32>())));

BENCHMARK_CASE("Vec4 MUL double packed", "[Vector]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4dp, 32>(), MakeVectorArray<Vec4dp, 32>())));


// Division
BENCHMARK_CASE("Vec2 DIV double vec", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2d, 32>(), MakeVectorArray<Vec2d, 32>())));

BENCHMARK_CASE("Vec3 DIV double vec", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3d, 32>(), MakeVectorArray<Vec3d, 32>())));

BENCHMARK_CASE("Vec4 DIV double vec", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4d, 32>(), MakeVectorArray<Vec4d, 32>())));

BENCHMARK_CASE("Vec2 DIV double packed", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec2dp, 32>(), MakeVectorArray<Vec2dp, 32>())));

BENCHMARK_CASE("Vec3 DIV double packed", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec3dp, 32>(), MakeVectorArray<Vec3dp, 32>())));

BENCHMARK_CASE("Vec4 DIV double packed", "[Vector]", 50, 64, opDiv, feedBinary,
			   (TuplizeArrays(MakeVectorArray<Vec4dp, 32>(), MakeVectorArray<Vec4dp, 32>())));