#include "Benchmark.hpp"
#include "Utility.hpp"

#include <Mathter/Matrix.hpp>


using namespace mathter;

//------------------------------------------------------------------------------
// Single precision
//------------------------------------------------------------------------------

// Addition
BENCHMARK_CASE("Mat22 ADD Mat22 float vec cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22f_fc, 32>(), MakeMatrixArray<Mat22f_fc, 32>())));

BENCHMARK_CASE("Mat33 ADD Mat33 float vec cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33f_fc, 32>(), MakeMatrixArray<Mat33f_fc, 32>())));

BENCHMARK_CASE("Mat44 ADD Mat44 float vec cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44f_fc, 32>(), MakeMatrixArray<Mat44f_fc, 32>())));

BENCHMARK_CASE("Mat22 ADD Mat22 float packed cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22fp_fc, 32>(), MakeMatrixArray<Mat22fp_fc, 32>())));

BENCHMARK_CASE("Mat33 ADD Mat33 float packed cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33fp_fc, 32>(), MakeMatrixArray<Mat33fp_fc, 32>())));

BENCHMARK_CASE("Mat44 ADD Mat44 float packed cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44fp_fc, 32>(), MakeMatrixArray<Mat44fp_fc, 32>())));

// Multiplication
BENCHMARK_CASE("Mat22 MUL Mat22 float vec cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22f_fc, 32>(), MakeMatrixArray<Mat22f_fc, 32>())));

BENCHMARK_CASE("Mat33 MUL Mat33 float vec cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33f_fc, 32>(), MakeMatrixArray<Mat33f_fc, 32>())));

BENCHMARK_CASE("Mat44 MUL Mat44 float vec cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44f_fc, 32>(), MakeMatrixArray<Mat44f_fc, 32>())));

BENCHMARK_CASE("Mat22 MUL Mat22 float vec rm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22f_fr, 32>(), MakeMatrixArray<Mat22f_fr, 32>())));

BENCHMARK_CASE("Mat33 MUL Mat33 float vec rm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33f_fr, 32>(), MakeMatrixArray<Mat33f_fr, 32>())));

BENCHMARK_CASE("Mat44 MUL Mat44 float vec rm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44f_fr, 32>(), MakeMatrixArray<Mat44f_fr, 32>())));

BENCHMARK_CASE("Mat22 MUL Mat22 float packed cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22fp_fc, 32>(), MakeMatrixArray<Mat22fp_fc, 32>())));

BENCHMARK_CASE("Mat33 MUL Mat33 float packed cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33fp_fc, 32>(), MakeMatrixArray<Mat33fp_fc, 32>())));

BENCHMARK_CASE("Mat44 MUL Mat44 float packed cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44fp_fc, 32>(), MakeMatrixArray<Mat44fp_fc, 32>())));


//------------------------------------------------------------------------------
// Single precision
//------------------------------------------------------------------------------

// Addition
BENCHMARK_CASE("Mat22 ADD Mat22 double vec cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22d_fc, 32>(), MakeMatrixArray<Mat22d_fc, 32>())));

BENCHMARK_CASE("Mat33 ADD Mat33 double vec cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33d_fc, 32>(), MakeMatrixArray<Mat33d_fc, 32>())));

BENCHMARK_CASE("Mat44 ADD Mat44 double vec cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44d_fc, 32>(), MakeMatrixArray<Mat44d_fc, 32>())));

BENCHMARK_CASE("Mat22 ADD Mat22 double packed cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22dp_fc, 32>(), MakeMatrixArray<Mat22dp_fc, 32>())));

BENCHMARK_CASE("Mat33 ADD Mat33 double packed cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33dp_fc, 32>(), MakeMatrixArray<Mat33dp_fc, 32>())));

BENCHMARK_CASE("Mat44 ADD Mat44 double packed cm", "[Matrix]", 50, 64, opAdd, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44dp_fc, 32>(), MakeMatrixArray<Mat44dp_fc, 32>())));

// Multiplication
BENCHMARK_CASE("Mat22 MUL Mat22 double vec cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22d_fc, 32>(), MakeMatrixArray<Mat22d_fc, 32>())));

BENCHMARK_CASE("Mat33 MUL Mat33 double vec cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33d_fc, 32>(), MakeMatrixArray<Mat33d_fc, 32>())));

BENCHMARK_CASE("Mat44 MUL Mat44 double vec cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44d_fc, 32>(), MakeMatrixArray<Mat44d_fc, 32>())));

BENCHMARK_CASE("Mat22 MUL Mat22 double vec rm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22d_fr, 32>(), MakeMatrixArray<Mat22d_fr, 32>())));

BENCHMARK_CASE("Mat33 MUL Mat33 double vec rm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33d_fr, 32>(), MakeMatrixArray<Mat33d_fr, 32>())));

BENCHMARK_CASE("Mat44 MUL Mat44 double vec rm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44d_fr, 32>(), MakeMatrixArray<Mat44d_fr, 32>())));

BENCHMARK_CASE("Mat22 MUL Mat22 double packed cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat22dp_fc, 32>(), MakeMatrixArray<Mat22dp_fc, 32>())));

BENCHMARK_CASE("Mat33 MUL Mat33 double packed cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat33dp_fc, 32>(), MakeMatrixArray<Mat33dp_fc, 32>())));

BENCHMARK_CASE("Mat44 MUL Mat44 double packed cm", "[Matrix]", 50, 64, opMul, feedBinary,
			   (TuplizeArrays(MakeMatrixArray<Mat44dp_fc, 32>(), MakeMatrixArray<Mat44dp_fc, 32>())));