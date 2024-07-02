add_library(Mathter INTERFACE)

target_sources(Mathter
	INTERFACE FILE_SET natvis TYPE HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/.." FILES
		"Mathter.natvis"
	INTERFACE FILE_SET headers TYPE HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/.." FILES
		"Geometry.hpp"
		"IoStream.hpp"
		"Matrix.hpp"
		"Quaternion.hpp"
		"Utility.hpp"
		"Vector.hpp"
		"Common/Approx.hpp"
		"Common/DeterministicInitializer.hpp"
		"Common/Functional.hpp"
		"Common/MathUtil.hpp"
		"Common/Range.hpp"
		"Common/Types.hpp"
		"Common/TypeTraits.hpp"
		"Decompositions/DecomposeLU.hpp"
		"Decompositions/DecomposeQR.hpp"
		"Decompositions/DecomposeSVD.hpp"
		"Matrix/MatrixArithmetic.hpp"
		"Matrix/MatrixCast.hpp"
		"Matrix/MatrixCompare.hpp"
		"Matrix/MatrixFunction.hpp"
		"Matrix/MatrixImpl.hpp"
		"Matrix/MatrixVectorArithmetic.hpp"
		"Quaternion/QuaternionArithmetic.hpp"
		"Quaternion/QuaternionCompare.hpp"
		"Quaternion/QuaternionFunction.hpp"
		"Quaternion/QuaternionImpl.hpp"
		"Quaternion/QuaternionLiterals.hpp"
		"Quaternion/QuaternionVectorArithmetic.hpp"
		"Swizzle/Swizzle_1.inc.hpp"
		"Swizzle/Swizzle_2.inc.hpp"
		"Swizzle/Swizzle_3.inc.hpp"
		"Swizzle/Swizzle_4.inc.hpp"
		"Transforms/IdentityBuilder.hpp"
		"Transforms/OrthographicBuilder.hpp"
		"Transforms/PerspectiveBuilder.hpp"
		"Transforms/Rotation2DBuilder.hpp"
		"Transforms/Rotation3DBuilder.hpp"
		"Transforms/ScaleBuilder.hpp"
		"Transforms/ShearBuilder.hpp"
		"Transforms/TranslationBuilder.hpp"
		"Transforms/ViewBuilder.hpp"
		"Transforms/ZeroBuilder.hpp"
		"Vector/VectorArithmetic.hpp"
		"Vector/VectorCompare.hpp"
		"Vector/VectorConcat.hpp"
		"Vector/VectorFunction.hpp"
		"Vector/VectorImpl.hpp"
)

set_target_properties(Mathter PROPERTIES VERIFY_INTERFACE_HEADER_SETS ON)

target_compile_features(Mathter INTERFACE cxx_std_17)

if (MATHTER_INIT_MODE STREQUAL "NULL")
	target_compile_definitions(Mathter INTERFACE MATHTER_NULL_INITIALIZE=1)
elseif(MATHTER_INIT_MODE STREQUAL "INVALID")
	target_compile_definitions(Mathter INTERFACE MATHTER_INVALID_INITIALIZE=1)
elseif(MATHTER_INIT_MODE STREQUAL "UNINITIALIZED")
	target_compile_definitions(Mathter INTERFACE MATHTER_DONT_INITIALIZE=1)
endif()

if (${MATHTER_ENABLE_SIMD})
	message("Vectorization: ON (using XSimd).")
	find_package(xsimd REQUIRED)
	target_link_libraries(Mathter INTERFACE xsimd)
	target_compile_definitions(Mathter INTERFACE MATHTER_ENABLE_SIMD=1)
else()
	message("Vectorization: OFF.")
endif()