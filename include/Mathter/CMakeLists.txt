add_library(Mathter INTERFACE)

target_sources(Mathter
	INTERFACE FILE_SET natvis TYPE HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/.." FILES
		"Mathter.natvis"
	INTERFACE FILE_SET headers TYPE HEADERS BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/.." FILES
		# Master headers
		"Geometry.hpp"
		"IoStream.hpp"
		"Matrix.hpp"
		"Quaternion.hpp"
		"Transforms.hpp"
		"Utility.hpp"
		"Vector.hpp"
		# Common
		"Common/DeterministicInitializer.hpp"
		"Common/Functional.hpp"
		"Common/LoopUtil.hpp"
		"Common/MathUtil.hpp"
		"Common/Range.hpp"
		"Common/Types.hpp"
		"Common/TypeTraits.hpp"
		# Decompositions
		"Decompositions/DecomposeLU.hpp"
		"Decompositions/DecomposeQR.hpp"
		# SVD will be updated separately: "Decompositions/DecomposeSVD.hpp"
		# Geometry
		"Geometry/BezierCurve.hpp"
		"Geometry/Hyperplane.hpp"
		"Geometry/Intersection.hpp"
		"Geometry/Line.hpp"
		"Geometry/LineSegment.hpp"
		"Geometry/Ray.hpp"
		"Geometry/Triangle.hpp"
		# Matrix
		"Matrix/Algorithm.hpp"
		"Matrix/Arithmetic.hpp"
		"Matrix/Cast.hpp"
		"Matrix/Comparison.hpp"
		"Matrix/Math.hpp"
		"Matrix/Matrix.hpp"
		# Quaternion
		"Quaternion/Arithmetic.hpp"
		"Quaternion/Comparison.hpp"
		"Quaternion/Literals.hpp"
		"Quaternion/Math.hpp"
		"Quaternion/Quaternion.hpp"
		"Quaternion/RotationArithmetic.hpp"
		# Transforms
		"Transforms/IdentityBuilder.hpp"
		"Transforms/OrthographicBuilder.hpp"
		"Transforms/PerspectiveBuilder.hpp"
		"Transforms/RandomBuilder.hpp"
		"Transforms/Rotation2DBuilder.hpp"
		"Transforms/Rotation3DBuilder.hpp"
		"Transforms/ScaleBuilder.hpp"
		"Transforms/ShearBuilder.hpp"
		"Transforms/TranslationBuilder.hpp"
		"Transforms/ViewBuilder.hpp"
		"Transforms/ZeroBuilder.hpp"
		# Vector
		"Vector/Arithmetic.hpp"
		"Vector/Comparison.hpp"
		"Vector/Concat.hpp"
		"Vector/Math.hpp"
		"Vector/OperationUtil.hpp"
		"Vector/SIMDUtil.hpp"
		"Vector/Swizzle.hpp"
		"Vector/Vector.hpp"
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