add_library(Mathter INTERFACE)

target_include_directories(Mathter
	INTERFACE 
		$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_sources(Mathter
	INTERFACE FILE_SET natvis TYPE HEADERS FILES
		"Mathter.natvis"
	INTERFACE FILE_SET headers TYPE HEADERS FILES
		"Geometry.hpp"
		"IoStream.hpp"
		"Matrix.hpp"
		"Quaternion.hpp"
		"Utility.hpp"
		"Vector.hpp"
		"Common/Approx.hpp"
		"Common/Definitions.hpp"
		"Common/DeterministicInitializer.hpp"
		"Common/MathUtil.hpp"
		"Common/Range.hpp"
		"Common/Traits.hpp"
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

if (${MATHTER_USE_XSIMD})
	find_package(xsimd)
	target_link_libraries(Mathter INTERFACE xsimd)
endif()