add_executable(UnitTest ${sources} ${headers})

target_sources(UnitTest
    PRIVATE
        "main.cpp"
        #"TestIntersection.cpp"
        "Matrix/TestMatrix.cpp"
        "Matrix/TestArithmetic.cpp"
        "Matrix/TestComparison.cpp"
        #"TestMatrix/TestMatrixVectorArithmetic.cpp"
        #"TestMatrix/TestMatrixDecomposition.cpp"
        #"TestMatrix/TestMatrixFunction.cpp"
        #"TestTransforms.cpp"
        #"TestQuaternion.cpp"
        "Vector/TestArithmetic.cpp"
        "Vector/TestComparison.cpp"
        "Vector/TestConcat.cpp"
        "Vector/TestMath.cpp"
        "Vector/TestSIMDUtil.cpp"
        "Vector/TestSwizzle.cpp"
        "Vector/TestVector.cpp"
        #"TestGenerators.hpp"
        #"TestGeometry.cpp"
)

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest Mathter)
target_link_libraries(UnitTest Catch2::Catch2)