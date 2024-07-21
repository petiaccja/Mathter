add_executable(UnitTest ${sources} ${headers})

target_sources(UnitTest
    PRIVATE
        "main.cpp"
        #"TestIntersection.cpp"
        #"TestMatrix/TestMatrixImpl.cpp"
        #"TestMatrix/TestMatrixArithmetic.cpp"
        #"TestMatrix/TestMatrixCompare.cpp"
        #"TestMatrix/TestMatrixVectorArithmetic.cpp"
        #"TestMatrix/TestMatrixDecomposition.cpp"
        #"TestMatrix/TestMatrixFunction.cpp"
        #"TestTransforms.cpp"
        #"TestQuaternion.cpp"
        "Vector/TestArithmetic.cpp"
        "Vector/TestComparison.cpp"
        "Vector/TestSIMDUtil.cpp"
        "Vector/TestSwizzle.cpp"
        "Vector/TestVector.cpp"
        "Vector/TestConcat.cpp"
        #"TestVector/TestVectorFunction.cpp"
        #"TestGenerators.hpp"
        #"TestGeometry.cpp"
)

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest Mathter)
target_link_libraries(UnitTest Catch2::Catch2)