add_executable(UnitTest ${sources} ${headers})

target_sources(UnitTest
    PRIVATE
        "main.cpp"
        "TestIntersection.cpp"
        "TestMatrix/TestMatrixImpl.cpp"
        "TestMatrix/TestMatrixArithmetic.cpp"
        "TestMatrix/TestMatrixVectorArithmetic.cpp"
        "TestMatrix/TestMatrixDecomposition.cpp"
        "TestMatrix/TestMatrixFunction.cpp"
        "TestTransforms.cpp"
        "TestQuaternion.cpp"
        "TestVector/TestVectorImpl.cpp"
        "TestVector/TestVectorFunction.cpp"
        "TestVector/TestVectorConcat.cpp"
        "TestVector/TestVectorArithmetic.cpp"
        "TestGenerators.hpp"
        "TestGeometry.cpp"
)       

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest Mathter)
target_link_libraries(UnitTest Catch2::Catch2)