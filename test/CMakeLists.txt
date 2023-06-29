add_executable(UnitTest ${sources} ${headers})

target_sources(UnitTest
    PRIVATE
        "main.cpp"
        "TestIntersection.cpp"
        "TestMatrixImpl.cpp"
        "TestMatrixArithmetic.cpp"
        "TestMatrixVectorArithmetic.cpp"
        "TestMatrixDecomposition.cpp"
        "TestMatrixFunction.cpp"
        "TestTransforms.cpp"
        "TestQuaternion.cpp"
        "TestVectorImpl.cpp"
        "TestVectorFunction.cpp"
        "TestVectorConcat.cpp"
        "TestVectorArithmetic.cpp"
        "TestGenerators.hpp"
        "TestGeometry.cpp"
)       

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest Mathter)
target_link_libraries(UnitTest Catch2::Catch2)