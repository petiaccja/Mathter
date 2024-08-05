add_executable(UnitTest ${sources} ${headers})

target_sources(UnitTest
    PRIVATE
        "main.cpp"
        "Geometry/TestBezierCurve.cpp"
        "Geometry/TestHyperplane.cpp"
        "Geometry/TestLine.cpp"
        "Geometry/TestLineSegment.cpp"
        "Geometry/TestRay.cpp"
        "Geometry/TestTriangle.cpp"
        "Matrix/TestMatrix.cpp"
        "Matrix/TestArithmetic.cpp"
        "Matrix/TestComparison.cpp"
        "Matrix/TestMath.cpp"
        "Quaternion/TestArithmetic.cpp"
        "Quaternion/TestComparison.cpp"
        "Quaternion/TestLiterals.cpp"
        "Quaternion/TestMath.cpp"
        "Quaternion/TestQuaternion.cpp"
        "Quaternion/TestRotationArithmetic.cpp"
        "Transforms/TestIdentityBuilder.cpp"
        "Transforms/TestOrthographicBuilder.cpp"
        "Transforms/TestPerspectiveBuilder.cpp"
        "Transforms/TestRandomBuilder.cpp"
        "Transforms/TestRotation2DBuilder.cpp"
        "Transforms/TestRotation3DBuilder.cpp"
        "Transforms/TestScaleBuilder.cpp"
        "Transforms/TestShearBuilder.cpp"
        "Transforms/TestTranslationBuilder.cpp"
        "Transforms/TestViewBuilder.cpp"
        "Transforms/TestZeroBuilder.cpp"
        "Vector/TestArithmetic.cpp"
        "Vector/TestComparison.cpp"
        "Vector/TestConcat.cpp"
        "Vector/TestMath.cpp"
        "Vector/TestSIMDUtil.cpp"
        "Vector/TestSwizzle.cpp"
        "Vector/TestVector.cpp"
        "TestUtils/TestRotation.cpp"
        "TestUtility.cpp"
)

find_package(Catch2 REQUIRED)

target_link_libraries(UnitTest Mathter)
target_link_libraries(UnitTest Catch2::Catch2)