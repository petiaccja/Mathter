# Enumerate sources
set(sources
        main.cpp
        TestIntersection.cpp
        TestMatrix.cpp
        TestQuaternion.cpp
        TestVectorImpl.cpp
        TestVectorFunction.cpp
        TestVectorConcat.cpp
        TestVectorArithmetic.cpp
        TestGenerators.hpp
        TestGeometry.cpp)


file(GLOB mathter_headers Mathter/*.hpp)
set (mathter_headers ${mathter_headers} Mathter/Mathter.natvis)

add_executable(UnitTest ${sources} ${headers})

# Source groups
source_group("Test" FILES ${sources} ${headers})
source_group("Mathter" FILES ${mathter_headers})