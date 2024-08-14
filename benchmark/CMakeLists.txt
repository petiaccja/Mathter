add_executable(Benchmark ${sources} ${headers})

target_sources(Benchmark
    PRIVATE
        "main.cpp"
        "Benchmark.hpp"
        "Benchmark.cpp"
) 

target_sources(Benchmark
    PRIVATE
        "Matrix/BenchmarkArithmetic.cpp"
        "Scalar/BenchmarkArithmetic.cpp"
        "Scalar/BenchmarkMath.cpp"
        "Vector/BenchmarkArithmetic.cpp"
)   

target_sources(Benchmark
    PRIVATE
        "z_MSVC_Repro/Matrix2x2Multiply.cpp"
)       

find_package(Catch2 REQUIRED)

target_link_libraries(Benchmark Mathter)
target_link_libraries(Benchmark Catch2::Catch2)