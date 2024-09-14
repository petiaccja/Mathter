// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include <Mathter/Vector.hpp>

#include <iostream>

#if MATHTER_ENABLE_SIMD
#include <xsimd/xsimd.hpp>
#endif

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>

using namespace std;
using namespace mathter;


void PrintArch(std::string_view name, unsigned supported) {
	std::cout << name << (supported ? "YES" : "NO") << std::endl;
}


template <class T, int Dim>
void PrintVectorType(std::string_view name) {
	std::cout << "  " << name << ": " << (IsBatched<T, Dim, false>() ? "YES" : "NO") << " - " << sizeof(Vector<T, Dim, false>) << " bytes" << std::endl;
}


void DisplayArchitectureInfo() {
#if MATHTER_ENABLE_SIMD
	std::cout << "Available on CPU: " << std::endl;
	const auto architectures = xsimd::available_architectures();

	PrintArch("sse2", architectures.sse2);
	PrintArch("sse3", architectures.sse3);
	PrintArch("ssse3", architectures.ssse3);
	PrintArch("sse4_1", architectures.sse4_1);
	PrintArch("sse4_2", architectures.sse4_2);
	PrintArch("fma3_sse42", architectures.fma3_sse42);
	PrintArch("fma4", architectures.fma4);
	PrintArch("avx", architectures.avx);
	PrintArch("fma3_avx", architectures.fma3_avx);
	PrintArch("avx2", architectures.avx2);
	PrintArch("avxvnni", architectures.avxvnni);
	PrintArch("fma3_avx2", architectures.fma3_avx2);
	PrintArch("avx512f", architectures.avx512f);
	PrintArch("avx512cd", architectures.avx512cd);
	PrintArch("avx512dq", architectures.avx512dq);
	PrintArch("avx512bw", architectures.avx512bw);
	PrintArch("avx512er", architectures.avx512er);
	PrintArch("avx512pf", architectures.avx512pf);
	PrintArch("avx512ifma", architectures.avx512ifma);
	PrintArch("avx512vbmi", architectures.avx512vbmi);
	PrintArch("avx512vnni_bw", architectures.avx512vnni_bw);
	PrintArch("avx512vnni_vbmi", architectures.avx512vnni_vbmi);
	PrintArch("neon", architectures.neon);
	PrintArch("neon64", architectures.neon64);
	PrintArch("i8mm_neon64", architectures.i8mm_neon64);
	PrintArch("sve", architectures.sve);
	PrintArch("rvv", architectures.rvv);
	PrintArch("wasm", architectures.wasm);

	std::cout << "Enabled in build: " << std::endl;
	xsimd::all_architectures::for_each([](const auto& arch) {
		if (arch.supported()) {
			std::cout << "  " << arch.name() << std::endl;
		}
	});
	std::cout << std::endl;
#else
	std::cout << "SIMD disabled." << std::endl
			  << std::endl;
#endif
}

int main(int argc, char* argv[]) {
	DisplayArchitectureInfo();

	std::cout << "SIMD support:" << std::endl;
	PrintVectorType<float, 2>("float2");
	PrintVectorType<float, 3>("float3");
	PrintVectorType<float, 4>("float4");
	PrintVectorType<double, 2>("double2");
	PrintVectorType<double, 3>("double3");
	PrintVectorType<double, 4>("double4");
	PrintVectorType<std::complex<float>, 2>("c_float2");
	PrintVectorType<std::complex<float>, 3>("c_float3");
	PrintVectorType<std::complex<float>, 4>("c_float4");
	PrintVectorType<std::complex<double>, 2>("c_double2");
	PrintVectorType<std::complex<double>, 3>("c_double3");
	PrintVectorType<std::complex<double>, 4>("c_double4");
	PrintVectorType<int32_t, 2>("i32_2");
	PrintVectorType<int32_t, 3>("i32_3");
	PrintVectorType<int32_t, 4>("i32_4");
	PrintVectorType<int64_t, 2>("i64_2");
	PrintVectorType<int64_t, 3>("i64_3");
	PrintVectorType<int64_t, 4>("i64_4");
	std::cout << std::endl;

	int ret = Catch::Session().run(argc, argv);
	return ret;
}