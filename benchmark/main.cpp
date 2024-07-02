// L=============================================================================
// L This software is distributed under the MIT license.
// L Copyright 2021 Péter Kardos
// L=============================================================================

#pragma warning(disable : 4244)

#include <Mathter/Vector.hpp>

#include <iostream>
#include <xsimd/xsimd.hpp>

#define CATCH_CONFIG_RUNNER
#include <catch2/catch_session.hpp>

using namespace std;
using namespace mathter;

int main(int argc, char* argv[]) {
	std::cout << "Available on CPU: " << std::endl;
	const auto architectures = xsimd::available_architectures();
	std::cout << "  sse2: " << (architectures.sse2 ? "YES" : "NO") << std::endl;
	std::cout << "  sse3: " << (architectures.sse3 ? "YES" : "NO") << std::endl;
	std::cout << "  ssse3: " << (architectures.ssse3 ? "YES" : "NO") << std::endl;
	std::cout << "  sse4_1: " << (architectures.sse4_1 ? "YES" : "NO") << std::endl;
	std::cout << "  sse4_2: " << (architectures.sse4_2 ? "YES" : "NO") << std::endl;
	std::cout << "  sse4a: " << (architectures.sse4a ? "YES" : "NO") << std::endl;
	std::cout << "  fma3_sse: " << (architectures.fma3_sse ? "YES" : "NO") << std::endl;
	std::cout << "  fma4: " << (architectures.fma4 ? "YES" : "NO") << std::endl;
	std::cout << "  xop: " << (architectures.xop ? "YES" : "NO") << std::endl;
	std::cout << "  avx: " << (architectures.avx ? "YES" : "NO") << std::endl;
	std::cout << "  fma3_avx: " << (architectures.fma3_avx ? "YES" : "NO") << std::endl;
	std::cout << "  avx2: " << (architectures.avx2 ? "YES" : "NO") << std::endl;
	std::cout << "  fma3_avx2: " << (architectures.fma3_avx2 ? "YES" : "NO") << std::endl;
	std::cout << "  avx512f: " << (architectures.avx512f ? "YES" : "NO") << std::endl;
	std::cout << "  avx512cd: " << (architectures.avx512cd ? "YES" : "NO") << std::endl;
	std::cout << "  avx512dq: " << (architectures.avx512dq ? "YES" : "NO") << std::endl;
	std::cout << "  avx512bw: " << (architectures.avx512bw ? "YES" : "NO") << std::endl;
	std::cout << "  neon: " << (architectures.neon ? "YES" : "NO") << std::endl;
	std::cout << "  neon64: " << (architectures.neon64 ? "YES" : "NO") << std::endl;
	std::cout << "  sve: " << (architectures.sve ? "YES" : "NO") << std::endl;
	std::cout << std::endl;

	std::cout << "Enabled in build: " << std::endl;
	xsimd::all_architectures::for_each([](const auto& arch) {
		if (arch.supported()) {
			std::cout << "  " << arch.name() << std::endl;
		}
	});
	std::cout << std::endl;

	std::cout << "SIMD support:" << std::endl;
	std::cout << "  float.2: " << (IsBatched<float, 2, false>() ? "YES" : "NO") << " - " << sizeof(Vector<float, 2, false>) << " bytes" << std::endl;
	std::cout << "  float.3: " << (IsBatched<float, 3, false>() ? "YES" : "NO") << " - " << sizeof(Vector<float, 3, false>) << " bytes" << std::endl;
	std::cout << "  float.4: " << (IsBatched<float, 4, false>() ? "YES" : "NO") << " - " << sizeof(Vector<float, 4, false>) << " bytes" << std::endl;
	std::cout << "  double.2: " << (IsBatched<double, 2, false>() ? "YES" : "NO") << " - " << sizeof(Vector<double, 2, false>) << " bytes" << std::endl;
	std::cout << "  double.3: " << (IsBatched<double, 3, false>() ? "YES" : "NO") << " - " << sizeof(Vector<double, 3, false>) << " bytes" << std::endl;
	std::cout << "  double.4: " << (IsBatched<double, 4, false>() ? "YES" : "NO") << " - " << sizeof(Vector<double, 4, false>) << " bytes" << std::endl;
	std::cout << std::endl;

	int ret = Catch::Session().run(argc, argv);
	return ret;
}