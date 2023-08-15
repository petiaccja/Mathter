#pragma once

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <tuple>
#include <vector>

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_X86))
#include <intrin.h>
#elif defined(__GLIBC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>
#else
#include <chrono>
#define MATHTER_TSC_USES_CHRONO
#endif

#ifdef _MSC_VER
#define MATHTER_NOINLINE __declspec(noinline)
#else
#define MATHTER_NOINLINE __attribute__((noinline))
#endif


namespace impl {


struct BenchmarkRecord {
	std::string name;
	float latency;
	float throughput;
};


extern std::vector<BenchmarkRecord> g_records;
extern std::mutex g_mutex;


inline int64_t ReadTSC() {
#ifdef MATHTER_TSC_USES_CHRONO
	using namespace std::chrono;
	const auto now = high_resolution_clock::now();
	return duration_cast<nanoseconds>(now.time_since_epoch()).count();
#else
	unsigned aux;
	return __rdtscp(&aux);
#endif
}


template <class T>
void DoNotOptimizeAway(T&& value) {
	// Idea from celero https://github.com/DigitalInBlue/Celero/blob/master/include/celero/Utilities.h
	if (ReadTSC() == 0) {
		alignas(T) thread_local std::array<volatile char, sizeof(T)> c;
		c = *reinterpret_cast<const std::array<volatile char, sizeof(T)>*>(&value);
	}
}


template <size_t... Indices, class Func>
void Unroll(Func func, std::index_sequence<Indices...>) {
	(..., func(Indices));
}


template <size_t Count, class Func>
void Unroll(Func func) {
	Unroll(func, std::make_index_sequence<Count>{});
}


template <class DoSample>
MATHTER_NOINLINE int64_t BestSample(DoSample&& doSample, int64_t samples) {
	int64_t bestTime = std::numeric_limits<int64_t>::max();
	for (size_t i = 0; i < samples; ++i) {
		bestTime = std::min(bestTime, doSample());
	}
	return bestTime;
}


template <class Operation, class Feed, class... Init, size_t N>
float Latency(int64_t samples, int64_t repeat, Operation&& operation, Feed&& feed, const std::array<std::tuple<Init...>, N>& init) {
	static constexpr size_t unrollCount = 16;

	auto doSample = [&]() {
		auto result = operation(init[0]);
		const auto seed = init[1 % N];

		const auto startTime = ReadTSC();
		for (int64_t i = 0; i < repeat; ++i) {
			Unroll<unrollCount>([&](size_t) {
				result = operation(feed(result, seed));
			});
		}
		DoNotOptimizeAway(result);

		const auto endTime = ReadTSC();
		return endTime - startTime;
	};

	return BestSample(doSample, samples) / float(repeat * unrollCount);
}


template <class Operation, class Feed, class... Init, size_t N>
float Throughput(int64_t samples, int64_t repeat, Operation&& operation, Feed&& feed, const std::array<std::tuple<Init...>, N>& init) {
	static constexpr size_t unrollCount = 16;

	auto doSample = [&]() {
		size_t index = 0;
		std::array<decltype(operation(init[0])), unrollCount> results;
		std::fill(results.begin(), results.end(), operation(init[0]));
		const auto seed = init[1 % N];

		const auto startTime = ReadTSC();

		for (int64_t i = 0; i < repeat; ++i) {
			Unroll<results.size()>([&](size_t unrollIdx) {
				results[unrollIdx] = operation(feed(results[unrollIdx], seed));
			});
		}
		DoNotOptimizeAway(results);

		const auto endTime = ReadTSC();
		return endTime - startTime;
	};

	return BestSample(doSample, samples) / float(repeat * unrollCount);
}


template <class Operation, class Feed, class... Init, size_t N>
void BenchmarkCase(std::string_view name, int64_t samples, int64_t repeat, Operation&& operation, Feed&& feed, const std::array<std::tuple<Init...>, N>& init) {
	const auto latency = Latency(samples, repeat, operation, feed, init);
	const auto throughput = Throughput(samples, repeat, operation, feed, init);
	std::lock_guard lk{ g_mutex };
	g_records.push_back(BenchmarkRecord{ std::string(name), latency, throughput });
}


} // namespace impl


#define BENCHMARK_CASE(NAME, TAG, SAMPLES, REPEAT, OPERATION, FEED, INIT)    \
	TEST_CASE(NAME, TAG) {                                                   \
		::impl::BenchmarkCase(NAME, SAMPLES, REPEAT, OPERATION, FEED, INIT); \
	}