#pragma once

#include <Mathter/Common/OptimizationUtil.hpp>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <vector>

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_X86))
#include <intrin.h>
#elif defined(__GLIBC__) && (defined(__x86_64__) || defined(__i386__))
#include <x86intrin.h>
#else
#include <chrono>
#define MATHTER_TSC_USES_CHRONO
#endif


namespace impl {


struct BenchmarkRecord {
	std::string name;
	double latency;
	double throughput;
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


template <class DoSample>
MATHTER_NOINLINE double BestSample(DoSample&& doSample, int64_t samples) {
	double bestTime = std::numeric_limits<double>::max();
	for (size_t i = 0; i < samples; ++i) {
		bestTime = std::min(bestTime, doSample());
	}
	return bestTime;
}


template <class Fixture, class FirstArg, class... Args>
double LatencySample(int64_t repeat, Fixture&& fixture, FirstArg&& arg, Args&&... args) {
	const auto startTime = ReadTSC();

	auto [result, count] = fixture.Latency(std::forward<FirstArg>(arg), args...);
	for (int64_t i = 0; i < repeat; ++i) {
		auto [result_, count_] = fixture.Latency(std::move(result), args...);
		result = std::move(result_);
		count += count_;
	}
	DoNotOptimizeAway(result);

	const auto endTime = ReadTSC();
	return (endTime - startTime) / double(count);
}


template <class Fixture, class FirstArg, class... Args>
double ThroughputSample(int64_t repeat, Fixture&& fixture, FirstArg&& arg, Args&&... args) {
	const auto startTime = ReadTSC();

	auto [lanes, count] = fixture.Throughput(std::forward<FirstArg>(arg), args...);

	for (int64_t i = 0; i < repeat; ++i) {
		auto [lanes_, count_] = fixture.Throughput(lanes, args...);
		count += count_;
		lanes = lanes_;
	}
	DoNotOptimizeAway(lanes);

	const auto endTime = ReadTSC();
	return (endTime - startTime) / double(count);
}


template <class Fixture, class FirstArg, class... Args>
double Latency(int64_t samples, int64_t repeat, Fixture&& fixture, FirstArg&& arg, Args&&... args) {
	return BestSample([&]() { return LatencySample(repeat, fixture, arg, args...); }, samples);
}


template <class Fixture, class FirstArg, class... Args>
double Throughput(int64_t samples, int64_t repeat, Fixture&& fixture, FirstArg&& arg, Args&&... args) {
	return BestSample([&]() { return ThroughputSample(repeat, fixture, arg, args...); }, samples);
}


template <class Fixture, class ArgLatency, class ArgThroughput, class... Args>
MATHTER_NOINLINE void BenchmarkCase(std::string_view name, int64_t samples, int64_t repeat, Fixture&& fixture, ArgLatency&& argLatency, ArgThroughput&& argThroutput, Args&&... args) {
	const auto latency = Latency(samples, repeat, fixture, argLatency, args...);
	const auto throughput = Throughput(samples, repeat, fixture, argThroutput, args...);
	std::lock_guard lk{ g_mutex };
	g_records.push_back(BenchmarkRecord{ std::string(name), latency, throughput });
}


} // namespace impl


template <class Op, size_t Count, class Lhs, class... Args>
MATHTER_FORCEINLINE static auto DependentLoop(Op op,
											  const Lhs& lhs,
											  const std::array<Args, Count>&... args) {
	auto result = op(lhs, args[0]...);
	for (size_t i = 1; i < Count; ++i) {
		result = op(result, args[0]...);
	}
	return result;
}


template <class Op, size_t Lanes, size_t Count, class Lhs, class... Args>
MATHTER_FORCEINLINE static auto IndependentLoop(Op op,
												const std::array<Lhs, Lanes>& lhs,
												const std::array<Args, Count>&... args) {
	using Result = std::invoke_result_t<Op, Lhs, Args...>;
	std::array<Result, Lanes> result;
	for (size_t lane = 0; lane < Lanes; ++lane) {
		result[lane] = op(lhs[lane], args[lane]...);
	}
	for (size_t i = Lanes; i < Count; i += Lanes) {
		for (size_t lane = 0; lane < Lanes; ++lane) {
			result[lane] = op(result[lane], args[i + lane]...);
		}
	}
	return result;
}


#define BENCHMARK_CASE(NAME, TAG, SAMPLES, REPEAT, FIXTURE, ARG1_LATENCY, ARG1_THROUGHPUT, ...)            \
	TEST_CASE(NAME, TAG) {                                                                                 \
		::impl::BenchmarkCase(NAME, SAMPLES, REPEAT, FIXTURE, ARG1_LATENCY, ARG1_THROUGHPUT, __VA_ARGS__); \
	}
