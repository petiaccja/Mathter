#pragma once

#include <Mathter/Common/OptimizationUtil.hpp>

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

	auto [result, count] = fixture.Throughput(std::forward<FirstArg>(arg), args...);

	for (int64_t i = 0; i < repeat; ++i) {
		auto [result_, count_] = fixture.Throughput(result[i % result.size()], args...);
		count += count_;
		result = result_;
	}
	DoNotOptimizeAway(result);

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


template <class Fixture, class Arg, class... Args>
void BenchmarkCase(std::string_view name, int64_t samples, int64_t repeat, Fixture&& fixture, Arg&& arg, Args&&... args) {
	const auto latency = Latency(samples, repeat, fixture, arg, args...);
	const auto throughput = Throughput(samples, repeat, fixture, arg, args...);
	std::lock_guard lk{ g_mutex };
	g_records.push_back(BenchmarkRecord{ std::string(name), latency, throughput });
}


} // namespace impl


template <class Lhs, class Rhs, class Op, size_t Count, size_t Index = 0>
MATHTER_FORCEINLINE static auto DependentUnroll(const Lhs& lhs,
												const std::array<Rhs, Count>& rhs,
												Op op,
												std::integral_constant<size_t, Index> = std::integral_constant<size_t, 0>{}) {
	if constexpr (Index < Count) {
		return DependentUnroll(op(lhs, rhs[Index]), rhs, op, std::integral_constant<size_t, Index + 1>{});
	}
	else {
		return lhs;
	}
}


template <class Out, class Lhs, class Rhs, class Op, size_t Count, size_t Index = 0>
MATHTER_FORCEINLINE static auto IndependentUnroll(std::array<Out, Count>& out,
												  const Lhs& lhs,
												  const std::array<Rhs, Count>& rhs,
												  Op op,
												  std::integral_constant<size_t, Index> = std::integral_constant<size_t, 0>{}) {
	if constexpr (Index < Count) {
		out[Index] = op(lhs, rhs[Index]);
		IndependentUnroll(out, lhs, rhs, op, std::integral_constant<size_t, Index + 1>{});
	}
}


template <class Lhs, class Rhs, class Op, size_t Count, size_t Index = 0>
MATHTER_FORCEINLINE static auto IndependentUnroll(const Lhs& lhs,
												  const std::array<Rhs, Count>& rhs,
												  Op op,
												  std::integral_constant<size_t, Index> = std::integral_constant<size_t, 0>{}) {
	using R = std::invoke_result_t<Op, Lhs, Rhs>;
	std::array<R, Count> out;
	IndependentUnroll(out, lhs, rhs, op, std::integral_constant<size_t, Index + 1>{});
	return out;
}


#define BENCHMARK_CASE(NAME, TAG, SAMPLES, REPEAT, FIXTURE, ARG1, ...)            \
	TEST_CASE(NAME, TAG) {                                                        \
		::impl::BenchmarkCase(NAME, SAMPLES, REPEAT, FIXTURE, ARG1, __VA_ARGS__); \
	}