#pragma once

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <intrin.h> // Not portable
#include <mutex>
#include <tuple>
#include <vector>


namespace impl {


struct BenchmarkRecord {
	std::string name;
	float latency;
	float throughput;
};


extern std::vector<BenchmarkRecord> g_records;
extern std::mutex g_mutex;


inline int64_t ReadTSC() {
	unsigned aux;
	return __rdtscp(&aux);
}


template <class T>
void DoNotOptimizeAway(T&& value) {
	// Idea from celero https://github.com/DigitalInBlue/Celero/blob/master/include/celero/Utilities.h
	if (ReadTSC() == 0) {
		volatile thread_local char c;
		c = *reinterpret_cast<const char*>(&value);
	}
}


template <class DoSample>
int64_t BestSample(DoSample&& doSample, int64_t samples) {
	int64_t bestTime = std::numeric_limits<int64_t>::max();
	for (size_t i = 0; i < samples; ++i) {
		bestTime = std::min(bestTime, doSample());
	}
	return bestTime;
}


template <class Operation, class Feed, class... Init>
float Latency(int64_t samples, int64_t repeat, Operation&& operation, Feed&& feed, const std::tuple<Init...>& init) {
	auto doSample = [&]() {
		auto v = init;

		const auto startTime = ReadTSC();
		for (int64_t i = 0; i < repeat; ++i) {
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);

			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);

			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);

			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
			v = feed(operation(v), v);
		}
		DoNotOptimizeAway(v);

		const auto endTime = ReadTSC();
		return endTime - startTime;
	};

	return BestSample(doSample, samples) / float(repeat * 20);
}


template <class Operation, class... Init, size_t N>
float Throughput(int64_t samples, int64_t repeat, Operation&& operation, const std::array<std::tuple<Init...>, N>& init) {
	auto doSample = [&]() {
		size_t index = 0;
		auto r = operation(init[0]);
		auto v = init[0];

		const auto startTime = ReadTSC();
		for (int64_t i = 0; i < repeat; ++i) {
			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;


			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;


			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;


			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;

			v = init[index];
			r = operation(v);
			index = (index + 1) % N;
		}
		DoNotOptimizeAway(r);

		const auto endTime = ReadTSC();
		return endTime - startTime;
	};

	return BestSample(doSample, samples) / float(repeat * 20);
}


template <class Operation, class Feed, class... Init, size_t N>
void BenchmarkCase(std::string_view name, int64_t samples, int64_t repeat, Operation&& operation, Feed&& feed, const std::array<std::tuple<Init...>, N>& init) {
	const auto latency = Latency(samples, repeat, operation, feed, init[0]);
	const auto throughput = Throughput(samples, repeat, operation, init);
	std::lock_guard lk{ g_mutex };
	g_records.push_back(BenchmarkRecord{ std::string(name), latency, throughput });
}


} // namespace impl


#define BENCHMARK_CASE(NAME, TAG, SAMPLES, REPEAT, OPERATION, FEED, INIT)    \
	TEST_CASE(NAME, TAG) {                                                   \
		::impl::BenchmarkCase(NAME, SAMPLES, REPEAT, OPERATION, FEED, INIT); \
	}