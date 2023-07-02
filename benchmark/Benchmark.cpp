#include "Benchmark.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>


namespace impl {

std::vector<BenchmarkRecord> g_records;
std::mutex g_mutex;

} // namespace impl


namespace {

#ifdef MATHTER_TSC_USES_CHRONO
#define MATHTER_TIME_MEASURE "ns"
#else
#define MATHTER_TIME_MEASURE "cycles"
#endif


void PrintCases() {
	constexpr std::string_view headerName = "Name";
	constexpr std::string_view headerLatency = "Latency (" MATHTER_TIME_MEASURE ")";
	constexpr std::string_view headerThroughput = "Throughput (" MATHTER_TIME_MEASURE ")";

	std::lock_guard lkg{ impl::g_mutex };

	const auto& records = impl::g_records;

	std::vector<std::string_view> names;
	std::vector<std::string> latencies;
	std::vector<std::string> throughputs;
	std::transform(records.begin(), records.end(), std::back_inserter(names), [](const auto& v) { return std::string_view(v.name); });
	std::transform(records.begin(), records.end(), std::back_inserter(latencies), [](const auto& v) { return std::to_string(v.latency); });
	std::transform(records.begin(), records.end(), std::back_inserter(throughputs), [](const auto& v) { return std::to_string(v.throughput); });

	auto sizeFun = [](const auto& v) { return std::size(v); };
	auto maxFun = [](const auto& a, const auto& b) { return std::max(a, b); };
	std::array colSizes = {
		std::transform_reduce(names.begin(), names.end(), std::size(headerName), maxFun, sizeFun),
		std::transform_reduce(names.begin(), names.end(), std::size(headerLatency), maxFun, sizeFun),
		std::transform_reduce(names.begin(), names.end(), std::size(headerThroughput), maxFun, sizeFun),
	};

	auto makeLine = [&](char fill, char column) {
		auto lineLength = std::reduce(colSizes.begin(), colSizes.end()) + 3 * colSizes.size() + 1;
		std::string line(lineLength, fill);
		auto offset = 0;
		line[0] = column;
		for (const auto& colSize : colSizes) {
			offset += colSize;
			offset += 3;
			line[offset] = column;
		}
		return line;
	};

	auto printSeparator = [&]() {
		std::cout << makeLine('-', '+') << std::endl;
	};

	auto printRecord = [&](std::string_view name, std::string_view latency, std::string_view throughput) {
		auto line = makeLine(' ', '|');
		auto columnBegin = std::find(line.rbegin(), line.rend(), '|');
		std::copy(throughput.begin(), throughput.end(), columnBegin.base() - throughput.size() - 2);
		columnBegin = std::find(columnBegin + 1, line.rend(), '|');
		std::copy(latency.begin(), latency.end(), columnBegin.base() - latency.size() - 2);
		columnBegin = std::find(columnBegin + 1, line.rend(), '|');
		columnBegin = std::find(columnBegin + 1, line.rend(), '|');
		std::copy(name.begin(), name.end(), columnBegin.base() + 1);
		std::cout << line << std::endl;
	};

	printSeparator();
	printRecord(headerName, headerLatency, headerThroughput);
	printSeparator();
	for (size_t i = 0; i < records.size(); ++i) {
		printRecord(names[i], latencies[i], throughputs[i]);
	}
	printSeparator();
}

volatile bool forcePrinting = [] {
	std::atexit([] { PrintCases(); });
	return true;
}();

} // namespace