/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include <string>
#include <iostream>
#include "MappingSettings.hpp"

#ifndef QMAP_MAPPINGRESULTS_HPP
#define QMAP_MAPPINGRESULTS_HPP

enum class Method {
	None, Exact, Heuristic
};
static std::string toString(const Method method) {
	switch (method) {
		case Method::None:
			return "none";
		case Method::Exact:
			return "exact";
		case Method::Heuristic:
			return "heuristic";
	}
	return " ";
}
// map Method values to JSON as strings
NLOHMANN_JSON_SERIALIZE_ENUM( Method, {
	{Method::None, "none"},
	{Method::Exact, "exact"},
	{Method::Heuristic, "heuristic"},
})

struct MappingResults {

	std::string input_name;
	unsigned long input_gates = 0;
	unsigned long input_singlequbitgates = 0;
	unsigned long input_cnots = 0;
	unsigned long input_qubits = 0;
	unsigned long input_layers = 0;

	std::string architecture;
	std::string calibration;

	Method method = Method::None;
	InitialLayoutStrategy initialLayoutStrategy = InitialLayoutStrategy::None;
	LayeringStrategy layeringStrategy = LayeringStrategy::None;

	double time = 0.0;
	bool timeout = true;
	std::string output_name;
	unsigned long output_gates = 0;
	unsigned long output_singlequbitgates = 0;
	unsigned long output_cnots = 0;
	unsigned long output_swaps = 0;
	unsigned long output_direction_reverse = 0;
	unsigned long output_qubits = 0;

	MappingResults() = default;
	virtual ~MappingResults() = default;

	virtual void copyInput(const MappingResults& mappingResults) {
		input_name = mappingResults.input_name;
		input_gates = mappingResults.input_gates;
		input_singlequbitgates = mappingResults.input_singlequbitgates;
		input_cnots = mappingResults.input_cnots;
		input_qubits = mappingResults.input_qubits;
		input_layers = mappingResults.input_layers;

		architecture = mappingResults.architecture;
		calibration = mappingResults.calibration;
		method = mappingResults.method;
		initialLayoutStrategy = mappingResults.initialLayoutStrategy;
		layeringStrategy = mappingResults.layeringStrategy;

		output_name = mappingResults.output_name;
		output_qubits = mappingResults.output_qubits;
	}

	virtual std::ostream& print(std::ostream& out, bool printStatistics) {
		out << "{\n";
		out << "\t\"circuit\": {\n";
		out << "\t\t\"name\": \"" << input_name << "\",\n";
		out << "\t\t\"qubits\": " << input_qubits << ",\n";
		out << "\t\t\"gates\": " << input_gates << ",\n";
		out << "\t\t\"singlequbitgates\": " << input_singlequbitgates << ",\n";
		out << "\t\t\"cnots\": " << input_cnots << ",\n";
		out << "\t\t\"layers\": " << input_layers << "\n";
		out << "\t},\n";
		out << "\t\"mapped_circuit\": {\n";
		out << "\t\t\"name\": \"" << output_name << "\",\n";
		out << "\t\t\"qubits\": " << output_qubits << ",\n";
		out << "\t\t\"gates\": " << output_gates << ",\n";
		out << "\t\t\"singlequbitgates\": " << output_singlequbitgates << ",\n";
		out << "\t\t\"cnots\": " << output_cnots << ",\n";
		out << "\t\t\"swaps\": " << output_swaps << ",\n";
		out << "\t\t\"direction_reverse\": " << output_direction_reverse << "\n";
		out << "\t}";
		if (printStatistics) {
			out << ",\n\t\"statistics\": {\n";
			out << "\t\t\"mapping_time\": " << (timeout? "\"timeout\"": std::to_string(time)) << ",\n";
			out << "\t\t\"additional_gates\": " << output_gates-input_gates << ",\n";
			out << "\t\t\"method\": \"" << toString(method) << "\",\n";

			if (layeringStrategy != LayeringStrategy::None) {
				out << "\t\t\"layeringStrategy\": \"" << toString(layeringStrategy) << "\",\n";
			}
			if (initialLayoutStrategy != InitialLayoutStrategy::None) {
				out << "\t\t\"initialLayoutStrategy\": \"" << toString(initialLayoutStrategy) << "\",\n";

			}
			out << "\t\t\"arch\": \"" << architecture << "\"";
			if (!calibration.empty()) {
				out << ",\n\t\t\"calibration\": \"" << calibration << "\"";
			}
			out << "\n\t}";
		}
		out << "\n}\n";

		return out;
	};

	virtual nlohmann::json produceJSON(bool statistics) {
		nlohmann::json resultJSON{};
		resultJSON["circuit"] = {};
		auto& circuit = resultJSON["circuit"];
		circuit["name"] = input_name;
		circuit["n_qubits"] = input_qubits;
		circuit["n_gates"] = input_gates;
		circuit["singlequbitgates"] = input_singlequbitgates;
		circuit["cnots"] = input_cnots;
		circuit["layers"] = input_layers;

		resultJSON["mapped_circuit"] = {};
		auto& mapped_circuit = resultJSON["mapped_circuit"];
		mapped_circuit["name"] = output_name;
		mapped_circuit["n_qubits"] = output_qubits;
		mapped_circuit["n_gates"] = output_gates;
		mapped_circuit["singlequbitgates"] = output_singlequbitgates;
		mapped_circuit["cnots"] = output_cnots;
		mapped_circuit["swaps"] = output_swaps;
		mapped_circuit["direction_reverse"] = output_direction_reverse;

		if (statistics) {
			resultJSON["statistics"] = {};
			auto& stats = resultJSON["statistics"];
			if (timeout)
				stats["timeout"] = timeout;
			stats["mapping_time"] = time;
			stats["additional_gates"] = output_gates-input_gates;
			stats["method"] = method;
			if (layeringStrategy != LayeringStrategy::None) {
				stats["layeringStrategy"] = layeringStrategy;
			}
			if (initialLayoutStrategy != InitialLayoutStrategy::None) {
				stats["initialLayoutStrategy"] = initialLayoutStrategy;
			}
			stats["arch"] = architecture;
			if (!calibration.empty())
				stats["calibration"] = calibration;
		}

		return resultJSON;
	}

	virtual std::string produceCSVEntry() {
		std::stringstream ss{};
		ss << input_name << ";" << input_qubits << ";" << input_gates << ";" << output_name << ";" << output_qubits << ";" << output_gates << ";" << toString(method) << ";";
		if (timeout) {
			ss << "TO";
		} else {
			ss << time;
		}
		ss << ";" << toString(initialLayoutStrategy) << ";" << toString(layeringStrategy);
		return ss.str();
	}
};

#endif //QMAP_MAPPINGRESULTS_HPP
