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
};

#endif //QMAP_MAPPINGRESULTS_HPP
