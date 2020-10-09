//
// Created by Lukas Burgholzer on 28.02.20.
//


#include <string>
#include <iostream>
#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;

#ifndef QMAP_MAPPINGRESULTS_HPP
#define QMAP_MAPPINGRESULTS_HPP

enum Method {
	None, Exact, Heuristic
};

static std::string toString(const Method method) {
	switch (method) {
		case None:
			return "none";
		case Exact:
			return "exact";
		case Heuristic:
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

	std::string architecture = "";
	std::string calibration = "";

	Method method = None;

	double time = 0.0;
	bool timeout = true;
	std::string output_name = "";
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

		output_name = mappingResults.output_name;
		output_qubits = mappingResults.output_qubits;
	}

	virtual std::ostream& print(std::ostream& out) {
		json j;

		j["circuit"]["name"] = input_name;
		j["circuit"]["qubits"] = input_qubits;
		j["circuit"]["gates"] = input_gates;
		j["circuit"]["singlequbitgates"] = input_singlequbitgates;
		j["circuit"]["cnots"] = input_cnots;
		j["circuit"]["layers"] = input_layers;

		j["mapped_circuit"]["name"] = output_name;
		j["mapped_circuit"]["qubits"] = output_qubits;
		j["mapped_circuit"]["gates"] = output_gates;
		j["mapped_circuit"]["singlequbitgates"] = output_singlequbitgates;
		j["mapped_circuit"]["cnots"] = output_cnots;
		j["mapped_circuit"]["swaps"] = output_swaps;
		j["mapped_circuit"]["direction_reverse"] = output_direction_reverse;

		if (timeout) {
			j["mapping_time"] = "timeout";
		} else {
			j["mapping_time"] = time;
		}
		j["method"] = toString(method);
		j["arch"] = architecture;
		j["calibration"] = calibration;

		out << j.dump(4) << "\n";

		return out;
	};
};

#endif //QMAP_MAPPINGRESULTS_HPP
