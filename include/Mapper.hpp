//
// Created by Lukas Burgholzer on 28.02.20.
//

#ifndef QMAP_MAPPER_HPP
#define QMAP_MAPPER_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <array>
#include <unordered_set>
#include <chrono>

#include "QuantumComputation.hpp"
#include "Architecture.hpp"
#include "MappingResults.hpp"
#include "MappingSettings.hpp"

constexpr short DEFAULT_POSITION = -1;
constexpr double INITIAL_FIDELITY = 1.0;
constexpr unsigned short MAX_DEVICE_QUBITS = 128;

class Mapper {
protected:
	Architecture architecture = {};

	// internal structures
	struct Gate {
		short control = -1;
		unsigned short target = 0;

		qc::Operation* op = nullptr;

		Gate (short c, unsigned short t): control(c), target(t) {};
		Gate (short c, unsigned short t, qc::Operation* op): control(c), target(t), op(op) {};

		bool singleQubit() const {
			return control == -1;
		}
	};

	qc::QuantumComputation qc;
	qc::QuantumComputation qcMapped;
	std::vector<std::vector<Gate>> layers{};

	std::array<short, MAX_DEVICE_QUBITS> qubits{};
	std::array<short, MAX_DEVICE_QUBITS> locations{};
	std::array<double, MAX_DEVICE_QUBITS> fidelities{};

	std::unordered_set<unsigned short> usedDeviceQubits{};

	MappingResults results{};
	MappingSettings settings{};

	explicit Mapper(const std::string& filename);
	//explicit Mapper(qc::QuantumComputation& qc);

	virtual void initResults();

	virtual void createLayers();

	virtual long getNextLayer(size_t idx);

public:
	Mapper(const std::string& filename, const std::string& cm_filename);
	Mapper(const std::string& filename, const std::string& cm_filename, const std::string& cal_filename);

	Mapper(const std::string& filename, unsigned short nQ, const CouplingMap& couplingMap);
	Mapper(const std::string& filename, unsigned short nQ, const CouplingMap& couplingMap, const std::vector<Architecture::CalibrationData>& calibrationData);

	//Mapper(qc::QuantumComputation& qc, const std::string& cm_filename);
	//Mapper(qc::QuantumComputation& qc, const std::string& cm_filename, const std::string& cal_filename);

	//Mapper(qc::QuantumComputation& qc, unsigned short nQ, const CouplingMap& couplingMap);
	//Mapper(qc::QuantumComputation& qc, unsigned short nQ, const CouplingMap& couplingMap, const std::vector<Architecture::CalibrationData>& calibrationData);

	virtual void map(const MappingSettings& ms) = 0;

	virtual void dumpResult(const std::string& outputFilename) {
		if (qcMapped.empty()) {
			std::cerr << "Mapped circuit is empty." << std::endl;
			return;
		}

		size_t dot = outputFilename.find_last_of('.');
		std::string extension = outputFilename.substr(dot + 1);
		std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) { return ::tolower(c); });
		if (extension == "real") {
			dumpResult(outputFilename, qc::Real);
		} else if (extension == "qasm") {
			dumpResult(outputFilename, qc::OpenQASM);
		} else {
			throw QMAPException("[dump] Extension " + extension + " not recognized/supported for dumping.");
		}

	}
	virtual void dumpResult(const std::string& outputFilename, qc::Format format) {
		size_t slash = outputFilename.find_last_of('/');
		size_t dot = outputFilename.find_last_of('.');
		results.output_name = outputFilename.substr(slash+1, dot-slash-1);
		qcMapped.dump(outputFilename, format);
	}


	virtual std::ostream& printResult(std::ostream& out) {
		return results.print(out);
	}

	std::ostream& printLayering(std::ostream& out) {
		out << "---------------- Layering -------------------" << std::endl;
		for(auto& layer: layers) {
			for (auto& gate: layer) {
				if (gate.singleQubit()) {
					out << "(" << gate.target << ") ";
				} else {
					out << "(" << gate.control << " " << gate.target << ") ";
				}
			}
			out << std::endl;
		}
		out << "---------------------------------------------" << std::endl;
		return out;
	}

	std::ostream& printLocations(std::ostream& out) {
		out << "---------------- Locations -------------------" << std::endl;
		for (unsigned short i = 0; i<qc.getNqubits(); ++i) {
			out << locations.at(i) << " ";
		}
		out << std::endl;
		out << "---------------------------------------------" << std::endl;
		return out;
	}
	std::ostream& printQubits(std::ostream& out) {
		out << "---------------- Qubits -------------------" << std::endl;
		for (unsigned short i=0; i<architecture.getNqubits(); ++i) {
			out << qubits.at(i) << " ";
		}
		out << std::endl;
		out << "---------------------------------------------" << std::endl;
		return out;
	}

	virtual void reset() {
		architecture.reset();
		qc.reset();
		layers.clear();
		qubits.fill(DEFAULT_POSITION);
		locations.fill(DEFAULT_POSITION);
		usedDeviceQubits.clear();

		results = MappingResults();
		settings = MappingSettings();
	}

};

#endif //QMAP_MAPPER_HPP
