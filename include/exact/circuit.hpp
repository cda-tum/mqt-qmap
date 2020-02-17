/*
Minimal Mapping of Quantum Circuits to IBM QX Architectures by JKU Linz, Austria

Developer: Robert Wille, Lukas Burgholzer, Alwin Zulehner

For more information, please visit http://iic.jku.at/eda/research/ibm_qx_mapping

If you have any questions feel free to contact us using
robert.wille@jku.at, lukas.burgholzer@jku.at or alwin.zulehner@jku.at

If you use the compiler for your research, we would be thankful if you referred to it
by citing the following publication:

@inproceedings{wille2019mapping,
    title={Mapping Quantum Circuits to {IBM QX} Architectures Using the Minimal Number of {SWAP} and {H} Operations},
    author={Wille, Robert and Burgholzer, Lukas and Zulehner, Alwin},
    booktitle={Design Automation Conference},
    year={2019}
}
*/

#ifndef Circuit_hpp
#define Circuit_hpp

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <chrono>

#include "z3++.h"
#include "parser.hpp"

using namespace std::chrono;
using namespace z3;

typedef std::set<std::pair<int,int>> CouplingMap;

/// Create a string representation of a given permutation
/// \param pi permutation
/// \return string representation of pi
std::string printPi(std::vector<int>& pi);

/// Iterating routine through all combinations
/// \tparam Iterator iterator type
/// \param first iterator to beginning
/// \param k current iterator
/// \param last iterator to end
/// \return true if another combination was found
template<typename Iterator> bool next_combination(Iterator first, Iterator k, Iterator last);

/// Simple depth-first-search implementation used to check whether a given subset of qubits is
/// connected on the given architecture
/// \param current index of current qubit
/// \param visited visited qubits
/// \param cm coupling map of architecture
void dfs(int current, std::set<int>& visited, CouplingMap& cm);

/// Structure encapsulating the mapping results
struct MappingResults {
	bool timeout;
	std::vector<std::vector<std::vector<int>>> X;
	std::vector<std::vector<int>> Y;
	std::vector<int> Z;
	unsigned long swapCost;
	unsigned long reverseCost;
	unsigned long totalCost;
	std::set<int> logicalQubits;
	std::set<int> usedPhysicalQubits;
	unsigned long nrGatesOriginalCurcuit;
	unsigned long nrGatesReducedCurcuit;
	unsigned long nrGatesMappedCurcuit;
	unsigned long nrIgnoredUnaryGates;
	unsigned long nrAllGatesMappedCurcuit;
	unsigned long nrLayersReducedCurcuit;

	MappingResults() {
		timeout = true;
		X = std::vector<std::vector<std::vector<int>>>();
		Y = std::vector<std::vector<int>>();
		Z = std::vector<int>();
		swapCost = ULONG_MAX;
		reverseCost = ULONG_MAX;
		totalCost = ULONG_MAX;
		logicalQubits = std::set<int>();
		usedPhysicalQubits = std::set<int>();
		nrGatesOriginalCurcuit = 0;
		nrGatesMappedCurcuit = 0;
		nrGatesReducedCurcuit = 0;
		nrIgnoredUnaryGates = 0;
		nrAllGatesMappedCurcuit = 0;
		nrLayersReducedCurcuit = 0;
	}

	/// Print results
	/// \param fullOutput If true, also output the assignment of all mapping variables
	void print(bool fullOutput = false) {

		if (timeout) {
			std::cout << "##############################################" << std::endl;
            std::cout << "nrGatesOriginalCurcuit: " << nrGatesOriginalCurcuit << std::endl;
            std::cout << "nrGatesReducedCurcuit: " << nrGatesReducedCurcuit << std::endl;
            std::cout << "nrLayersReducedCurcuit: " << nrLayersReducedCurcuit << std::endl;
            std::cout << "nrIgnoredUnaryGates: " << nrIgnoredUnaryGates << std::endl;
            std::cout << "nrLogicalQubits: " << logicalQubits.size() << std::endl;
            std::cout << "logicalQubits: ";
            for(auto q: logicalQubits)
                std::cout << q << " ";
            std::cout << "\n";
            std::cout << "nrUsedPhysicalQubits: " << usedPhysicalQubits.size() << std::endl;
            std::cout << "usedPhysicalQubits: ";
            for(auto q: usedPhysicalQubits)
                std::cout << q << " ";
            std::cout << "\n";
			std::cout << "##############################################" << std::endl;
			std::cout << "timeout:" << std::endl;
			return;
		}
		std::cout << "##############################################" << std::endl;
        std::cout << "nrGatesOriginalCurcuit: " << nrGatesOriginalCurcuit << std::endl;
		std::cout << "nrGatesReducedCurcuit: " << nrGatesReducedCurcuit << std::endl;
        std::cout << "nrLayersReducedCurcuit: " << nrLayersReducedCurcuit << std::endl;
        std::cout << "nrIgnoredUnaryGates: " << nrIgnoredUnaryGates << std::endl;
        std::cout << "nrLogicalQubits: " << logicalQubits.size() << std::endl;
        std::cout << "logicalQubits: ";
		for(auto q: logicalQubits)
			std::cout << q << " ";
		std::cout << "\n";
        std::cout << "nrUsedPhysicalQubits: " << usedPhysicalQubits.size() << std::endl;
        std::cout << "usedPhysicalQubits: ";
		for(auto q: usedPhysicalQubits)
			std::cout << q << " ";
		std::cout << "\n";

		std::cout << "swapCost: " << swapCost << std::endl;
		std::cout << "reverseCost: " << reverseCost << std::endl;
		std::cout << "totalCost: " << totalCost << std::endl;
		std::cout << "nrGatesMappedCurcuit: " << nrGatesMappedCurcuit << std::endl;
        std::cout << "nrAllGatesMappedCurcuit: " << nrAllGatesMappedCurcuit << std::endl;
        if (fullOutput) {
			std::cout << "----------------------------------------------" << std::endl;
			for(auto& gate: X) {
				for(auto& row: gate) {
					for (auto& col: row) {
						std::cout << col;
					}
					std::cout << "\n";
				}
				std::cout << "\n";
			}
			std::cout << "----------------------------------------------" << std::endl;
			for(auto& gate: Y) {
				std::cout << printPi(gate) << "\n";
			}
			std::cout << "----------------------------------------------" << std::endl;
			for(auto& gate: Z) {
				std::cout << gate << "\n";
			}
		}
		std::cout << "##############################################" << std::endl;
	}

	/// Print minimal mapping results
	void printOptimum() {
		if (timeout) {
			std::cout << "Optimum:######################################" << std::endl;
            std::cout << "nrGatesOriginalCurcuit: " << nrGatesOriginalCurcuit << std::endl;
            std::cout << "nrGatesReducedCurcuit: " << nrGatesReducedCurcuit << std::endl;
            std::cout << "nrLayersReducedCurcuit: " << nrLayersReducedCurcuit << std::endl;
            std::cout << "nrIgnoredUnaryGates: " << nrIgnoredUnaryGates << std::endl;
            std::cout << "nrLogicalQubits: " << logicalQubits.size() << std::endl;
            std::cout << "logicalQubits: ";
            for(auto q: logicalQubits)
                std::cout << q << " ";
            std::cout << "\n";
            std::cout << "nrUsedPhysicalQubits: " << usedPhysicalQubits.size() << std::endl;
            std::cout << "usedPhysicalQubits: ";
            for(auto q: usedPhysicalQubits)
                std::cout << q << " ";
            std::cout << "\n";
			std::cout << "timeout:" << std::endl;
			std::cout << "##############################################" << std::endl;
			return;
		}
        std::cout << "Optimum:######################################" << std::endl;
        std::cout << "nrGatesOriginalCurcuit: " << nrGatesOriginalCurcuit << std::endl;
        std::cout << "nrGatesReducedCurcuit: " << nrGatesReducedCurcuit << std::endl;
        std::cout << "nrLayersReducedCurcuit: " << nrLayersReducedCurcuit << std::endl;
        std::cout << "nrIgnoredUnaryGates: " << nrIgnoredUnaryGates << std::endl;
        std::cout << "nrLogicalQubits: " << logicalQubits.size() << std::endl;
        std::cout << "logicalQubits: ";
        for(auto q: logicalQubits)
            std::cout << q << " ";
        std::cout << "\n";
        std::cout << "nrUsedPhysicalQubits: " << usedPhysicalQubits.size() << std::endl;
        std::cout << "usedPhysicalQubits: ";
        for(auto q: usedPhysicalQubits)
            std::cout << q << " ";
        std::cout << "\n";

        std::cout << "swapCost: " << swapCost << std::endl;
        std::cout << "reverseCost: " << reverseCost << std::endl;
        std::cout << "totalCost: " << totalCost << std::endl;
        std::cout << "nrGatesMappedCurcuit: " << nrGatesMappedCurcuit << std::endl;
        std::cout << "nrAllGatesMappedCurcuit: " << nrAllGatesMappedCurcuit << std::endl;
		std::cout << "----------------------------------------------" << std::endl;
		for (auto& gate: X) {
			for (auto& row: gate) {
				for (auto& col: row) {
					std::cout << col;
				}
				std::cout << "\n";
			}
			std::cout << "\n";
		}
		std::cout << "----------------------------------------------" << std::endl;
		for (auto& gate: Y) {
			std::cout << printPi(gate) << "\n";
		}
		std::cout << "----------------------------------------------" << std::endl;
		for (auto& gate: Z) {
			std::cout << gate << "\n";
		}
		std::cout << "##############################################" << std::endl;
	}

};

/// Structure encapsulating the strategies used during the mapping
struct MappingSettings {
    // Using subsets of qubits
    bool useMinimumSetOfQubits = false;

    // Clustering
    bool exactStrategy = true;
    bool disjointQubitsStrategy = false;
    bool oddGatesStrategy = false;
    bool qubitTriangleStrategy = false;

    MappingSettings() {};
    void considerQubitSubsets() { useMinimumSetOfQubits=true; }
    void considerAllQubits() { useMinimumSetOfQubits=false; }
    void useExactStrategy() {
	    exactStrategy = true;
	    disjointQubitsStrategy = false;
	    oddGatesStrategy = false;
	    qubitTriangleStrategy = false;
    }
    void useDisjointQubitsStrategy() {
    	exactStrategy=false;
		disjointQubitsStrategy = true;
	    oddGatesStrategy = false;
	    qubitTriangleStrategy = false;
    }
    void useOddGatesStrategy() {
	    exactStrategy = false;
	    disjointQubitsStrategy = false;
	    oddGatesStrategy = true;
	    qubitTriangleStrategy = false;
    }
    void useQubitTriangleStrategy() {
	    exactStrategy = false;
	    disjointQubitsStrategy = false;
	    oddGatesStrategy = false;
	    qubitTriangleStrategy = true;
    }


};

/// Main structure representing the circuit and mapping functionality
struct Circuit {
	unsigned long nrLogicalQubits;
	unsigned long nrUsedPhysicalQubits;
	unsigned long nrPhysicalQubits;
	unsigned long nrGates;
	unsigned long nrLayers;
	std::vector<QASMparser::gate> gates;
	std::vector<std::vector<QASMparser::gate>> layers;
	std::set<int> logicalQubits;
	std::set<int> usedPhysicalQubits;
	std::vector<int> physicalQubits;
	unsigned int timeout;

	/// Constructor
	/// \param logicalQubits Logical qubits to use
	/// \param usedPhysicalQubits Physical qubits to use
	/// \param physicalQubits  All available physical qubits
	/// \param layers Gate groupings
	/// \param timeout Timeout threshold
    Circuit(std::set<int>& logicalQubits, std::set<int>& usedPhysicalQubits, std::vector<int>& physicalQubits, std::vector<std::vector<QASMparser::gate>>& layers, unsigned int timeout) {
        this->nrLogicalQubits = logicalQubits.size();
        this->nrPhysicalQubits = physicalQubits.size();
        this->nrUsedPhysicalQubits = usedPhysicalQubits.size();
        this->logicalQubits = logicalQubits;
        this->usedPhysicalQubits = usedPhysicalQubits;
        this->physicalQubits = physicalQubits;
        this->timeout = timeout;
        this->layers = layers;
        this->nrLayers = layers.size();

        this->gates = std::vector<QASMparser::gate>();
        this->nrGates = 0;
        for(auto& layer: layers) {
            for (auto& gate: layer) {
                gates.push_back(gate);
            }
        }
        nrGates = gates.size();
    }

    /// Static driver routine
    /// \param filename QASM-file to read circuit from
    /// \param timeout Timeout threshold
    /// \param cm Coupling map
    /// \param physicalQubits Available physical qubits
    /// \param cost Cost function
    /// \param settings Mapping settings
    /// \return Mapping results
    static MappingResults run(std::string& filename, unsigned int timeout, const CouplingMap& cm, std::vector<int>& physicalQubits, const std::function< unsigned int(std::vector<int>&)>& cost, MappingSettings settings);

    /// Core mapping routine
    /// \param cm Coupling map
    /// \param cost Cost function
    /// \return Mapping results
	MappingResults mapping(const CouplingMap& cm, const std::function< unsigned int(std::vector<int>&)>& cost);

	/// Helper function returning correct index in 1D array for specific gate and logical/physical qubit
	/// \param k gate index
	/// \param i physical qubit index
	/// \param j logical qubit index
	/// \return index in 1D array
    inline
    unsigned long idx(int k, int i, int j) {
        int counti = 0;
        for (int usedPhysicalQubit : usedPhysicalQubits) {
            if (usedPhysicalQubit == i) break;
            counti++;
        }
        int countj = 0;
        for (int logicalQubit : logicalQubits) {
            if (logicalQubit == j) break;
            countj++;
        }

        return k*nrLogicalQubits*nrUsedPhysicalQubits + counti*nrLogicalQubits + countj;
    }

    /// Computes n! recursively
    /// \param n interger to compute factorial of
    /// \return n!
	static inline
    unsigned long factorial(unsigned long n) {
        if (n == 1)
            return 1;
        else
            return n * factorial(n - 1);
    }
};

int exact_mapping(std::string filename);

#endif /* Circuit_hpp */