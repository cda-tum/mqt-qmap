//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "DataLogger.hpp"

#include "nlohmann/json.hpp"

#include <filesystem>

void DataLogger::initLog() {
  if (dataLoggingPath.back() != '/') {
    dataLoggingPath += '/';
  }
  const std::filesystem::path dirPath(dataLoggingPath);
  if (!std::filesystem::exists(dirPath)) {
    std::filesystem::create_directory(dirPath);
  }
  clearLog();
};

void DataLogger::clearLog() {
  for (const auto& entry :
       std::filesystem::directory_iterator(dataLoggingPath)) {
    std::filesystem::remove_all(entry.path());
  }
};

void DataLogger::logArchitecture() {
  if (deactivated) {
    return;
  }

  auto of = std::ofstream(dataLoggingPath + "architecture.json");
  if (!of.good()) {
    deactivated = true;
    std::cerr << "[data-logging] Error opening file: " << dataLoggingPath
              << "architecture.json" << std::endl;
    return;
  }
  nlohmann::json json;
  json["name"]         = architecture->getName();
  json["nqubits"]      = architecture->getNqubits();
  json["coupling_map"] = architecture->getCouplingMap();
  json["distances"]    = architecture->getDistanceTable();
  if (architecture->isFidelityAvailable()) {
    auto& fidelity = json["fidelity"];
    fidelity["single_qubit_fidelities"] =
        architecture->getSingleQubitFidelities();
    fidelity["two_qubit_fidelities"] = architecture->getFidelityTable();
    fidelity["single_qubit_fidelity_costs"] =
        architecture->getSingleQubitFidelityCosts();
    fidelity["two_qubit_fidelity_costs"] =
        architecture->getTwoQubitFidelityCosts();
    fidelity["swap_fidelity_costs"] = architecture->getSwapFidelityCosts();
    fidelity["fidelity_distances"]  = architecture->getFidelityDistanceTables();
  }
  of << json.dump(2);
  of.close();
};

void DataLogger::openNewLayer(std::size_t layerIndex) {
  if (deactivated) {
    return;
  }

  for (std::size_t i = searchNodesLogFiles.size(); i <= layerIndex; ++i) {
    searchNodesLogFiles.emplace_back(dataLoggingPath + "nodes_layer_" +
                                     std::to_string(i) + ".csv");
    if (!searchNodesLogFiles.at(i).good()) {
      deactivated = true;
      std::cerr << "[data-logging] Error opening file: " << dataLoggingPath
                << "layer_" << i << ".json" << std::endl;
      return;
    }
  }
};

void DataLogger::logFinalizeLayer(
    std::size_t layerIndex, const qc::CompoundOperation& ops,
    const std::vector<std::uint16_t>& singleQubitMultiplicity,
    const std::map<std::pair<std::uint16_t, std::uint16_t>,
                   std::pair<std::uint16_t, std::uint16_t>>&
                                                       twoQubitMultiplicity,
    const std::array<std::int16_t, MAX_DEVICE_QUBITS>& initialLayout,
    std::size_t finalNodeId, double finalCostFixed, double finalCostHeur,
    double                                             finalLookaheadPenalty,
    const std::array<std::int16_t, MAX_DEVICE_QUBITS>& finalLayout,
    const std::vector<Exchange>& finalSwaps, std::size_t finalSearchDepth) {
  if (deactivated) {
    return;
  }

  if (!searchNodesLogFiles.at(layerIndex).is_open()) {
    std::cerr << "[data-logging] Error: layer " << layerIndex
              << " has already been finalized" << std::endl;
    return;
  }
  searchNodesLogFiles.at(layerIndex).close();

  auto of = std::ofstream(dataLoggingPath + "layer_" +
                          std::to_string(layerIndex) + ".json");
  if (!of.good()) {
    deactivated = true;
    std::cerr << "[data-logging] Error opening file: " << dataLoggingPath
              << "layer_" << layerIndex << ".json" << std::endl;
    return;
  }
  nlohmann::json    json;
  std::stringstream qasmStream;
  ops.dumpOpenQASM(qasmStream, qregs, cregs);
  json["qasm"] = qasmStream.str();
  if (twoQubitMultiplicity.empty()) {
    json["two_qubit_multiplicity"] = nlohmann::json::array();
  } else {
    auto&       twoMultJSON = json["two_qubit_multiplicity"];
    std::size_t j           = 0;
    for (const auto& [qubits, multiplicity] : twoQubitMultiplicity) {
      twoMultJSON[j]["q1"]       = qubits.first;
      twoMultJSON[j]["q2"]       = qubits.second;
      twoMultJSON[j]["forward"]  = multiplicity.first;
      twoMultJSON[j]["backward"] = multiplicity.second;
      ++j;
    }
  }
  json["single_qubit_multiplicity"] = singleQubitMultiplicity;
  auto& initialLayoutJSON           = json["initial_layout"];
  for (std::size_t i = 0; i < nqubits; ++i) {
    initialLayoutJSON[i] = initialLayout.at(i);
  }
  json["final_node_id"]           = finalNodeId;
  json["final_cost_fixed"]        = finalCostFixed;
  json["final_cost_heur"]         = finalCostHeur;
  json["final_lookahead_penalty"] = finalLookaheadPenalty;
  auto& finalLayoutJSON           = json["final_layout"];
  for (std::size_t i = 0; i < nqubits; ++i) {
    finalLayoutJSON[i] = finalLayout.at(i);
  }
  if (finalSwaps.empty()) {
    json["final_swaps"] = nlohmann::json::array();
  } else {
    auto&       finalSwapsJSON = json["final_swaps"];
    std::size_t i              = 0;
    for (const auto& swap : finalSwaps) {
      finalSwapsJSON[i][0] = swap.first;
      finalSwapsJSON[i][1] = swap.second;
      ++i;
    }
  }
  json["final_search_depth"] = finalSearchDepth;
  of << json.dump(2);
  of.close();
};

void DataLogger::splitLayer() {
  if (deactivated) {
    return;
  }

  const std::size_t layerIndex = searchNodesLogFiles.size() - 1;
  if (searchNodesLogFiles.at(layerIndex).is_open()) {
    std::cerr << "[data-logging] Error: layer " << layerIndex
              << " has not been finalized before splitting" << std::endl;
    return;
  }
  searchNodesLogFiles.pop_back();
  std::size_t splitIndex = 0;
  while (std::filesystem::exists(dataLoggingPath + "nodes_layer_" +
                                 std::to_string(layerIndex) + ".presplit-" +
                                 std::to_string(splitIndex) + ".csv")) {
    ++splitIndex;
  }
  std::filesystem::rename(
      dataLoggingPath + "nodes_layer_" + std::to_string(layerIndex) + ".csv",
      dataLoggingPath + "nodes_layer_" + std::to_string(layerIndex) +
          ".presplit-" + std::to_string(splitIndex) + ".csv");
  std::filesystem::rename(
      dataLoggingPath + "layer_" + std::to_string(layerIndex) + ".json",
      dataLoggingPath + "layer_" + std::to_string(layerIndex) + ".presplit-" +
          std::to_string(splitIndex) + ".json");
}

void DataLogger::logSearchNode(
    std::size_t layerIndex, std::size_t nodeId, std::size_t parentId,
    double costFixed, double costHeur, double lookaheadPenalty,
    const std::array<std::int16_t, MAX_DEVICE_QUBITS>& qubits,
    bool validMapping, const std::vector<Exchange>& swaps, std::size_t depth) {
  if (deactivated) {
    return;
  }

  if (layerIndex >= searchNodesLogFiles.size()) {
    openNewLayer(layerIndex);
  }

  auto& of = searchNodesLogFiles.at(layerIndex);
  if (!of.is_open()) {
    deactivated = true;
    std::cerr << "[data-logging] Error: layer " << layerIndex
              << " has already been finalized" << std::endl;
    return;
  }
  of << nodeId << ";" << parentId << ";" << costFixed << ";" << costHeur << ";"
     << lookaheadPenalty << ";" << validMapping << ";" << depth << ";";
  for (std::size_t i = 0; i < nqubits; ++i) {
    of << qubits.at(i) << ",";
  }
  if (nqubits > 0) {
    of.seekp(-1, std::ios_base::cur); // remove last comma
  }
  of << ";";
  for (const auto& s : swaps) {
    of << s.first << " " << s.second;
    if (s.op != qc::OpType::SWAP) {
      of << " " << s.op;
      if (s.middleAncilla !=
          std::numeric_limits<decltype(s.middleAncilla)>::max()) {
        of << " " << s.middleAncilla;
      }
    }
    of << ",";
  }
  if (!swaps.empty()) {
    of.seekp(-1, std::ios_base::cur); // remove last comma
  }
  of << std::endl;
};

void DataLogger::logMappingResult(MappingResults& result) {
  if (deactivated) {
    return;
  }

  // load output file
  auto of = std::ofstream(dataLoggingPath + "mapping_result.json");
  if (!of.good()) { // if loading failed, output warning and deactivate logging
    deactivated = true;
    std::cerr << "[data-logging] Error opening file: " << dataLoggingPath
              << "mapping_result.json" << std::endl;
    return;
  }

  // prepare json data
  nlohmann::json json      = result.json();
  auto&          stats     = json["statistics"];
  auto&          benchmark = stats["benchmark"];
  for (std::size_t i = 0; i < result.layerHeuristicBenchmark.size(); ++i) {
    auto& layerBenchmark                 = result.layerHeuristicBenchmark.at(i);
    auto& jsonLayerBenchmark             = benchmark["layers"][i];
    jsonLayerBenchmark["expanded_nodes"] = layerBenchmark.expandedNodes;
    jsonLayerBenchmark["generated_nodes"] = layerBenchmark.generatedNodes;
    jsonLayerBenchmark["solution_depth"]  = layerBenchmark.solutionDepth;
    jsonLayerBenchmark["time_per_node"]   = layerBenchmark.timePerNode;
    jsonLayerBenchmark["average_branching_factor"] =
        layerBenchmark.averageBranchingFactor;
    jsonLayerBenchmark["effective_branching_factor"] =
        layerBenchmark.effectiveBranchingFactor;
  }

  // write json data to file
  of << json.dump(2);
  of.close();
};

void DataLogger::close() {
  for (std::size_t i = 0; i < searchNodesLogFiles.size(); ++i) {
    if (searchNodesLogFiles.at(i).is_open()) {
      std::cerr << "[data-logging] Error: layer " << i << " was not finalized"
                << std::endl;
      searchNodesLogFiles.at(i).close();
    }
  }
  deactivated = true;
}
