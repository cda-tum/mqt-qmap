//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "heuristic/HeuristicMapper.hpp"
#include "nlohmann/json.hpp"

#include "gtest/gtest.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stack>
#include <string>

TEST(Functionality, NodeCostCalculation) {
  const double                  tolerance = 1e-6;
  const CouplingMap             cm        = {{0, 1}, {1, 2}, {3, 1}, {4, 3}};
  Architecture                  arch{5, cm};
  const SingleQubitMultiplicity empty1Mult                  = {};
  const std::unordered_set<std::uint16_t>& consideredQubits = {0, 1, 2, 3, 5};
  const TwoQubitMultiplicity               multiplicity     = {{{0, 1}, {5, 2}},
                                                               {{2, 3}, {0, 1}}};
  const std::array<std::int16_t, MAX_DEVICE_QUBITS> qubits  = {4, 3, 1, 2, 0};
  const std::array<std::int16_t, MAX_DEVICE_QUBITS> locations = {4, 2, 3, 1, 0};

  const std::vector<std::vector<Exchange>> swaps = {
      {Exchange(0, 1, qc::OpType::Teleportation)},
      {Exchange(1, 2, qc::OpType::SWAP)}};

  HeuristicMapper::Node node(0, 0, qubits, locations, swaps, 5., 0, false,
                             false);
  node.updateHeuristicCost(arch, empty1Mult, multiplicity, consideredQubits);
  EXPECT_NEAR(node.costHeur,
              COST_UNIDIRECTIONAL_SWAP * 14 + COST_DIRECTION_REVERSE * 3,
              tolerance);

  node =
      HeuristicMapper::Node(0, 0, qubits, locations, swaps, 5., 0, false, true);
  EXPECT_NEAR(node.costFixed, 5., tolerance);
  node.updateHeuristicCost(arch, empty1Mult, multiplicity, consideredQubits);
  EXPECT_NEAR(node.costHeur,
              COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE, tolerance);
  node.applySWAP({3, 4}, arch, empty1Mult, multiplicity);
  node.updateHeuristicCost(arch, empty1Mult, multiplicity, consideredQubits);
  EXPECT_NEAR(node.costFixed, 5. + COST_UNIDIRECTIONAL_SWAP, tolerance);
  EXPECT_NEAR(node.costHeur, COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE,
              tolerance);
  node.lookaheadPenalty = 0.;
  EXPECT_NEAR(node.getTotalCost(),
              5. + COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalFixedCost(), 5. + COST_UNIDIRECTIONAL_SWAP,
              tolerance);
  node.lookaheadPenalty = 2.;
  EXPECT_NEAR(node.getTotalCost(),
              7. + COST_UNIDIRECTIONAL_SWAP * 2 + COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalFixedCost(), 7. + COST_UNIDIRECTIONAL_SWAP,
              tolerance);
  node.recalculateFixedCost(arch, empty1Mult, multiplicity);
  EXPECT_NEAR(node.costFixed, COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 2,
              tolerance);
  EXPECT_NEAR(node.costHeur, COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalCost(),
              2. + COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 3 +
                  COST_DIRECTION_REVERSE,
              tolerance);
  EXPECT_NEAR(node.getTotalFixedCost(),
              2. + COST_TELEPORTATION + COST_UNIDIRECTIONAL_SWAP * 2,
              tolerance);
}

TEST(Functionality, HeuristicBenchmark) {
  /*
      3
     / \
    4   2
    |   |
    0---1
  */
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
                          {3, 2}, {3, 4}, {4, 3}, {4, 0}, {0, 4}};
  architecture.loadCouplingMap(5, cm);

  qc::QuantumComputation qc{5, 5};
  qc.cx(qc::Control{4}, 2);
  qc.cx(qc::Control{3}, 1);
  qc.cx(qc::Control{4}, 1);

  qc.barrier({0, 1, 2, 3, 4});
  for (size_t i = 0; i < 5; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  const auto    mapper = std::make_unique<HeuristicMapper>(qc, architecture);
  Configuration settings{};
  settings.admissibleHeuristic      = true;
  settings.layering                 = Layering::DisjointQubits;
  settings.initialLayout            = InitialLayout::Identity;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  settings.lookahead                = false;
  settings.debug                    = true;
  mapper->map(settings);
  auto& result = mapper->getResults();

  /*
  generated nodes (unit of costs: COST_BIDIRECTIONAL_SWAP):
  layer 1:
    0: {swaps: {}, cost: 0, heur: 1, total: 1}
  --- priority queue: [0] -> expand: 0
    1: {swaps: {{0, 1}}, cost: 1, heur: 1, total: 2}
    2: {swaps: {{1, 2}}, cost: 1, heur: 1, total: 2}
    3: {swaps: {{2, 3}}, cost: 1, heur: 0, total: 1}
    4: {swaps: {{3, 4}}, cost: 1, heur: 1, total: 2}
    5: {swaps: {{4, 0}}, cost: 1, heur: 1, total: 2}
  --- priority queue: [3,1,2,4,5] -> done: 3
  layer 2:
    0: {swaps: {}, cost: 0, heur: 1, total: 1}
  --- priority queue: [0] -> expand: 0
    1: {swaps: {{0, 1}}, cost: 1, heur: 0, total: 1}
    2: {swaps: {{1, 2}}, cost: 1, heur: 1, total: 2}
    3: {swaps: {{3, 4}}, cost: 1, heur: 1, total: 2}
    4: {swaps: {{4, 0}}, cost: 1, heur: 0, total: 1}
  --- priority queue: [1,4,2,3] -> done: 1
  */

  EXPECT_EQ(result.layerHeuristicBenchmark.size(), 2);
  const auto& layerResults0 = result.layerHeuristicBenchmark[0];
  EXPECT_EQ(layerResults0.solutionDepth, 1);
  EXPECT_EQ(layerResults0.generatedNodes, 6);
  EXPECT_EQ(layerResults0.expandedNodes, 1);
  EXPECT_NEAR(layerResults0.averageBranchingFactor, 5.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_NEAR(layerResults0.effectiveBranchingFactor, 1.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  const auto& layerResults1 = result.layerHeuristicBenchmark[1];
  EXPECT_EQ(layerResults1.solutionDepth, 1);
  EXPECT_EQ(layerResults1.generatedNodes, 5);
  EXPECT_EQ(layerResults1.expandedNodes, 1);
  EXPECT_NEAR(layerResults1.averageBranchingFactor, 4.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_NEAR(layerResults1.effectiveBranchingFactor, 1.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_EQ(result.heuristicBenchmark.generatedNodes, 11);
  EXPECT_EQ(result.heuristicBenchmark.expandedNodes, 2);
  EXPECT_NEAR(result.heuristicBenchmark.averageBranchingFactor, 4.5,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
  EXPECT_NEAR(result.heuristicBenchmark.effectiveBranchingFactor, 1.,
              HeuristicMapper::EFFECTIVE_BRANCH_RATE_TOLERANCE);
}

TEST(Functionality, EmptyDump) {
  qc::QuantumComputation qc{1};
  qc.x(0);
  Architecture    arch{1, {}};
  HeuristicMapper mapper(qc, arch);
  mapper.dumpResult("test.qasm");
  mapper.map({});
  EXPECT_NO_THROW(mapper.dumpResult("test.qasm"););
  EXPECT_NO_THROW(mapper.dumpResult("test.real"););
  EXPECT_THROW(mapper.dumpResult("test.dummy"), QMAPException);
}

TEST(Functionality, InvalidSettings) {
  qc::QuantumComputation qc{1};
  qc.x(0);
  Architecture arch{1, {}};
  auto         props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.1);
  arch.loadProperties(props);
  HeuristicMapper mapper(qc, arch);
  auto            config = Configuration{};
  config.method          = Method::Heuristic;
  config.layering        = Layering::OddGates;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.layering = Layering::QubitTriangle;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.layering         = Layering::IndividualGates;
  config.considerFidelity = true;
  config.lookahead        = true;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.lookahead           = false;
  config.teleportationQubits = 2;
  EXPECT_THROW(mapper.map(config), QMAPException);
  config.teleportationQubits = 0;
  EXPECT_NO_THROW(mapper.map(config));
}

TEST(Functionality, NoMeasurmentsAdded) {
  using namespace qc::literals;
  // construct circuit
  qc::QuantumComputation qc{4U};
  qc.cx(0_pc, 1);
  qc.cx(2_pc, 1);
  qc.cx(3_pc, 1);

  // load architecture
  Architecture arch{};
  arch.loadCouplingMap(AvailableArchitecture::IbmqLondon);

  // create heuristic mapper
  HeuristicMapper mapper(qc, arch);

  // configure to not include measurements after mapping
  auto config                           = Configuration{};
  config.addMeasurementsToMappedCircuit = false;

  // perform the mapping
  mapper.map(config);

  // get the resulting circuit
  auto              qcMapped = qc::QuantumComputation();
  std::stringstream qasm{};
  mapper.dumpResult(qasm, qc::Format::OpenQASM);
  qcMapped.import(qasm, qc::Format::OpenQASM);

  // check no measurements were added
  EXPECT_EQ(qcMapped.getNops(), 3U);
  EXPECT_NE(qcMapped.back()->getType(), qc::Measure);
}

TEST(Functionality, InvalidCircuits) {
  Configuration config{};
  config.method = Method::Heuristic;

  // architecture not connected
  qc::QuantumComputation qc{2U};
  qc.cx({0}, 1);
  Architecture    arch{2U, {}};
  HeuristicMapper mapper(qc, arch);
  EXPECT_THROW(mapper.map(config), QMAPException);

  // gate with >1 control
  qc::QuantumComputation      qc3{3U};
  const qc::StandardOperation op =
      qc::StandardOperation(3, {{0}, {1}}, qc::Qubit{2});
  qc3.emplace_back(op.clone());
  Architecture    arch2{3U, {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 0}, {0, 2}}};
  HeuristicMapper mapper3(qc3, arch2);
  EXPECT_THROW(mapper3.map(config), QMAPException);
}

TEST(Functionality, HeuristicAdmissibility) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
                          {3, 2}, {3, 4}, {4, 3}, {4, 5}, {5, 4}};
  architecture.loadCouplingMap(6, cm);
  const std::vector<Edge> perms{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}};

  const SingleQubitMultiplicity            empty1Mult       = {};
  const std::unordered_set<std::uint16_t>& consideredQubits = {0, 1, 2,
                                                               3, 4, 5};
  const TwoQubitMultiplicity               multiplicity     = {
      {{0, 4}, {1, 0}}, {{1, 3}, {1, 0}}, {{2, 5}, {1, 0}}};

  // perform depth-limited depth first search
  const std::size_t                  depthLimit = 6;
  std::vector<HeuristicMapper::Node> nodeStack{};
  nodeStack.reserve(depthLimit);
  std::stack<std::size_t> permStack{};

  auto initNode = HeuristicMapper::Node(
      0, 0, {0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5}, {}, 0., 0, false, true);
  initNode.recalculateFixedCost(architecture, empty1Mult, multiplicity);
  initNode.updateHeuristicCost(architecture, empty1Mult, multiplicity,
                               consideredQubits);
  nodeStack.emplace_back(initNode);
  permStack.emplace(perms.size());

  while (!nodeStack.empty()) {
    const auto& node = nodeStack.back();
    if (node.done) {
      // check if all nodes in stack have lower or equal cost
      for (const auto& prevNode : nodeStack) {
        EXPECT_LE(prevNode.getTotalCost(), node.getTotalCost());
      }
    }
    if (node.done || nodeStack.size() >= depthLimit || permStack.top() == 0) {
      nodeStack.pop_back();
      permStack.pop();
      continue;
    }
    --permStack.top();
    const auto perm = perms[permStack.top()];
    auto       newNode =
        HeuristicMapper::Node(1, 0, node.qubits, node.locations, node.swaps,
                              node.costFixed, node.depth + 1, false, true);
    newNode.applySWAP(perm, architecture, empty1Mult, multiplicity);
    newNode.updateHeuristicCost(architecture, empty1Mult, multiplicity,
                                consideredQubits);
    nodeStack.emplace_back(newNode);
    permStack.emplace(perms.size());
  }
}

TEST(Functionality, DataLoggerAfterClose) {
  const std::string      dataLoggingPath = "test_log/";
  qc::QuantumComputation qc{3};
  qc.x(0);
  Architecture                arch{3, {}};
  std::unique_ptr<DataLogger> dataLogger =
      std::make_unique<DataLogger>(dataLoggingPath, arch, qc);
  const qc::CompoundOperation compOp(3);
  Exchange                    teleport(0, 2, 1, qc::OpType::Teleportation);
  dataLogger->logSearchNode(0, 0, 0, 0., 0., 0., {}, false, {{teleport}}, 0);
  dataLogger->logSearchNode(1, 0, 0, 0., 0., 0., {}, false, {}, 0);
  dataLogger->splitLayer();
  dataLogger->logFinalizeLayer(0, compOp, {}, {}, {}, 0, 0., 0., 0., {}, {}, 0);
  dataLogger->logFinalizeLayer(0, compOp, {}, {}, {}, 0, 0., 0., 0., {}, {}, 0);
  dataLogger->logSearchNode(0, 0, 0, 0., 0., 0., {}, false, {}, 0);
  dataLogger->close();
  dataLogger->clearLog();

  dataLogger->logArchitecture();
  dataLogger->logInputCircuit(qc);
  dataLogger->logOutputCircuit(qc);
  dataLogger->logSearchNode(0, 0, 0, 0., 0., 0., {}, false, {}, 0);
  dataLogger->logFinalizeLayer(0, compOp, {}, {}, {}, 0, 0., 0., 0., {}, {}, 0);
  dataLogger->splitLayer();
  MappingResults result;
  dataLogger->logMappingResult(result);

  // count files and subdirectories in data logging path
  std::size_t fileCount = 0;
  for ([[maybe_unused]] const auto& _ :
       std::filesystem::directory_iterator(dataLoggingPath)) {
    ++fileCount;
  }
  EXPECT_EQ(fileCount, 0);
}

TEST(Functionality, DataLogger) {
  // setting up example architecture and circuit
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {1, 3}, {3, 1}};
  architecture.loadCouplingMap(4, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.1);
  props.setSingleQubitErrorRate(1, "x", 0.2);
  props.setSingleQubitErrorRate(2, "x", 0.3);
  props.setSingleQubitErrorRate(3, "x", 0.4);

  props.setTwoQubitErrorRate(0, 1, 0.5);
  props.setTwoQubitErrorRate(1, 0, 0.5);
  props.setTwoQubitErrorRate(1, 2, 0.6);
  props.setTwoQubitErrorRate(2, 1, 0.6);
  props.setTwoQubitErrorRate(1, 3, 0.7);
  props.setTwoQubitErrorRate(3, 1, 0.7);

  architecture.loadProperties(props);
  architecture.setName("test_architecture");

  qc::QuantumComputation qc{4, 4};
  qc.cx(1, 0);
  qc.cx(3, 2);
  qc.setName("test_circ");

  Configuration settings{};
  settings.admissibleHeuristic      = true;
  settings.layering                 = Layering::IndividualGates;
  settings.initialLayout            = InitialLayout::Identity;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  settings.lookahead                = true;
  settings.nrLookaheads             = 1;
  settings.firstLookaheadFactor     = 0.5;
  settings.lookaheadFactor          = 0.9;
  settings.debug                    = true;
  settings.considerFidelity         = false;
  settings.useTeleportation         = false;
  // setting data logging path to enable data logging
  settings.dataLoggingPath = "test_log";

  // remove directory at data logging path if it already exists
  if (std::filesystem::exists(settings.dataLoggingPath)) {
    std::filesystem::remove_all(settings.dataLoggingPath);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);
  mapper->map(settings);
  mapper->printResult(std::cout);
  MappingResults& results = mapper->getResults();

  // comparing logged architecture information with original architecture object
  auto archFile =
      std::ifstream(settings.dataLoggingPath + "/architecture.json");
  if (!archFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/architecture.json";
  }
  const auto archJson = nlohmann::json::parse(archFile);
  EXPECT_EQ(archJson["name"], architecture.getName());
  EXPECT_EQ(archJson["nqubits"], architecture.getNqubits());
  EXPECT_EQ(archJson["distances"], architecture.getDistanceTable());
  EXPECT_EQ(archJson["coupling_map"], architecture.getCouplingMap());
  const auto& fidelityJson = archJson["fidelity"];
  EXPECT_EQ(fidelityJson["fidelity_distances"],
            architecture.getFidelityDistanceTables());
  EXPECT_EQ(fidelityJson["single_qubit_fidelities"],
            architecture.getSingleQubitFidelities());
  EXPECT_EQ(fidelityJson["two_qubit_fidelities"],
            architecture.getFidelityTable());
  // json does not support inf values, instead nlohmann::json replaces inf
  // with null
  const auto& singleQubitFidelityCosts =
      architecture.getSingleQubitFidelityCosts();
  for (std::size_t i = 0; i < singleQubitFidelityCosts.size(); ++i) {
    if (std::isinf(singleQubitFidelityCosts[i])) {
      EXPECT_TRUE(fidelityJson["single_qubit_fidelity_costs"][i].is_null());
    } else {
      EXPECT_EQ(fidelityJson["single_qubit_fidelity_costs"][i],
                singleQubitFidelityCosts[i]);
    }
  }
  const auto& twoQubitFidelityCosts = architecture.getTwoQubitFidelityCosts();
  for (std::size_t i = 0; i < twoQubitFidelityCosts.size(); ++i) {
    for (std::size_t j = 0; j < twoQubitFidelityCosts[i].size(); ++j) {
      if (std::isinf(twoQubitFidelityCosts[i][j])) {
        EXPECT_TRUE(fidelityJson["two_qubit_fidelity_costs"][i][j].is_null());
      } else {
        EXPECT_EQ(fidelityJson["two_qubit_fidelity_costs"][i][j],
                  twoQubitFidelityCosts[i][j]);
      }
    }
  }
  const auto& swapFidelityCosts = architecture.getSwapFidelityCosts();
  for (std::size_t i = 0; i < swapFidelityCosts.size(); ++i) {
    for (std::size_t j = 0; j < swapFidelityCosts[i].size(); ++j) {
      if (std::isinf(swapFidelityCosts[i][j])) {
        EXPECT_TRUE(fidelityJson["swap_fidelity_costs"][i][j].is_null());
      } else {
        EXPECT_EQ(fidelityJson["swap_fidelity_costs"][i][j],
                  swapFidelityCosts[i][j]);
      }
    }
  }

  // comparing logged mapping result with mapping result object
  auto resultFile =
      std::ifstream(settings.dataLoggingPath + "/mapping_result.json");
  if (!resultFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/mapping_result.json";
  }
  const auto  resultJson = nlohmann::json::parse(resultFile);
  const auto& configJson = resultJson["config"];
  EXPECT_EQ(configJson["add_measurements_to_mapped_circuit"],
            settings.addMeasurementsToMappedCircuit);
  EXPECT_EQ(configJson["debug"], settings.debug);
  EXPECT_EQ(configJson["verbose"], settings.verbose);
  EXPECT_EQ(configJson["layering_strategy"], toString(settings.layering));
  EXPECT_EQ(configJson["method"], toString(settings.method));
  EXPECT_EQ(configJson["post_mapping_optimizations"],
            settings.postMappingOptimizations);
  EXPECT_EQ(configJson["pre_mapping_optimizations"],
            settings.preMappingOptimizations);
  const auto& heuristicSettingsJson = configJson["settings"];
  EXPECT_EQ(heuristicSettingsJson["admissible_heuristic"],
            settings.admissibleHeuristic);
  EXPECT_EQ(heuristicSettingsJson["consider_fidelity"],
            settings.considerFidelity);
  EXPECT_EQ(heuristicSettingsJson["initial_layout"],
            toString(settings.initialLayout));
  if (settings.lookahead) {
    const auto& lookaheadJson = heuristicSettingsJson["lookahead"];
    EXPECT_EQ(lookaheadJson["factor"], settings.lookaheadFactor);
    EXPECT_EQ(lookaheadJson["first_factor"], settings.firstLookaheadFactor);
    EXPECT_EQ(lookaheadJson["lookaheads"], settings.nrLookaheads);
  }
  if (settings.useTeleportation) {
    const auto& teleportationJson = heuristicSettingsJson["teleportation"];
    EXPECT_EQ(teleportationJson["qubits"], settings.teleportationQubits);
    EXPECT_EQ(teleportationJson["seed"], settings.teleportationSeed);
    EXPECT_EQ(teleportationJson["fake"], settings.teleportationFake);
  }

  const auto& inCircJson = resultJson["circuit"];
  EXPECT_EQ(inCircJson["cnots"], results.input.cnots);
  EXPECT_EQ(inCircJson["gates"], results.input.gates);
  EXPECT_EQ(inCircJson["name"], results.input.name);
  EXPECT_EQ(inCircJson["qubits"], results.input.qubits);
  EXPECT_EQ(inCircJson["single_qubit_gates"], results.input.singleQubitGates);

  const auto& outCircJson = resultJson["mapped_circuit"];
  EXPECT_EQ(outCircJson["cnots"], results.output.cnots);
  EXPECT_EQ(outCircJson["gates"], results.output.gates);
  EXPECT_EQ(outCircJson["name"], results.output.name);
  EXPECT_EQ(outCircJson["qubits"], results.output.qubits);
  EXPECT_EQ(outCircJson["single_qubit_gates"], results.output.singleQubitGates);

  const auto& statJson = resultJson["statistics"];
  EXPECT_EQ(statJson["additional_gates"],
            (static_cast<std::make_signed_t<decltype(results.output.gates)>>(
                 results.output.gates) -
             static_cast<std::make_signed_t<decltype(results.input.gates)>>(
                 results.input.gates)));
  EXPECT_EQ(statJson["layers"], results.input.layers);
  EXPECT_EQ(statJson["mapping_time"], results.time);
  EXPECT_EQ(statJson["swaps"], results.output.swaps);
  EXPECT_EQ(statJson["teleportations"], results.output.teleportations);
  EXPECT_EQ(statJson["total_fidelity"], results.output.totalFidelity);
  EXPECT_EQ(statJson["total_log_fidelity"], results.output.totalLogFidelity);

  const auto& benchmarkJson = statJson["benchmark"];
  EXPECT_EQ(benchmarkJson["average_branching_factor"],
            results.heuristicBenchmark.averageBranchingFactor);
  EXPECT_EQ(benchmarkJson["effective_branching_factor"],
            results.heuristicBenchmark.effectiveBranchingFactor);
  EXPECT_EQ(benchmarkJson["expanded_nodes"],
            results.heuristicBenchmark.expandedNodes);
  EXPECT_EQ(benchmarkJson["generated_nodes"],
            results.heuristicBenchmark.generatedNodes);
  EXPECT_EQ(benchmarkJson["time_per_node"],
            results.heuristicBenchmark.timePerNode);
  const auto& benchmarkLayersJson = benchmarkJson["layers"];
  EXPECT_EQ(benchmarkLayersJson.size(), results.layerHeuristicBenchmark.size());
  for (std::size_t i = 0; i < results.layerHeuristicBenchmark.size(); ++i) {
    EXPECT_EQ(benchmarkLayersJson[i]["average_branching_factor"],
              results.layerHeuristicBenchmark.at(i).averageBranchingFactor);
    EXPECT_EQ(benchmarkLayersJson[i]["effective_branching_factor"],
              results.layerHeuristicBenchmark.at(i).effectiveBranchingFactor);
    EXPECT_EQ(benchmarkLayersJson[i]["expanded_nodes"],
              results.layerHeuristicBenchmark.at(i).expandedNodes);
    EXPECT_EQ(benchmarkLayersJson[i]["generated_nodes"],
              results.layerHeuristicBenchmark.at(i).generatedNodes);
    EXPECT_EQ(benchmarkLayersJson[i]["solution_depth"],
              results.layerHeuristicBenchmark.at(i).solutionDepth);
    EXPECT_EQ(benchmarkLayersJson[i]["time_per_node"],
              results.layerHeuristicBenchmark.at(i).timePerNode);
  }

  // comparing logged input and output circuits with input circuit object and
  // mapped circuit object
  auto inputQasmFile = std::ifstream(settings.dataLoggingPath + "/input.qasm");
  if (!inputQasmFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/input.qasm";
  }
  std::stringstream inputFileBuffer;
  inputFileBuffer << inputQasmFile.rdbuf();
  std::stringstream inputQasmBuffer;
  qc.dumpOpenQASM(inputQasmBuffer);
  EXPECT_EQ(inputFileBuffer.str(), inputQasmBuffer.str());

  auto outputQasmFile =
      std::ifstream(settings.dataLoggingPath + "/output.qasm");
  if (!outputQasmFile.is_open()) {
    FAIL() << "Could not open file " << settings.dataLoggingPath
           << "/output.qasm";
  }
  std::stringstream outputFileBuffer;
  outputFileBuffer << outputQasmFile.rdbuf();
  std::stringstream outputQasmBuffer;
  mapper->dumpResult(outputQasmBuffer, qc::Format::OpenQASM);
  EXPECT_EQ(outputFileBuffer.str(), outputQasmBuffer.str());

  // checking logged search graph info against known values (correct qubit
  // number, valid layouts, correct data types in all csv fields, etc.)
  for (std::size_t i = 0; i < results.input.layers; ++i) {
    auto layerFile = std::ifstream(settings.dataLoggingPath + "/layer_" +
                                   std::to_string(i) + ".json");
    if (!layerFile.is_open()) {
      FAIL() << "Could not open file " << settings.dataLoggingPath << "/layer_"
             << i << ".json";
    }
    const auto        layerJson   = nlohmann::json::parse(layerFile);
    const std::size_t finalNodeId = layerJson["final_node_id"];
    EXPECT_EQ(layerJson["initial_layout"].size(), architecture.getNqubits());
    EXPECT_EQ(layerJson["single_qubit_multiplicity"].size(),
              architecture.getNqubits());

    auto layerNodeFile =
        std::ifstream(settings.dataLoggingPath + "/nodes_layer_" +
                      std::to_string(i) + ".csv");
    if (!layerNodeFile.is_open()) {
      FAIL() << "Could not open file " << settings.dataLoggingPath
             << "/nodes_layer_" << i << ".csv";
    }
    std::string           line;
    bool                  foundFinalNode = false;
    std::set<std::size_t> nodeIds;
    while (std::getline(layerNodeFile, line)) {
      if (line.empty()) {
        continue;
      }
      std::string                                        col;
      std::size_t                                        nodeId           = 0;
      std::size_t                                        parentId         = 0;
      std::size_t                                        depth            = 0;
      std::size_t                                        isValidMapping   = 0;
      double                                             costFixed        = 0.;
      double                                             costHeur         = 0.;
      double                                             lookaheadPenalty = 0.;
      std::vector<std::int32_t>                          layout{};
      std::vector<std::pair<std::int16_t, std::int16_t>> swaps{};
      std::stringstream                                  lineStream(line);
      if (std::getline(lineStream, col, ';')) {
        nodeId = std::stoull(col);
        nodeIds.insert(nodeId);
      } else {
        FAIL() << "Missing value for node id in " << settings.dataLoggingPath
               << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        parentId = std::stoull(col);
        if (nodeId != 0) {
          EXPECT_TRUE(nodeIds.count(parentId) > 0);
        }
      } else {
        FAIL() << "Missing value for parent node id in "
               << settings.dataLoggingPath << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        costFixed = std::stod(col);
      } else {
        FAIL() << "Missing value for fixed cost in " << settings.dataLoggingPath
               << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        costHeur = std::stod(col);
      } else {
        FAIL() << "Missing value for heuristic cost in "
               << settings.dataLoggingPath << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        lookaheadPenalty = std::stod(col);
      } else {
        FAIL() << "Missing value for lookahead penalty in "
               << settings.dataLoggingPath << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        isValidMapping = std::stoull(col);
        if (isValidMapping > 1) {
          FAIL() << "Non-boolean value " << isValidMapping
                 << " for isValidMapping in " << settings.dataLoggingPath
                 << "/nodes_layer_" << i << ".csv";
        }
      } else {
        FAIL() << "Missing value for isValidMapping in "
               << settings.dataLoggingPath << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        depth = std::stoull(col);
      } else {
        FAIL() << "Missing value for depth in " << settings.dataLoggingPath
               << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        std::stringstream qubitMapBuffer(col);
        std::string       entry;
        while (std::getline(qubitMapBuffer, entry, ',')) {
          const std::int32_t qubit = std::stoi(entry);
          layout.push_back(qubit);
          EXPECT_TRUE(-1 <= qubit && qubit < architecture.getNqubits());
        }
        EXPECT_EQ(layout.size(), architecture.getNqubits());
      } else {
        FAIL() << "Missing value for layout in " << settings.dataLoggingPath
               << "/nodes_layer_" << i << ".csv";
      }
      if (std::getline(lineStream, col, ';')) {
        std::stringstream swapBuffer(col);
        std::string       entry;
        while (std::getline(swapBuffer, entry, ',')) {
          std::int32_t q1 = 0;
          std::int32_t q2 = 0;
          std::stringstream(entry) >> q1 >> q2;
          EXPECT_TRUE(0 <= q1 && q1 < architecture.getNqubits());
          EXPECT_TRUE(0 <= q2 && q2 < architecture.getNqubits());
          swaps.emplace_back(q1, q2);
        }
      }

      if (nodeId == finalNodeId) {
        foundFinalNode = true;
        EXPECT_EQ(layerJson["final_cost_fixed"], costFixed);
        EXPECT_EQ(layerJson["final_cost_heur"], costHeur);
        EXPECT_EQ(layerJson["final_layout"], layout);
        EXPECT_EQ(layerJson["final_lookahead_penalty"], lookaheadPenalty);
        EXPECT_EQ(layerJson["final_search_depth"], depth);
        EXPECT_EQ(layerJson["final_swaps"], swaps);
        EXPECT_EQ(isValidMapping, 1);
      }
    }
    if (!foundFinalNode) {
      FAIL() << "Could not find final node in " << settings.dataLoggingPath
             << "/nodes_layer_" << i << ".csv";
    }
  }

  // checking if files for non-existing layers are not created
  auto afterLastLayerFile =
      std::ifstream(settings.dataLoggingPath + "/layer_" +
                    std::to_string(results.input.layers) + ".json");
  if (afterLastLayerFile.is_open()) {
    FAIL() << "File " << settings.dataLoggingPath << "/layer_"
           << results.input.layers
           << ".json should not exist, as there are not that many layers";
  }
  auto afterLastLayerNodesFile =
      std::ifstream(settings.dataLoggingPath + "/nodes_layer_" +
                    std::to_string(results.input.layers) + ".csv");
  if (afterLastLayerNodesFile.is_open()) {
    FAIL() << "File " << settings.dataLoggingPath << "/nodes_layer_"
           << results.input.layers
           << ".csv should not exist, as there are not that many layers";
  }
}

class HeuristicTest5Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";
  std::string testCalibrationDir  = "../extern/calibration/";

  qc::QuantumComputation           qc{};
  Architecture                     ibmqYorktown{};
  Architecture                     ibmqLondon{};
  std::unique_ptr<HeuristicMapper> ibmqYorktownMapper;
  std::unique_ptr<HeuristicMapper> ibmqLondonMapper;
  Configuration                    settings{};

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    ibmqYorktown.loadCouplingMap(AvailableArchitecture::IbmqYorktown);
    ibmqLondon.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    ibmqLondon.loadProperties(testCalibrationDir + "ibmq_london.csv");
    ibmqYorktownMapper = std::make_unique<HeuristicMapper>(qc, ibmqYorktown);
    ibmqLondonMapper   = std::make_unique<HeuristicMapper>(qc, ibmqLondon);
    settings.debug     = true;
    settings.verbose   = true;

    settings.iterativeBidirectionalRouting       = true;
    settings.iterativeBidirectionalRoutingPasses = 3;
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest5Q,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84",
                    "4mod5-v0_20", "mod5d1_63"),
    [](const testing::TestParamInfo<HeuristicTest5Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest5Q, Identity) {
  settings.initialLayout = InitialLayout::Identity;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_heuristic_qx4_identity.qasm");
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_heuristic_london_identity.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Static) {
  settings.initialLayout = InitialLayout::Static;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_heuristic_qx4_static.qasm");
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_heuristic_london_static.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Dynamic) {
  settings.initialLayout = InitialLayout::Dynamic;
  ibmqYorktownMapper->map(settings);
  ibmqYorktownMapper->dumpResult(GetParam() + "_heuristic_qx4_dynamic.qasm");
  ibmqYorktownMapper->printResult(std::cout);

  ibmqLondonMapper->map(settings);
  ibmqLondonMapper->dumpResult(GetParam() + "_heuristic_london_dynamic.qasm");
  ibmqLondonMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest16Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     ibmQX5{};
  std::unique_ptr<HeuristicMapper> ibmQX5Mapper;
  Configuration                    settings{};

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    ibmQX5.loadCouplingMap(AvailableArchitecture::IbmQx5);
    ibmQX5Mapper   = std::make_unique<HeuristicMapper>(qc, ibmQX5);
    settings.debug = true;
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest16Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16"),
    [](const testing::TestParamInfo<HeuristicTest16Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest16Q, Dynamic) {
  settings.initialLayout = InitialLayout::Dynamic;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->dumpResult(GetParam() + "_heuristic_qx5_dynamic.qasm");
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest16Q, Disjoint) {
  settings.layering = Layering::DisjointQubits;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->dumpResult(GetParam() + "_heuristic_qx5_disjoint.qasm");
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest16Q, Disjoint2qBlocks) {
  settings.layering = Layering::Disjoint2qBlocks;
  ibmQX5Mapper->map(settings);
  ibmQX5Mapper->dumpResult(GetParam() + "_heuristic_qx5_disjoint_2q.qasm");
  ibmQX5Mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20Q : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  std::unique_ptr<HeuristicMapper> tokyoMapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IbmqTokyo);
    tokyoMapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTest20Q,
    testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16",
                    "z4_268"),
    [](const testing::TestParamInfo<HeuristicTest20Q::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTest20Q, Dynamic) {
  Configuration settings{};
  settings.initialLayout = InitialLayout::Dynamic;
  settings.debug         = true;
  tokyoMapper->map(settings);
  tokyoMapper->dumpResult(GetParam() + "_heuristic_tokyo_dynamic.qasm");
  tokyoMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTest20QTeleport
    : public testing::TestWithParam<std::tuple<std::uint64_t, std::string>> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  std::unique_ptr<HeuristicMapper> tokyoMapper;

  void SetUp() override {
    qc.import(testExampleDir + std::get<1>(GetParam()) + ".qasm");
    arch.loadCouplingMap(AvailableArchitecture::IbmqTokyo);
    tokyoMapper = std::make_unique<HeuristicMapper>(qc, arch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    HeuristicTeleport, HeuristicTest20QTeleport,
    testing::Combine(testing::Values(0, 1, 2, 3, 1337, 1338, 3147),
                     testing::Values("ising_model_10", "rd73_140", "cnt3-5_179",
                                     "qft_16", "z4_268")),
    [](const testing::TestParamInfo<HeuristicTest20QTeleport::ParamType>& inf) {
      std::string name = std::get<1>(inf.param);
      std::replace(name.begin(), name.end(), '-', '_');
      std::stringstream ss{};
      ss << name << "_seed" << std::get<0>(inf.param);
      return ss.str();
    });

TEST_P(HeuristicTest20QTeleport, Teleportation) {
  Configuration settings{};
  settings.initialLayout       = InitialLayout::Dynamic;
  settings.debug               = true;
  settings.verbose             = true;
  settings.teleportationQubits = std::min(
      (arch.getNqubits() - qc.getNqubits()) & ~1U, static_cast<std::size_t>(8));
  settings.teleportationSeed = std::get<0>(GetParam());
  tokyoMapper->map(settings);
  tokyoMapper->dumpResult(std::get<1>(GetParam()) +
                          "_heuristic_tokyo_teleport.qasm");
  tokyoMapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

class HeuristicTestFidelity : public testing::TestWithParam<std::string> {
protected:
  std::string testExampleDir      = "../examples/";
  std::string testArchitectureDir = "../extern/architectures/";
  std::string testCalibrationDir  = "../extern/calibration/";

  qc::QuantumComputation           qc{};
  Architecture                     arch{};
  Architecture                     nonFidelityArch{};
  std::unique_ptr<HeuristicMapper> mapper;
  std::unique_ptr<HeuristicMapper> nonFidelityMapper;

  void SetUp() override {
    qc.import(testExampleDir + GetParam() + ".qasm");
    arch.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
    arch.loadProperties(testCalibrationDir + "ibmq_london.csv");
    mapper = std::make_unique<HeuristicMapper>(qc, arch);
    nonFidelityArch.loadCouplingMap(AvailableArchitecture::IbmqYorktown);
    nonFidelityMapper = std::make_unique<HeuristicMapper>(qc, nonFidelityArch);
  }
};

INSTANTIATE_TEST_SUITE_P(
    Heuristic, HeuristicTestFidelity,
    testing::Values("3_17_13", "ex-1_166", "ham3_102", "miller_11", "4gt11_84",
                    "4mod5-v0_20", "mod5d1_63"),
    [](const testing::TestParamInfo<HeuristicTestFidelity::ParamType>& inf) {
      std::string name = inf.param;
      std::replace(name.begin(), name.end(), '-', '_');
      return name;
    });

TEST_P(HeuristicTestFidelity, Identity) {
  Configuration settings{};
  settings.layering         = Layering::DisjointQubits;
  settings.initialLayout    = InitialLayout::Identity;
  settings.considerFidelity = true;
  settings.lookahead        = false;
  mapper->map(settings);
  mapper->dumpResult(GetParam() + "_heuristic_london_fidelity_identity.qasm");
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, Static) {
  Configuration settings{};
  settings.layering         = Layering::DisjointQubits;
  settings.initialLayout    = InitialLayout::Static;
  settings.considerFidelity = true;
  settings.lookahead        = false;
  mapper->map(settings);
  mapper->dumpResult(GetParam() + "_heuristic_london_fidelity_static.qasm");
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, Dynamic) {
  Configuration settings{};
  settings.layering         = Layering::DisjointQubits;
  settings.initialLayout    = InitialLayout::Dynamic;
  settings.considerFidelity = true;
  settings.lookahead        = false;
  mapper->map(settings);
  mapper->dumpResult(GetParam() + "_heuristic_london_fidelity_static.qasm");
  mapper->printResult(std::cout);
  SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTestFidelity, NoFidelity) {
  Configuration settings{};
  settings.layering         = Layering::DisjointQubits;
  settings.initialLayout    = InitialLayout::Static;
  settings.considerFidelity = true;
  settings.lookahead        = false;
  EXPECT_THROW(nonFidelityMapper->map(settings), QMAPException);
}

TEST(HeuristicTestFidelity, RemapSingleQubit) {
  Architecture      architecture{};
  const CouplingMap cm = {
      {0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3},
      {3, 2}, {3, 4}, {4, 3}, {4, 5}, {5, 4},
  };
  architecture.loadCouplingMap(6, cm);

  const double e5 = 0.99;
  const double e4 = 0.9;
  const double e3 = 0.5;
  const double e1 = 0.1;
  const double e0 = 0.01;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e4);
  props.setSingleQubitErrorRate(4, "x", e4);
  props.setSingleQubitErrorRate(5, "x", e1);

  props.setTwoQubitErrorRate(0, 1, e1);
  props.setTwoQubitErrorRate(1, 0, e1);
  props.setTwoQubitErrorRate(1, 2, e3);
  props.setTwoQubitErrorRate(2, 1, e3);
  props.setTwoQubitErrorRate(2, 3, e5);
  props.setTwoQubitErrorRate(3, 2, e5);
  props.setTwoQubitErrorRate(3, 4, e0);
  props.setTwoQubitErrorRate(4, 3, e0);
  props.setTwoQubitErrorRate(4, 5, e0);
  props.setTwoQubitErrorRate(5, 4, e0);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{6, 6};
  for (std::size_t i = 0; i < 5; ++i) {
    qc.cx(2, 0);
    qc.x(3);
  }

  for (size_t i = 0; i < 6; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.admissibleHeuristic      = true;
  settings.layering                 = Layering::Disjoint2qBlocks;
  settings.initialLayout            = InitialLayout::Identity;
  settings.considerFidelity         = true;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  settings.verbose                  = true;
  settings.swapOnFirstLayer         = true;
  settings.lookahead                = false;
  mapper->map(settings);
  mapper->dumpResult("remap_single_qubit_mapped.qasm");
  mapper->printResult(std::cout);

  // 0 --e1-- 1 --e3-- 2 --e5-- 3 --e0-- 4 --e0-- 5
  // e5       e5       e5       e4       e4       e1

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  /*
  expected output (order of gates may vary):
  SWAP(1,2)
  SWAP(3,4)
  SWAP(4,5)
  CX(0,1) [x5]
  X(5) [x5]
  */
  EXPECT_EQ(result.output.swaps, 3);

  const double c3 = -std::log2(1 - e3);
  const double c1 = -std::log2(1 - e1);
  const double c0 = -std::log2(1 - e0);

  const double expectedFidelity = 3 * c3 + 3 * c0 + 3 * c0 + // SWAPs
                                  5 * c1 +                   // Xs
                                  5 * c1;                    // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, QubitRideAlong) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2},
                          {1, 4}, {4, 1}, {2, 5}, {5, 2}, {5, 6}, {6, 5}};
  architecture.loadCouplingMap(7, cm);

  const double e5 = 0.99;
  const double e4 = 0.9;
  const double e3 = 0.5;
  const double e1 = 0.1;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e5);
  props.setSingleQubitErrorRate(4, "x", e5);
  props.setSingleQubitErrorRate(5, "x", e5);
  props.setSingleQubitErrorRate(6, "x", e5);

  props.setTwoQubitErrorRate(0, 1, e4);
  props.setTwoQubitErrorRate(1, 0, e4);
  props.setTwoQubitErrorRate(1, 2, e3);
  props.setTwoQubitErrorRate(2, 1, e3);
  props.setTwoQubitErrorRate(2, 3, e1);
  props.setTwoQubitErrorRate(3, 2, e1);
  props.setTwoQubitErrorRate(1, 4, e1);
  props.setTwoQubitErrorRate(4, 1, e1);
  props.setTwoQubitErrorRate(2, 5, e3);
  props.setTwoQubitErrorRate(5, 2, e3);
  props.setTwoQubitErrorRate(5, 6, e4);
  props.setTwoQubitErrorRate(6, 5, e4);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{7, 7};
  for (std::size_t i = 0; i < 5; ++i) {
    qc.cx(3, 0);
    qc.cx(6, 4);
  }

  for (size_t i = 0; i < 7; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.admissibleHeuristic      = true;
  settings.layering                 = Layering::Disjoint2qBlocks;
  settings.initialLayout            = InitialLayout::Identity;
  settings.considerFidelity         = true;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  settings.verbose                  = true;
  settings.swapOnFirstLayer         = true;
  settings.lookahead                = false;
  mapper->map(settings);
  mapper->dumpResult("qubit_ride_along_mapped.qasm");
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  /*
  expected output (order of gates may vary):
  SWAP(5,6)
  SWAP(2,5)
  SWAP(0,1)
  SWAP(1,2)
  CX(2,3) [x5]
  CX(4,1) [x5]
  */
  EXPECT_EQ(result.output.swaps, 4);

  const double c4 = -std::log2(1 - e4);
  const double c3 = -std::log2(1 - e3);
  const double c1 = -std::log2(1 - e1);

  const double expectedFidelity = 3 * c4 + 3 * c3 + 3 * c4 + 3 * c3 + // SWAPs
                                  5 * c1 + 5 * c1;                    // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, SingleQubitsCompete) {
  Architecture      architecture{};
  const CouplingMap cm = {{0, 1}, {1, 0}, {1, 2}, {2, 1}};
  architecture.loadCouplingMap(3, cm);

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", 0.8);
  props.setSingleQubitErrorRate(1, "x", 0.1);
  props.setSingleQubitErrorRate(2, "x", 0.8);

  props.setTwoQubitErrorRate(0, 1, 0.1);
  props.setTwoQubitErrorRate(1, 0, 0.1);
  props.setTwoQubitErrorRate(1, 2, 0.1);
  props.setTwoQubitErrorRate(2, 1, 0.1);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{3, 3};
  qc.x(0);
  qc.x(2);

  for (size_t i = 0; i < 3; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.admissibleHeuristic      = true;
  settings.layering                 = Layering::Disjoint2qBlocks;
  settings.initialLayout            = InitialLayout::Identity;
  settings.considerFidelity         = true;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  settings.verbose                  = true;
  settings.swapOnFirstLayer         = true;
  settings.lookahead                = false;
  mapper->map(settings);
  mapper->dumpResult("single_qubits_compete.qasm");
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 1);
  EXPECT_EQ(result.output.swaps, 1);

  const double expectedFidelity = -3 * std::log2(1 - 0.1) // SWAPs
                                  - std::log2(1 - 0.8) -
                                  std::log2(1 - 0.1); // Xs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);
}

TEST(HeuristicTestFidelity, LayerSplitting) {
  Architecture      architecture{};
  const CouplingMap cm = {
      {0, 1}, {1, 0}, {1, 2},  {2, 1},  {2, 3},   {3, 2},

      {0, 4}, {4, 0}, {1, 5},  {5, 1},  {2, 6},   {6, 2},  {3, 7},  {7, 3},

      {4, 5}, {5, 4}, {5, 6},  {6, 5},  {6, 7},   {7, 6},

      {4, 8}, {8, 4}, {5, 9},  {9, 5},  {6, 10},  {10, 6}, {7, 11}, {11, 7},

      {8, 9}, {9, 8}, {9, 10}, {10, 9}, {10, 11}, {11, 10}};
  architecture.loadCouplingMap(12, cm);

  const double e5 = 0.99;
  const double e4 = 0.9;
  const double e3 = 0.5;
  const double e2 = 0.4;
  const double e1 = 0.1;
  const double e0 = 0.01;

  auto props = Architecture::Properties();
  props.setSingleQubitErrorRate(0, "x", e5);
  props.setSingleQubitErrorRate(1, "x", e5);
  props.setSingleQubitErrorRate(2, "x", e5);
  props.setSingleQubitErrorRate(3, "x", e5);
  props.setSingleQubitErrorRate(4, "x", e3);
  props.setSingleQubitErrorRate(5, "x", e3);
  props.setSingleQubitErrorRate(6, "x", e5);
  props.setSingleQubitErrorRate(7, "x", e3);
  props.setSingleQubitErrorRate(8, "x", e2);
  props.setSingleQubitErrorRate(9, "x", e1);
  props.setSingleQubitErrorRate(10, "x", e5);
  props.setSingleQubitErrorRate(11, "x", e2);

  props.setTwoQubitErrorRate(0, 1, e4);
  props.setTwoQubitErrorRate(1, 0, e4);
  props.setTwoQubitErrorRate(1, 2, e4);
  props.setTwoQubitErrorRate(2, 1, e4);
  props.setTwoQubitErrorRate(2, 3, e1);
  props.setTwoQubitErrorRate(3, 2, e1);

  props.setTwoQubitErrorRate(0, 4, e5);
  props.setTwoQubitErrorRate(4, 0, e5);
  props.setTwoQubitErrorRate(1, 5, e5);
  props.setTwoQubitErrorRate(5, 1, e5);
  props.setTwoQubitErrorRate(2, 6, e4);
  props.setTwoQubitErrorRate(6, 2, e4);
  props.setTwoQubitErrorRate(3, 7, e5);
  props.setTwoQubitErrorRate(7, 3, e5);

  props.setTwoQubitErrorRate(4, 5, e3);
  props.setTwoQubitErrorRate(5, 4, e3);
  props.setTwoQubitErrorRate(5, 6, e5);
  props.setTwoQubitErrorRate(6, 5, e5);
  props.setTwoQubitErrorRate(6, 7, e5);
  props.setTwoQubitErrorRate(7, 6, e5);

  props.setTwoQubitErrorRate(4, 8, e0);
  props.setTwoQubitErrorRate(8, 4, e0);
  props.setTwoQubitErrorRate(5, 9, e3);
  props.setTwoQubitErrorRate(9, 5, e3);
  props.setTwoQubitErrorRate(6, 10, e1);
  props.setTwoQubitErrorRate(10, 6, e1);
  props.setTwoQubitErrorRate(7, 11, e3);
  props.setTwoQubitErrorRate(11, 7, e3);

  props.setTwoQubitErrorRate(8, 9, e5);
  props.setTwoQubitErrorRate(9, 8, e5);
  props.setTwoQubitErrorRate(9, 10, e4);
  props.setTwoQubitErrorRate(10, 9, e4);
  props.setTwoQubitErrorRate(10, 11, e5);
  props.setTwoQubitErrorRate(11, 10, e5);

  architecture.loadProperties(props);

  qc::QuantumComputation qc{12, 12};

  for (std::size_t i = 0; i < 50; ++i) {
    qc.x(4);
  }
  qc.x(5);
  qc.x(7);
  for (std::size_t i = 0; i < 5; ++i) {
    qc.cx(qc::Control{3}, 0);
    qc.cx(qc::Control{9}, 2);
  }

  for (size_t i = 0; i < 12; ++i) {
    qc.measure(static_cast<qc::Qubit>(i), i);
  }

  auto mapper = std::make_unique<HeuristicMapper>(qc, architecture);

  Configuration settings{};
  settings.verbose                  = true;
  settings.admissibleHeuristic      = true;
  settings.layering                 = Layering::Disjoint2qBlocks;
  settings.initialLayout            = InitialLayout::Identity;
  settings.considerFidelity         = true;
  settings.preMappingOptimizations  = false;
  settings.postMappingOptimizations = false;
  settings.swapOnFirstLayer         = true;
  settings.lookahead                = false;
  settings.automaticLayerSplits     = true;
  settings.automaticLayerSplitsNodeLimit =
      1; // force splittings after 1st expanded node until layers are
         // unsplittable
  settings.dataLoggingPath = "test_log/";
  mapper->map(settings);
  mapper->dumpResult("simple_grid_mapped.qasm");
  mapper->printResult(std::cout);

  auto& result = mapper->getResults();
  EXPECT_EQ(result.input.layers, 5); // originally 1 but split into 5 during A*
  /*
  expected output:
  === layout (placing of logical qubits in physical grid):
       0  1  2  3
       4  5  6  7
       8  9 10 11
  SWAP (4,5)
  SWAP (5,9)
  === layout:
       0  1  2  3
       5  9  6  7
       8  4 10 11
  X(9) [x50]   (originally X(4))

  === layout:
       0  1  2  3
       5  9  6  7
       8  4 10 11
  X(7)   (originally X(7))

  SWAP (5,9)
  SWAP (9,10)
  SWAP (2,6)
  === layout:
       0  1  6  3
       5  4  2  7
       8 10  9 11
  CX(6,10) [x5]   (originally CX(2,9))

  SWAP (4,8)
  === layout:
       0  1  6  3
       8  4  2  7
       5 10  9 11
  X(8)   (originally X(5))

  SWAP (0,1)
  SWAP (1,2)
  === layout:
       1  6  0  3
       8  4  2  7
       5 10  9 11
  CX(2,3) [x5]   (originally CX(0,3))
  */
  EXPECT_EQ(result.output.swaps, 8);

  const double c4 = -std::log2(1 - e4);
  const double c3 = -std::log2(1 - e3);
  const double c2 = -std::log2(1 - e2);
  const double c1 = -std::log2(1 - e1);
  const double c0 = -std::log2(1 - e0);

  const double expectedFidelity =
      3 * (c3 + c3 + c3 + c4 + c4 + c0 + c4 + c4) + // SWAPs
      50 * c1 + c3 + c2 +                           // Xs
      5 * c1 + 5 * c1;                              // CXs
  EXPECT_NEAR(result.output.totalLogFidelity, expectedFidelity, 1e-6);

  // check data log
  const std::array<std::string, 4> layerNodeFilePaths = {
      "nodes_layer_0.presplit-0.csv", "nodes_layer_0.presplit-1.csv",
      "nodes_layer_1.presplit-0.csv", "nodes_layer_3.presplit-0.csv"};
  for (const auto& path : layerNodeFilePaths) {
    auto layerNodeFile = std::ifstream(settings.dataLoggingPath + path);
    if (!layerNodeFile.is_open()) {
      FAIL() << "Could not open file " << settings.dataLoggingPath << path;
    }
    std::string line;
    while (std::getline(layerNodeFile, line)) {
      if (line.empty()) {
        continue;
      }
      std::string       col;
      std::stringstream lineStream(line);
      if (!std::getline(lineStream, col, ';')) {
        FAIL() << "Missing value for node id in " << settings.dataLoggingPath
               << path;
      }
      if (std::getline(lineStream, col, ';')) {
        if (std::stoull(col) != 0) {
          // should only contain root node and its direct children
          FAIL() << "Unexpected value for parent node id in "
                 << settings.dataLoggingPath << path;
        }
      } else {
        FAIL() << "Missing value for parent node id in "
               << settings.dataLoggingPath << path;
      }
    }
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_0.presplit-0.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_0.presplit-0.json"
           << " does not exist";
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_0.presplit-1.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_0.presplit-1.json"
           << " does not exist";
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_1.presplit-0.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_1.presplit-0.json"
           << " does not exist";
  }
  if (!std::filesystem::exists(settings.dataLoggingPath +
                               "layer_3.presplit-0.json")) {
    FAIL() << "File " << settings.dataLoggingPath << "layer_3.presplit-0.json"
           << " does not exist";
  }
}
