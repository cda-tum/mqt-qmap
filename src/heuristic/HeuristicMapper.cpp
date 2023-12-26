//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "heuristic/HeuristicMapper.hpp"

#include "utils.hpp"

#include <chrono>

void HeuristicMapper::map(const Configuration& configuration) {
  if (configuration.dataLoggingEnabled()) {
    dataLogger = std::make_unique<DataLogger>(configuration.dataLoggingPath,
                                              *architecture, qc);
  }
  results        = MappingResults{};
  results.config = configuration;
  auto& config   = results.config;
  if (config.layering == Layering::OddGates ||
      config.layering == Layering::QubitTriangle) {
    throw QMAPException("Layering strategy " + toString(config.layering) +
                        " not suitable for heuristic mapper!");
  }
  if (config.considerFidelity && !architecture->isFidelityAvailable()) {
    throw QMAPException("No calibration data available for this architecture!");
  }
  if (config.considerFidelity && config.lookahead) {
    throw QMAPException("Lookahead is not yet supported for heuristic mapper "
                        "using fidelity-aware mapping!");
  }
  if (config.considerFidelity && config.teleportationQubits > 0) {
    throw QMAPException("Teleportation is not yet supported for heuristic "
                        "mapper using fidelity-aware mapping!");
  }
  const auto start = std::chrono::steady_clock::now();
  initResults();

  // perform pre-mapping optimizations
  preMappingOptimizations(config);

  createLayers();
  if (config.verbose) {
    std::clog << "Teleportation qubits: " << config.teleportationQubits << "\n";
    printLayering(std::clog);
  }

  createInitialMapping();
  if (config.verbose) {
    printLocations(std::clog);
    printQubits(std::clog);
  }

  for (std::size_t i = 0; i < config.iterativeBidirectionalRoutingPasses; ++i) {
    if (config.verbose) {
      std::clog << std::endl
                << "Iterative bidirectional routing (forward pass " << i
                << "):" << std::endl;
    }
    pseudoRouteCircuit(false);
    if (config.verbose) {
      std::clog << std::endl
                << "Iterative bidirectional routing (backward pass " << i
                << "):" << std::endl;
    }
    pseudoRouteCircuit(true);

    if (config.verbose) {
      std::clog << std::endl << "Main routing:" << std::endl;
    }
  }

  routeCircuit();

  postMappingOptimizations(config);
  countGates(qcMapped, results.output);
  finalizeMappedCircuit();

  const auto                          end  = std::chrono::steady_clock::now();
  const std::chrono::duration<double> diff = end - start;
  results.time                             = diff.count();
  results.timeout                          = false;

  if (config.dataLoggingEnabled()) {
    dataLogger->logOutputCircuit(qcMapped);
    dataLogger->logMappingResult(results);
    dataLogger->close();
  }
}

void HeuristicMapper::staticInitialMapping() {
  for (const auto& gate : layers.at(0U)) {
    if (gate.singleQubit()) {
      continue;
    }

    for (const auto& [q0, q1] : architecture->getCouplingMap()) {
      if (qubits.at(q0) == DEFAULT_POSITION &&
          qubits.at(q1) == DEFAULT_POSITION) {
        qubits.at(q0) = gate.control;
        qubits.at(q1) = static_cast<std::int16_t>(gate.target);
        locations.at(static_cast<std::uint16_t>(gate.control)) =
            static_cast<std::int16_t>(q0);
        locations.at(gate.target)     = static_cast<std::int16_t>(q1);
        qcMapped.initialLayout.at(q0) = static_cast<qc::Qubit>(gate.control);
        qcMapped.initialLayout.at(q1) = static_cast<qc::Qubit>(gate.target);
        qcMapped.outputPermutation.at(q0) =
            static_cast<qc::Qubit>(gate.control);
        qcMapped.outputPermutation.at(q1) = static_cast<qc::Qubit>(gate.target);
        break;
      }
    }
  }

  // assign remaining logical qubits
  for (qc::Qubit i = 0U; i < architecture->getNqubits(); ++i) {
    if (qc.initialLayout.count(i) > 0 && locations.at(i) == DEFAULT_POSITION) {
      for (qc::Qubit j = 0U; j < architecture->getNqubits(); ++j) {
        if (qubits.at(j) == DEFAULT_POSITION) {
          locations.at(i)                  = static_cast<std::int16_t>(j);
          qubits.at(j)                     = static_cast<std::int16_t>(i);
          qcMapped.initialLayout.at(j)     = i;
          qcMapped.outputPermutation.at(j) = i;
          break;
        }
      }
    }
  }
}

void HeuristicMapper::createInitialMapping() {
  auto& config = results.config;

  if (layers.empty()) {
    return;
  }

  if (config.teleportationQubits > 0) {
    std::mt19937_64 mt;
    if (config.teleportationSeed == 0) {
      std::array<std::mt19937_64::result_type, std::mt19937_64::state_size>
                         randomData{};
      std::random_device rd;
      std::generate(std::begin(randomData), std::end(randomData),
                    [&rd]() { return rd(); });
      std::seed_seq seeds(std::begin(randomData), std::end(randomData));
      mt.seed(seeds);
    } else {
      mt.seed(config.teleportationSeed);
    }

    std::uniform_int_distribution<> dis(0, architecture->getNqubits() - 1);

    for (std::size_t i = 0; i < config.teleportationQubits; i += 2) {
      Edge e{};
      do { // NOLINT(cppcoreguidelines-avoid-do-while)
        auto it = std::begin(architecture->getCouplingMap());
        std::advance(it, dis(mt));
        e = *it;
      } while (qubits.at(e.first) != -1 || qubits.at(e.second) != -1);
      const auto teleportationQubit    = qc.getNqubits() + i;
      locations.at(teleportationQubit) = static_cast<std::int16_t>(e.first);
      locations.at(teleportationQubit + 1) =
          static_cast<std::int16_t>(e.second);
      qubits.at(e.first)  = static_cast<std::int16_t>(teleportationQubit);
      qubits.at(e.second) = static_cast<std::int16_t>(teleportationQubit + 1);
    }

    if (config.teleportationFake) {
      config.teleportationQubits = 0;
    }
  }

  switch (config.initialLayout) {
  case InitialLayout::Identity:
    for (qc::Qubit i = 0; i < architecture->getNqubits(); ++i) {
      if (qc.initialLayout.count(i) > 0) {
        locations.at(i) = static_cast<std::int16_t>(i);
        qubits.at(i)    = static_cast<std::int16_t>(i);
      }
    }
    break;
  case InitialLayout::Static:
    staticInitialMapping();
    break;
  case InitialLayout::Dynamic:
    // nothing to be done here
    break;

    // TODO: Design strategy that maps most used qubit to most connected qubits
    // on architecture
  }
}

void HeuristicMapper::mapUnmappedGates(std::size_t layer) {
  if (results.config.considerFidelity) {
    for (std::size_t q = 0; q < singleQubitMultiplicities.at(layer).size();
         ++q) {
      if (singleQubitMultiplicities.at(layer).at(q) == 0) {
        continue;
      }
      if (locations.at(q) == DEFAULT_POSITION) {
        // TODO: consider fidelity
        // map to first free physical qubit
        for (std::uint16_t physQbit = 0; physQbit < architecture->getNqubits();
             ++physQbit) {
          if (qubits.at(physQbit) == -1) {
            locations.at(q)     = static_cast<std::int16_t>(physQbit);
            qubits.at(physQbit) = static_cast<std::int16_t>(q);
            break;
          }
        }
      }
    }
  }

  for (const auto& [logEdge, _] : twoQubitMultiplicities.at(layer)) {
    const auto& [q1, q2] = logEdge;

    auto q1Location = locations.at(q1);
    auto q2Location = locations.at(q2);

    if (q1Location == DEFAULT_POSITION && q2Location == DEFAULT_POSITION) {
      std::set<Edge> possibleEdges{};
      // gather all edges in the architecture for which both qubits are unmapped
      for (const auto& edge : architecture->getCouplingMap()) {
        if (qubits.at(edge.first) == DEFAULT_POSITION &&
            qubits.at(edge.second) == DEFAULT_POSITION) {
          possibleEdges.emplace(edge);
        }
      }
      std::pair<std::uint16_t, std::uint16_t> chosenEdge;

      if (possibleEdges.empty()) {
        // map to 2 qubits with minimal distance
        double bestScore = std::numeric_limits<int>::max();

        for (std::uint16_t i = 0; i < architecture->getNqubits(); i++) {
          for (std::uint16_t j = i + 1; j < architecture->getNqubits(); j++) {
            if (qubits.at(i) == DEFAULT_POSITION &&
                qubits.at(j) == DEFAULT_POSITION) {
              const double dist = architecture->distance(i, j);
              if (dist < bestScore) {
                bestScore  = dist;
                chosenEdge = std::make_pair(i, j);
              }
            }
          }
        }
      } else {
        chosenEdge = *possibleEdges.begin();
      }
      // TODO: Consider fidelity here if available. The best available edge
      // should be chosen
      locations.at(q1) = static_cast<std::int16_t>(chosenEdge.first);
      locations.at(q2) = static_cast<std::int16_t>(chosenEdge.second);
      qubits.at(chosenEdge.first)  = static_cast<std::int16_t>(q1);
      qubits.at(chosenEdge.second) = static_cast<std::int16_t>(q2);
      qc::QuantumComputation::findAndSWAP(q1, chosenEdge.first,
                                          qcMapped.initialLayout);
      qc::QuantumComputation::findAndSWAP(q2, chosenEdge.second,
                                          qcMapped.initialLayout);
      qc::QuantumComputation::findAndSWAP(q1, chosenEdge.first,
                                          qcMapped.outputPermutation);
      qc::QuantumComputation::findAndSWAP(q2, chosenEdge.second,
                                          qcMapped.outputPermutation);
    } else if (q1Location == DEFAULT_POSITION) {
      mapToMinDistance(q2, q1);
    } else if (q2Location == DEFAULT_POSITION) {
      mapToMinDistance(q1, q2);
    }
  }
}

void HeuristicMapper::mapToMinDistance(const std::uint16_t source,
                                       const std::uint16_t target) {
  auto                         min = std::numeric_limits<double>::max();
  std::optional<std::uint16_t> pos = std::nullopt;
  for (std::uint16_t i = 0; i < architecture->getNqubits(); ++i) {
    if (qubits.at(i) == DEFAULT_POSITION) {
      // TODO: Consider fidelity here if available
      auto distance = distanceOnArchitectureOfPhysicalQubits(
          static_cast<std::uint16_t>(locations.at(source)), i);
      if (distance < min) {
        min = distance;
        pos = i;
      }
    }
  }
  assert(pos.has_value());
  qubits.at(*pos)      = static_cast<std::int16_t>(target);
  locations.at(target) = static_cast<std::int16_t>(*pos);
  qc::QuantumComputation::findAndSWAP(target, *pos, qcMapped.initialLayout);
  qc::QuantumComputation::findAndSWAP(target, *pos, qcMapped.outputPermutation);
}

void HeuristicMapper::pseudoRouteCircuit(bool reverse) {
  // save original global data for restoring it later
  const auto originalResults                   = results;
  const auto originalLayers                    = layers;
  const auto originalSingleQubitMultiplicities = singleQubitMultiplicities;
  const auto originalTwoQubitMultiplicities    = twoQubitMultiplicities;
  const auto originalActiveQubits              = activeQubits;
  const auto originalActiveQubits1QGates       = activeQubits1QGates;
  const auto originalActiveQubits2QGates       = activeQubits2QGates;

  auto& config           = results.config;
  config.dataLoggingPath = ""; // disable data logging for pseudo routing
  config.debug           = false;

  for (std::size_t i = 0; i < layers.size(); ++i) {
    auto       layerIndex = (reverse ? layers.size() - i - 1 : i);
    const Node result     = aStarMap(layerIndex, reverse);

    qubits    = result.qubits;
    locations = result.locations;

    if (config.verbose) {
      printLocations(std::clog);
      printQubits(std::clog);
    }
  }

  // restore original global data
  results                   = originalResults;
  layers                    = originalLayers;
  singleQubitMultiplicities = originalSingleQubitMultiplicities;
  twoQubitMultiplicities    = originalTwoQubitMultiplicities;
  activeQubits              = originalActiveQubits;
  activeQubits1QGates       = originalActiveQubits1QGates;
  activeQubits2QGates       = originalActiveQubits2QGates;
}

void HeuristicMapper::routeCircuit() {
  auto& config = results.config;

  std::size_t              gateidx = 0;
  std::vector<std::size_t> gatesToAdjust{};
  results.output.gates = 0U;
  for (std::size_t layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
    const Node result = aStarMap(layerIndex, false);

    qubits    = result.qubits;
    locations = result.locations;

    if (config.verbose) {
      printLocations(std::clog);
      printQubits(std::clog);
    }

    if (layerIndex != 0 && config.addBarriersBetweenLayers) {
      qcMapped.barrier();
      gateidx++;
    }

    // initial layer needs no swaps
    if (layerIndex != 0 || config.swapOnFirstLayer) {
      for (const auto& swaps : result.swaps) {
        for (const auto& swap : swaps) {
          if (swap.op == qc::SWAP) {
            if (config.verbose) {
              std::clog << "SWAP: " << swap.first << " <-> " << swap.second
                        << "\n";
            }
            if (!architecture->isEdgeConnected({swap.first, swap.second}) &&
                !architecture->isEdgeConnected({swap.second, swap.first})) {
              throw QMAPException(
                  "Invalid SWAP: " + std::to_string(swap.first) + "<->" +
                  std::to_string(swap.second));
            }
            qcMapped.swap(swap.first, swap.second);
            results.output.swaps++;
          } else if (swap.op == qc::Teleportation) {
            if (config.verbose) {
              std::clog << "TELE: " << swap.first << " <-> " << swap.second
                        << "\n";
            }
            qcMapped.emplace_back<qc::StandardOperation>(
                qcMapped.getNqubits(),
                qc::Targets{static_cast<qc::Qubit>(swap.first),
                            static_cast<qc::Qubit>(swap.second),
                            static_cast<qc::Qubit>(swap.middleAncilla)},
                qc::Teleportation);
            results.output.teleportations++;
          }
          gateidx++;
        }
      }
    }

    // add gates of the layer to circuit
    for (const auto& gate : layers.at(layerIndex)) {
      auto* op = dynamic_cast<qc::StandardOperation*>(gate.op);
      if (op == nullptr) {
        throw QMAPException(
            "Cast to StandardOperation not possible during mapping. Check that "
            "circuit contains only StandardOperations");
      }

      if (gate.singleQubit()) {
        if (locations.at(gate.target) == DEFAULT_POSITION) {
          qcMapped.emplace_back<qc::StandardOperation>(
              qcMapped.getNqubits(), gate.target, op->getType(),
              op->getParameter());
          gatesToAdjust.push_back(gateidx);
          gateidx++;
        } else {
          qcMapped.emplace_back<qc::StandardOperation>(
              qcMapped.getNqubits(), locations.at(gate.target), op->getType(),
              op->getParameter());
          gateidx++;
        }
      } else {
        const Edge cnot = {
            locations.at(static_cast<std::uint16_t>(gate.control)),
            locations.at(gate.target)};
        if (!architecture->isEdgeConnected(cnot)) {
          const Edge reversed = {cnot.second, cnot.first};
          if (!architecture->isEdgeConnected(reversed)) {
            throw QMAPException(
                "Invalid CNOT: " + std::to_string(reversed.first) + "-" +
                std::to_string(reversed.second));
          }
          qcMapped.h(reversed.first);
          qcMapped.h(reversed.second);
          qcMapped.cx(qc::Control{static_cast<qc::Qubit>(reversed.first)},
                      reversed.second);
          qcMapped.h(reversed.second);
          qcMapped.h(reversed.first);

          results.output.directionReverse++;
          gateidx += 5;
        } else {
          qcMapped.cx(qc::Control{static_cast<qc::Qubit>(cnot.first)},
                      cnot.second);
          gateidx++;
        }
      }
    }
  }

  if (config.debug && results.heuristicBenchmark.expandedNodes > 0) {
    auto& benchmark = results.heuristicBenchmark;
    benchmark.timePerNode /= static_cast<double>(benchmark.expandedNodes);
    benchmark.averageBranchingFactor =
        static_cast<double>(benchmark.generatedNodes - layers.size()) /
        static_cast<double>(benchmark.expandedNodes);
    for (const auto& layer : results.layerHeuristicBenchmark) {
      benchmark.effectiveBranchingFactor +=
          layer.effectiveBranchingFactor *
          (static_cast<double>(layer.expandedNodes) /
           static_cast<double>(benchmark.expandedNodes));
    }
  }

  // infer output permutation from qubit locations
  qcMapped.outputPermutation.clear();
  for (std::size_t i = 0U; i < architecture->getNqubits(); ++i) {
    if (const auto lq = qubits.at(i); lq != -1) {
      const auto logicalQubit = static_cast<qc::Qubit>(lq);
      // check whether this is a qubit from the original circuit
      if (logicalQubit < qc.getNqubits()) {
        qcMapped.outputPermutation[static_cast<qc::Qubit>(i)] =
            static_cast<qc::Qubit>(qubits.at(i));
      } else {
        qcMapped.setLogicalQubitGarbage(logicalQubit);
      }
    }
  }

  // fix single qubit gates
  if (!gatesToAdjust.empty()) {
    gateidx--; // index of last operation
    for (auto it = qcMapped.rbegin(); it != qcMapped.rend(); ++it, --gateidx) {
      auto* op = dynamic_cast<qc::StandardOperation*>(it->get());
      if (op == nullptr) {
        throw QMAPException(
            "Cast to StandardOperation not possible during mapping. Check that "
            "circuit contains only StandardOperations");
      }
      if (op->getType() == qc::SWAP) {
        const auto q0                     = qubits.at(op->getTargets().at(0));
        const auto q1                     = qubits.at(op->getTargets().at(1));
        qubits.at(op->getTargets().at(0)) = q1;
        qubits.at(op->getTargets().at(1)) = q0;

        if (q0 != DEFAULT_POSITION) {
          locations.at(static_cast<std::size_t>(q0)) =
              static_cast<std::int16_t>(op->getTargets().at(1));
        }
        if (q1 != DEFAULT_POSITION) {
          locations.at(static_cast<std::size_t>(q1)) =
              static_cast<std::int16_t>(op->getTargets().at(0));
        }
      }
      if (!gatesToAdjust.empty() && gatesToAdjust.back() == gateidx) {
        gatesToAdjust.pop_back();
        auto target         = op->getTargets().at(0);
        auto targetLocation = locations.at(target);

        if (targetLocation == -1) {
          // qubit only occurs in single qubit gates, can be mapped to an
          // arbitrary free qubit
          std::uint16_t loc = 0;
          while (qubits.at(loc) != DEFAULT_POSITION) {
            ++loc;
          }
          locations.at(target) = static_cast<std::int16_t>(loc);
          qubits.at(loc)       = static_cast<std::int16_t>(target);
          op->setTargets({static_cast<qc::Qubit>(loc)});
          qcMapped.initialLayout.at(target)                       = loc;
          qcMapped.outputPermutation[static_cast<qc::Qubit>(loc)] = target;
          qcMapped.garbage.at(loc)                                = false;
        } else {
          op->setTargets({static_cast<qc::Qubit>(targetLocation)});
        }
      }
    }
  }

  // mark every qubit that is not mapped to a logical qubit as garbage
  std::size_t count = 0U;
  for (std::size_t i = 0U; i < architecture->getNqubits(); ++i) {
    if (const auto lq = qubits.at(i); lq == -1) {
      qcMapped.setLogicalQubitGarbage(
          static_cast<qc::Qubit>(qc.getNqubits() + count));
      ++count;
    }
  }
}

HeuristicMapper::Node HeuristicMapper::aStarMap(size_t layer, bool reverse) {
  auto&      config           = results.config;
  const bool considerFidelity = config.considerFidelity;
  nextNodeId                  = 0;

  const std::unordered_set<std::uint16_t>& consideredQubits =
      (considerFidelity ? activeQubits[layer] : activeQubits2QGates[layer]);
  const SingleQubitMultiplicity& singleQubitMultiplicity =
      singleQubitMultiplicities.at(layer);
  const TwoQubitMultiplicity& twoQubitMultiplicity =
      twoQubitMultiplicities.at(layer);
  Node node(nextNodeId++, considerFidelity, config.admissibleHeuristic);
  Node bestDoneNode(0);
  bool done = false;

  mapUnmappedGates(layer);

  node.locations = locations;
  node.qubits    = qubits;
  node.recalculateFixedCost(*architecture, singleQubitMultiplicity,
                            twoQubitMultiplicity);
  node.updateHeuristicCost(*architecture, singleQubitMultiplicity,
                           twoQubitMultiplicity, consideredQubits);

  if (config.dataLoggingEnabled()) {
    dataLogger->logSearchNode(layer, node.id, node.parent, node.costFixed,
                              node.costHeur, node.lookaheadPenalty, node.qubits,
                              node.done, node.swaps, node.depth);
  }
  nodes.push(node);

  const auto  start         = std::chrono::steady_clock::now();
  std::size_t expandedNodes = 0;

  const bool splittable =
      config.automaticLayerSplits ? isLayerSplittable(layer) : false;

  while (!nodes.empty() && (!done || nodes.top().getTotalCost() <
                                         bestDoneNode.getTotalFixedCost())) {
    if (splittable && expandedNodes >= config.automaticLayerSplitsNodeLimit) {
      if (config.dataLoggingEnabled()) {
        qc::CompoundOperation compOp(architecture->getNqubits());
        for (const auto& gate : layers.at(layer)) {
          compOp.emplace_back(gate.op->clone());
        }

        dataLogger->logFinalizeLayer(layer, compOp, singleQubitMultiplicity,
                                     twoQubitMultiplicity, qubits, 0, 0, 0, 0,
                                     {}, {}, 0);
        dataLogger->splitLayer();
      }
      splitLayer(layer, *architecture);
      if (config.verbose) {
        std::clog << "Split layer" << std::endl;
      }
      // recursively restart search with newly split layer
      // (step to the end of the circuit, if reverse mapping is active, since
      // the split layer is inserted in this direction, otherwise 1 layer would
      // be skipped)
      return aStarMap(reverse ? layer + 1 : layer, reverse);
    }
    Node current = nodes.top();
    if (current.done) {
      if (!done ||
          current.getTotalFixedCost() < bestDoneNode.getTotalFixedCost()) {
        bestDoneNode = current;
      }
      done = true;
      if (!considerFidelity) {
        break;
      }
    }
    nodes.pop();
    expandNode(consideredQubits, current, layer);
    ++expandedNodes;
  }

  if (!done) {
    throw QMAPException("No viable mapping found.");
  }

  Node result = bestDoneNode;
  if (config.debug) {
    const auto end = std::chrono::steady_clock::now();
    results.layerHeuristicBenchmark.emplace_back();
    auto layerResultsIt           = results.layerHeuristicBenchmark.rbegin();
    layerResultsIt->expandedNodes = expandedNodes;
    results.heuristicBenchmark.expandedNodes += expandedNodes;

    layerResultsIt->solutionDepth = result.depth;

    const std::chrono::duration<double> diff = end - start;
    results.heuristicBenchmark.timePerNode += diff.count();

    layerResultsIt->generatedNodes =
        layerResultsIt->expandedNodes + nodes.size();
    results.heuristicBenchmark.generatedNodes += layerResultsIt->generatedNodes;

    if (layerResultsIt->expandedNodes > 0) {
      layerResultsIt->timePerNode =
          diff.count() / static_cast<double>(layerResultsIt->expandedNodes);
      layerResultsIt->averageBranchingFactor =
          static_cast<double>(layerResultsIt->generatedNodes - 1) /
          static_cast<double>(layerResultsIt->expandedNodes);
    }

    layerResultsIt->effectiveBranchingFactor = computeEffectiveBranchingRate(
        layerResultsIt->expandedNodes + 1, result.depth);
  }

  if (config.dataLoggingEnabled()) {
    qc::CompoundOperation compOp(architecture->getNqubits());
    for (const auto& gate : layers.at(layer)) {
      compOp.emplace_back(gate.op->clone());
    }

    dataLogger->logFinalizeLayer(
        layer, compOp, singleQubitMultiplicities.at(layer),
        twoQubitMultiplicities.at(layer), qubits, result.id, result.costFixed,
        result.costHeur, result.lookaheadPenalty, result.qubits, result.swaps,
        result.depth);
  }

  // clear nodes
  while (!nodes.empty()) {
    nodes.pop();
  }

  return result;
}

void HeuristicMapper::expandNode(
    const std::unordered_set<std::uint16_t>& consideredQubits, Node& node,
    std::size_t layer) {
  std::vector<std::vector<bool>> usedSwaps;
  usedSwaps.reserve(architecture->getNqubits());
  for (int p = 0; p < architecture->getNqubits(); ++p) {
    usedSwaps.emplace_back(architecture->getNqubits());
  }

  // set up new teleportation qubits
  std::set<Edge> perms = architecture->getCouplingMap();
  architecture->getCurrentTeleportations().clear();
  architecture->getTeleportationQubits().clear();
  for (std::size_t i = 0; i < results.config.teleportationQubits; i += 2) {
    architecture->getTeleportationQubits().emplace_back(
        node.locations.at(qc.getNqubits() + i),
        node.locations.at(qc.getNqubits() + i + 1));
    Edge e;
    for (auto const& g : architecture->getCouplingMap()) {
      if (g.first == node.locations.at(qc.getNqubits() + i) &&
          g.second != node.locations.at(qc.getNqubits() + i + 1)) {
        e.first  = g.second;
        e.second = static_cast<std::uint16_t>(
            node.locations.at(qc.getNqubits() + i + 1));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.second == node.locations.at(qc.getNqubits() + i) &&
          g.first != node.locations.at(qc.getNqubits() + i + 1)) {
        e.first  = g.first;
        e.second = static_cast<std::uint16_t>(
            node.locations.at(qc.getNqubits() + i + 1));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.first == node.locations.at(qc.getNqubits() + i + 1) &&
          g.second != node.locations.at(qc.getNqubits() + i)) {
        e.first = g.second;
        e.second =
            static_cast<std::uint16_t>(node.locations.at(qc.getNqubits() + i));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.second == node.locations.at(qc.getNqubits() + i + 1) &&
          g.first != node.locations.at(qc.getNqubits() + i)) {
        e.first = g.first;
        e.second =
            static_cast<std::uint16_t>(node.locations.at(qc.getNqubits() + i));
        architecture->getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
    }
  }

  for (const auto& q : consideredQubits) {
    for (const auto& edge : perms) {
      if (edge.first == node.locations.at(q) ||
          edge.second == node.locations.at(q)) {
        auto q1 = node.qubits.at(edge.first);
        auto q2 = node.qubits.at(edge.second);
        if (q2 == -1 || q1 == -1) {
          expandNodeAddOneSwap(edge, node, layer, consideredQubits);
        } else if (!usedSwaps.at(static_cast<std::size_t>(q1))
                        .at(static_cast<std::size_t>(q2))) {
          usedSwaps.at(static_cast<std::size_t>(q1))
              .at(static_cast<std::size_t>(q2)) = true;
          usedSwaps.at(static_cast<std::size_t>(q2))
              .at(static_cast<std::size_t>(q1)) = true;
          expandNodeAddOneSwap(edge, node, layer, consideredQubits);
        }
      }
    }
  }
}

void HeuristicMapper::expandNodeAddOneSwap(
    const Edge& swap, Node& node, const std::size_t layer,
    const std::unordered_set<std::uint16_t>& consideredQubits) {
  const auto& config = results.config;

  Node newNode = Node(nextNodeId++, node.id, node.qubits, node.locations,
                      node.swaps, node.costFixed, node.depth + 1,
                      node.considerFidelity, node.admissibleHeuristic);

  if (architecture->isEdgeConnected(swap) ||
      architecture->isEdgeConnected(Edge{swap.second, swap.first})) {
    newNode.applySWAP(swap, *architecture, singleQubitMultiplicities.at(layer),
                      twoQubitMultiplicities.at(layer), consideredQubits);
  } else {
    newNode.applyTeleportation(swap, *architecture);
  }

  newNode.updateHeuristicCost(
      *architecture, singleQubitMultiplicities.at(layer),
      twoQubitMultiplicities.at(layer), consideredQubits);

  // calculate heuristics for the cost of the following layers
  if (config.lookahead) {
    lookahead(getNextLayer(layer), newNode);
  }

  nodes.push(newNode);
  if (results.config.dataLoggingEnabled()) {
    dataLogger->logSearchNode(layer, newNode.id, newNode.parent,
                              newNode.costFixed, newNode.costHeur,
                              newNode.lookaheadPenalty, newNode.qubits,
                              newNode.done, newNode.swaps, newNode.depth);
  }
}

void HeuristicMapper::lookahead(const std::size_t      layer,
                                HeuristicMapper::Node& node) {
  const auto& config    = results.config;
  auto        nextLayer = layer;
  double      factor    = config.firstLookaheadFactor;

  for (std::size_t i = 0; i < config.nrLookaheads; ++i) {
    if (nextLayer == std::numeric_limits<std::size_t>::max()) {
      break;
    }

    double penalty = 0.;
    for (const auto& gate : layers.at(nextLayer)) {
      if (gate.singleQubit()) {
        continue;
      }

      auto loc1 = node.locations.at(static_cast<std::uint16_t>(gate.control));
      auto loc2 = node.locations.at(gate.target);
      if (loc1 == DEFAULT_POSITION && loc2 == DEFAULT_POSITION) {
        // no penalty
      } else if (loc1 == DEFAULT_POSITION) {
        auto min = std::numeric_limits<double>::max();
        for (std::uint16_t j = 0; j < architecture->getNqubits(); ++j) {
          if (node.qubits.at(j) == DEFAULT_POSITION) {
            // TODO: Consider fidelity here if available
            min = std::min(min, distanceOnArchitectureOfPhysicalQubits(
                                    j, static_cast<std::uint16_t>(
                                           node.locations.at(gate.target))));
          }
        }
        penalty = heuristicAddition(penalty, min);
      } else if (loc2 == DEFAULT_POSITION) {
        auto min = std::numeric_limits<double>::max();
        for (std::uint16_t j = 0; j < architecture->getNqubits(); ++j) {
          if (node.qubits.at(j) == DEFAULT_POSITION) {
            // TODO: Consider fidelity here if available
            min = std::min(min,
                           distanceOnArchitectureOfPhysicalQubits(
                               static_cast<std::uint16_t>(node.locations.at(
                                   static_cast<std::uint16_t>(gate.control))),
                               j));
          }
        }
        penalty = heuristicAddition(penalty, min);
      } else {
        auto cost = architecture->distance(
            static_cast<std::uint16_t>(
                node.locations.at(static_cast<std::uint16_t>(gate.control))),
            static_cast<std::uint16_t>(node.locations.at(gate.target)));
        penalty = heuristicAddition(penalty, cost);
      }
    }

    node.lookaheadPenalty += factor * penalty;
    factor *= config.lookaheadFactor;
    nextLayer = getNextLayer(nextLayer); // TODO: consider single qubits here
                                         // for better fidelity lookahead
  }
}

void HeuristicMapper::Node::applySWAP(
    const Edge& swap, Architecture& arch,
    const SingleQubitMultiplicity& singleQubitGateMultiplicity,
    const TwoQubitMultiplicity&    twoQubitGateMultiplicity,
    const std::unordered_set<std::uint16_t>& consideredQubits) {
  ++nswaps;
  swaps.emplace_back();
  const auto q1 = qubits.at(swap.first);
  const auto q2 = qubits.at(swap.second);
  
  if (consideredQubits.find(swap.first) != consideredQubits.end() &&
      consideredQubits.find(swap.second) != consideredQubits.end()) {
    ++sharedSwaps;
  }

  qubits.at(swap.first)  = q2;
  qubits.at(swap.second) = q1;

  if (q1 != -1) {
    locations.at(static_cast<std::size_t>(q1)) =
        static_cast<std::int16_t>(swap.second);
  }
  if (q2 != -1) {
    locations.at(static_cast<std::size_t>(q2)) =
        static_cast<std::int16_t>(swap.first);
  }

  if (arch.isEdgeConnected(swap) ||
      arch.isEdgeConnected(Edge{swap.second, swap.first})) {
    swaps.back().emplace_back(swap.first, swap.second, qc::SWAP);
  } else {
    throw QMAPException("Something wrong in applySWAP.");
  }

  if (considerFidelity) {
    std::uint16_t q1Mult = 0;
    std::uint16_t q2Mult = 0;
    if (q1 != -1) {
      q1Mult = singleQubitGateMultiplicity.at(static_cast<std::size_t>(q1));
    }
    if (q2 != -1) {
      q2Mult = singleQubitGateMultiplicity.at(static_cast<std::size_t>(q2));
    }
    // accounting for fidelity difference of single qubit gates (two qubit
    // gates are handled in the heuristic)
    costFixed +=
        ((q2Mult - q1Mult) * arch.getSingleQubitFidelityCost(swap.first) +
         (q1Mult - q2Mult) * arch.getSingleQubitFidelityCost(swap.second));
    // adding cost of the swap gate itself
    costFixed += arch.getSwapFidelityCost(swap.first, swap.second);
    // add cost of newly validly mapped gates and
    // remove cost of now no longer validly mapped gates
    for (const auto& [edge, mult] : twoQubitGateMultiplicity) {
      auto [q3, q4] = edge;
      if (q3 == q1 || q3 == q2 || q4 == q1 || q4 == q2) {
        auto physQ3 = static_cast<std::uint16_t>(locations.at(q3));
        auto physQ4 = static_cast<std::uint16_t>(locations.at(q4));
        if (arch.isEdgeConnected(Edge{physQ3, physQ4}) ||
            arch.isEdgeConnected(Edge{physQ4, physQ3})) {
          // validly mapped now
          if (validMappedTwoQubitGates.find(edge) ==
              validMappedTwoQubitGates.end()) { // not mapped validly before
            costFixed +=
                mult.first * arch.getTwoQubitFidelityCost(physQ3, physQ4) +
                mult.second * arch.getTwoQubitFidelityCost(physQ4, physQ3);
            validMappedTwoQubitGates.emplace(edge);
          }
        } else { // not mapped validly now
          if (validMappedTwoQubitGates.find(edge) !=
              validMappedTwoQubitGates.end()) { // mapped validly before
            auto prevPhysQ3 = physQ3;
            if (prevPhysQ3 == swap.first) {
              prevPhysQ3 = swap.second;
            } else if (prevPhysQ3 == swap.second) {
              prevPhysQ3 = swap.first;
            }
            auto prevPhysQ4 = physQ4;
            if (prevPhysQ4 == swap.first) {
              prevPhysQ4 = swap.second;
            } else if (prevPhysQ4 == swap.second) {
              prevPhysQ4 = swap.first;
            }
            costFixed -= mult.first * arch.getTwoQubitFidelityCost(prevPhysQ3,
                                                                   prevPhysQ4) +
                         mult.second * arch.getTwoQubitFidelityCost(prevPhysQ4,
                                                                    prevPhysQ3);
            validMappedTwoQubitGates.erase(edge);
          }
        }
      }
    }
  } else {
    if (arch.bidirectional()) {
      costFixed += COST_BIDIRECTIONAL_SWAP;
    } else {
      costFixed += COST_UNIDIRECTIONAL_SWAP;
    }
  }
}

void HeuristicMapper::Node::applyTeleportation(const Edge&   swap,
                                               Architecture& arch) {
  nswaps++;
  swaps.emplace_back();
  const auto q1 = qubits.at(swap.first);
  const auto q2 = qubits.at(swap.second);

  qubits.at(swap.first)  = q2;
  qubits.at(swap.second) = q1;

  if (q1 != -1) {
    locations.at(static_cast<std::size_t>(q1)) =
        static_cast<std::int16_t>(swap.second);
  }
  if (q2 != -1) {
    locations.at(static_cast<std::size_t>(q2)) =
        static_cast<std::int16_t>(swap.first);
  }

  std::uint16_t middleAnc = std::numeric_limits<decltype(middleAnc)>::max();
  for (const auto& qpair : arch.getTeleportationQubits()) {
    if (swap.first == qpair.first || swap.second == qpair.first) {
      middleAnc = static_cast<std::uint16_t>(qpair.second);
    } else if (swap.first == qpair.second || swap.second == qpair.second) {
      middleAnc = static_cast<std::uint16_t>(qpair.first);
    }
  }

  if (middleAnc == std::numeric_limits<decltype(middleAnc)>::max()) {
    throw QMAPException("Teleportation between seemingly wrong qubits: " +
                        std::to_string(swap.first) + " <--> " +
                        std::to_string(swap.second));
  }

  std::uint16_t source = std::numeric_limits<decltype(source)>::max();
  std::uint16_t target = std::numeric_limits<decltype(target)>::max();
  if (arch.isEdgeConnected({swap.first, middleAnc}) ||
      arch.isEdgeConnected({middleAnc, swap.first})) {
    source = swap.first;
    target = swap.second;
  } else {
    source = swap.second;
    target = swap.first;
  }

  if (source == middleAnc || target == middleAnc) {
    std::clog << "FAIL: TELE " << source << " -(" << middleAnc << ")-> "
              << target << "\n";
    throw QMAPException("Overlap between source/target and middle "
                        "ancillary in teleportation.");
  }

  swaps.back().emplace_back(source, target, middleAnc, qc::Teleportation);

  costFixed += COST_TELEPORTATION;
}

void HeuristicMapper::Node::recalculateFixedCost(
    const Architecture&            arch,
    const SingleQubitMultiplicity& singleQubitGateMultiplicity,
    const TwoQubitMultiplicity&    twoQubitGateMultiplicity) {
  costFixed = 0;
  if (considerFidelity) {
    // adding costs of single qubit gates
    for (std::uint16_t i = 0U; i < arch.getNqubits(); ++i) {
      if (singleQubitGateMultiplicity.at(i) == 0) {
        continue;
      }
      costFixed += singleQubitGateMultiplicity.at(i) *
                   arch.getSingleQubitFidelityCost(
                       static_cast<std::uint16_t>(locations.at(i)));
    }
    // adding cost of the swap gates
    for (auto& swapNode : swaps) {
      for (auto& swap : swapNode) {
        if (swap.op == qc::SWAP) {
          costFixed += arch.getSwapFidelityCost(swap.first, swap.second);
        } else if (swap.op == qc::Teleportation) {
          throw QMAPException("Teleportation currently not supported for "
                              "noise-aware mapping");
        }
      }
    }
    validMappedTwoQubitGates.clear();
    // adding cost of two qubit gates that are already mapped next to each other
    for (const auto& edgeMultiplicity : twoQubitGateMultiplicity) {
      const auto& q1                   = edgeMultiplicity.first.first;
      const auto& q2                   = edgeMultiplicity.first.second;
      const auto& straightMultiplicity = edgeMultiplicity.second.first;
      const auto& reverseMultiplicity  = edgeMultiplicity.second.second;

      if (arch.isEdgeConnected(
              {static_cast<std::uint16_t>(locations.at(q1)),
               static_cast<std::uint16_t>(locations.at(q2))}) ||
          arch.isEdgeConnected({static_cast<std::uint16_t>(locations.at(q2)),
                                static_cast<std::uint16_t>(
                                    locations.at(q1))})) { // validly mapped
        costFixed += (straightMultiplicity *
                          arch.getTwoQubitFidelityCost(
                              static_cast<std::uint16_t>(locations.at(q1)),
                              static_cast<std::uint16_t>(locations.at(q2))) +
                      reverseMultiplicity *
                          arch.getTwoQubitFidelityCost(
                              static_cast<std::uint16_t>(locations.at(q2)),
                              static_cast<std::uint16_t>(locations.at(q1))));
        validMappedTwoQubitGates.emplace(q1, q2);
      }
    }
    // 2-qubit-gates not yet mapped next to each other are handled in the
    // heuristic
  } else {
    for (auto& swapNode : swaps) {
      for (auto& swap : swapNode) {
        if (swap.op == qc::SWAP) {
          if (arch.bidirectional()) {
            costFixed += COST_BIDIRECTIONAL_SWAP;
          } else {
            costFixed += COST_UNIDIRECTIONAL_SWAP;
          }
        } else if (swap.op == qc::Teleportation) {
          costFixed += COST_TELEPORTATION;
        }
      }
    }
  }
}

void HeuristicMapper::Node::updateHeuristicCost(
    const Architecture&                      arch,
    const SingleQubitMultiplicity&           singleQubitGateMultiplicity,
    const TwoQubitMultiplicity&              twoQubitGateMultiplicity,
    const std::unordered_set<std::uint16_t>& consideredQubits) {
  costHeur = 0.;
  done     = true;

  // single qubit gate savings potential by moving them to different physical
  // qubits with higher fidelity
  double savingsPotential = 0.;
  if (considerFidelity) {
    for (std::uint16_t logQbit = 0U; logQbit < arch.getNqubits(); ++logQbit) {
      if (singleQubitGateMultiplicity.at(logQbit) == 0) {
        continue;
      }
      double       qbitSavings  = 0;
      const double currFidelity = arch.getSingleQubitFidelityCost(
          static_cast<std::uint16_t>(locations.at(logQbit)));
      for (std::uint16_t physQbit = 0U; physQbit < arch.getNqubits();
           ++physQbit) {
        if (arch.getSingleQubitFidelityCost(physQbit) >= currFidelity) {
          continue;
        }
        const double curSavings =
            singleQubitGateMultiplicity.at(logQbit) *
                (currFidelity - arch.getSingleQubitFidelityCost(physQbit)) -
            arch.fidelityDistance(
                static_cast<std::uint16_t>(locations.at(logQbit)), physQbit,
                consideredQubits.size());
        qbitSavings = std::max(qbitSavings, curSavings);
      }
      savingsPotential += qbitSavings;
    }
  }

  // iterating over all virtual qubit pairs, that share a gate on the
  // current layer
  for (const auto& [edge, multiplicity] : twoQubitGateMultiplicity) {
    const auto& [q1, q2] = edge;

    const auto& [straightMultiplicity, reverseMultiplicity] = multiplicity;

    const bool edgeDone =
        (arch.isEdgeConnected({static_cast<std::uint16_t>(locations.at(q1)),
                               static_cast<std::uint16_t>(locations.at(q2))}) ||
         arch.isEdgeConnected({static_cast<std::uint16_t>(locations.at(q2)),
                               static_cast<std::uint16_t>(locations.at(q1))}));
    // only if all qubit pairs are mapped next to each other the mapping
    // is complete
    if (!edgeDone) {
      done = false;
    }

    if (!considerFidelity) {
      costHeur += arch.distance(static_cast<std::uint16_t>(locations.at(q1)),
                        static_cast<std::uint16_t>(locations.at(q2)));
    } else {
      // find the optimal edge, to which to remap the given virtual qubit
      // pair and take the cost of moving it there via swaps plus the
      // fidelity cost  of executing all their shared gates on that edge
      // as the qubit pairs cost
      double swapCost = std::numeric_limits<double>::max();
      for (const auto& [q3, q4] : arch.getCouplingMap()) {
        swapCost = std::min(
            swapCost,
            straightMultiplicity * arch.getTwoQubitFidelityCost(q3, q4) +
                reverseMultiplicity * arch.getTwoQubitFidelityCost(q4, q3) +
                arch.fidelityDistance(
                    static_cast<std::uint16_t>(locations.at(q1)), q3,
                    consideredQubits.size()) +
                arch.fidelityDistance(
                    static_cast<std::uint16_t>(locations.at(q2)), q4,
                    consideredQubits.size()));
        swapCost = std::min(
            swapCost,
            straightMultiplicity * arch.getTwoQubitFidelityCost(q4, q3) +
                reverseMultiplicity * arch.getTwoQubitFidelityCost(q3, q4) +
                arch.fidelityDistance(
                    static_cast<std::uint16_t>(locations.at(q2)), q3,
                    consideredQubits.size()) +
                arch.fidelityDistance(
                    static_cast<std::uint16_t>(locations.at(q1)), q4,
                    consideredQubits.size()));
      }

      if (edgeDone) {
        const double currEdgeCost =
            (straightMultiplicity *
                 arch.getTwoQubitFidelityCost(
                     static_cast<std::uint16_t>(locations.at(q1)),
                     static_cast<std::uint16_t>(locations.at(q2))) +
             reverseMultiplicity *
                 arch.getTwoQubitFidelityCost(
                     static_cast<std::uint16_t>(locations.at(q2)),
                     static_cast<std::uint16_t>(locations.at(q1))));
        savingsPotential += (currEdgeCost - swapCost);
      } else {
        costHeur += swapCost;
      }
    }
  }
  
  if(!considerFidelity && admissibleHeuristic){
    auto n = consideredQubits.size();
    if (arch.bidirectional()) {
      costHeur -= ((n-1)*n/2 - sharedSwaps)*COST_BIDIRECTIONAL_SWAP;
    } else {
      costHeur -= ((n-1)*n/2 - sharedSwaps)*COST_UNIDIRECTIONAL_SWAP;
    }
  }
  
  costHeur -= savingsPotential;
}
