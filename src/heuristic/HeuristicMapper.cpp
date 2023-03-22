//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "heuristic/HeuristicMapper.hpp"

#include <chrono>

void HeuristicMapper::map(const Configuration& configuration) {
  results.config = configuration;
  auto& config   = results.config;
  if (config.layering == Layering::OddGates ||
      config.layering == Layering::QubitTriangle) {
    std::cerr << "Layering strategy " << toString(config.layering)
              << " not suitable for heuristic mapper!" << std::endl;
    return;
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

  std::size_t              gateidx = 0;
  std::vector<std::size_t> gatesToAdjust{};
  results.output.gates = 0U;
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const Node result = aStarMap(i);

    qubits    = result.qubits;
    locations = result.locations;

    if (config.verbose) {
      printLocations(std::clog);
      printQubits(std::clog);
    }

    // initial layer needs no swaps
    if (i != 0) {
      for (const auto& swaps : result.swaps) {
        for (const auto& swap : swaps) {
          if (swap.op == qc::SWAP) {
            if (config.verbose) {
              std::clog << "SWAP: " << swap.first << " <-> " << swap.second
                        << "\n";
            }
            if (architecture.getCouplingMap().find({swap.first, swap.second}) ==
                    architecture.getCouplingMap().end() &&
                architecture.getCouplingMap().find({swap.second, swap.first}) ==
                    architecture.getCouplingMap().end()) {
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
    for (const auto& gate : layers.at(i)) {
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
        if (architecture.getCouplingMap().find(cnot) ==
            architecture.getCouplingMap().end()) {
          const Edge reverse = {cnot.second, cnot.first};
          if (architecture.getCouplingMap().find(reverse) ==
              architecture.getCouplingMap().end()) {
            throw QMAPException(
                "Invalid CNOT: " + std::to_string(reverse.first) + "-" +
                std::to_string(reverse.second));
          }
          qcMapped.h(reverse.first);
          qcMapped.h(reverse.second);
          qcMapped.x(reverse.second,
                     qc::Control{static_cast<qc::Qubit>(reverse.first)});
          qcMapped.h(reverse.second);
          qcMapped.h(reverse.first);

          results.output.directionReverse++;
          gateidx += 5;
        } else {
          qcMapped.x(cnot.second,
                     qc::Control{static_cast<qc::Qubit>(cnot.first)});
          gateidx++;
        }
      }
    }
  }

  // infer output permutation from qubit locations
  qcMapped.outputPermutation.clear();
  std::size_t count = 0U;
  for (std::size_t i = 0U; i < architecture.getNqubits(); ++i) {
    if (qubits.at(i) != -1) {
      qcMapped.outputPermutation[static_cast<qc::Qubit>(i)] =
          static_cast<qc::Qubit>(qubits.at(i));
    } else {
      qcMapped.setLogicalQubitGarbage(
          static_cast<qc::Qubit>(qc.getNqubits() + count));
      ++count;
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
      if (gatesToAdjust.back() == gateidx) {
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

  postMappingOptimizations(config);
  countGates(qcMapped, results.output);
  finalizeMappedCircuit();

  const auto                          end  = std::chrono::steady_clock::now();
  const std::chrono::duration<double> diff = end - start;
  results.time                             = diff.count();
  results.timeout                          = false;
}

void HeuristicMapper::staticInitialMapping() {
  for (const auto& gate : layers.at(0U)) {
    if (gate.singleQubit()) {
      continue;
    }

    for (const auto& [q0, q1] : architecture.getCouplingMap()) {
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
  for (qc::Qubit i = 0U; i < architecture.getNqubits(); ++i) {
    if (qc.initialLayout.count(i) > 0 && locations.at(i) == DEFAULT_POSITION) {
      for (qc::Qubit j = 0U; j < architecture.getNqubits(); ++j) {
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

    std::uniform_int_distribution<> dis(0, architecture.getNqubits() - 1);

    for (std::size_t i = 0; i < config.teleportationQubits; i += 2) {
      Edge e{};
      do {
        auto it = std::begin(architecture.getCouplingMap());
        std::advance(it, dis(mt));
        e = *it;
      } while (qubits.at(e.first) != -1 || qubits.at(e.second) != -1);
      locations.at(qc.getNqubits() + i) = static_cast<std::int16_t>(e.first);
      locations.at(qc.getNqubits() + i + 1) =
          static_cast<std::int16_t>(e.second);
      qubits.at(e.first)  = static_cast<std::int16_t>(qc.getNqubits() + i);
      qubits.at(e.second) = static_cast<std::int16_t>(qc.getNqubits() + i + 1);
    }

    if (config.teleportationFake) {
      config.teleportationQubits = 0;
    }
  }

  switch (config.initialLayout) {
  case InitialLayout::Identity:
    for (qc::Qubit i = 0; i < architecture.getNqubits(); ++i) {
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
  case InitialLayout::None:
    // nothing to be done here
    break;

    // TODO: Design strategy that maps most used qubit to most connected qubits
    // on architecture
  }
}

void HeuristicMapper::mapUnmappedGates(
    const TwoQubitMultiplicity& twoQubitGateMultiplicity) {
  for (const auto& [logEdge, _] : twoQubitGateMultiplicity) {
    const auto& [q1, q2] = logEdge;

    auto q1Location = locations.at(q1);
    auto q2Location = locations.at(q2);

    if (q1Location == DEFAULT_POSITION && q2Location == DEFAULT_POSITION) {
      std::set<Edge> possibleEdges{};
      // gather all edges in the architecture for which both qubits are unmapped
      for (const auto& edge : architecture.getCouplingMap()) {
        if (qubits.at(edge.first) == DEFAULT_POSITION &&
            qubits.at(edge.second) == DEFAULT_POSITION) {
          possibleEdges.emplace(edge);
        }
      }
      std::pair<std::uint16_t, std::uint16_t> chosenEdge;

      if (possibleEdges.empty()) {
        // map to 2 qubits with minimal distance
        double bestScore = std::numeric_limits<int>::max();

        for (std::uint16_t i = 0; i < architecture.getNqubits(); i++) {
          for (std::uint16_t j = i + 1; j < architecture.getNqubits(); j++) {
            if (qubits.at(i) == DEFAULT_POSITION &&
                qubits.at(j) == DEFAULT_POSITION) {
              const double dist = architecture.distance(i, j);
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
  for (std::uint16_t i = 0; i < architecture.getNqubits(); ++i) {
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

HeuristicMapper::Node HeuristicMapper::aStarMap(size_t layer) {
  std::unordered_set<std::uint16_t> consideredQubits{};
  Node                              node{};
  TwoQubitMultiplicity              twoQubitGateMultiplicity{};

  for (const auto& gate : layers.at(layer)) {
    if (!gate.singleQubit()) {
      consideredQubits.emplace(gate.control);
      consideredQubits.emplace(gate.target);
      if (gate.control >= gate.target) {
        const auto edge =
            std::pair(gate.target, static_cast<std::uint16_t>(gate.control));
        if (twoQubitGateMultiplicity.find(edge) ==
            twoQubitGateMultiplicity.end()) {
          twoQubitGateMultiplicity[edge] = {0, 1};
        } else {
          twoQubitGateMultiplicity[edge].second++;
        }
      } else {
        const auto edge =
            std::pair(static_cast<std::uint16_t>(gate.control), gate.target);
        if (twoQubitGateMultiplicity.find(edge) ==
            twoQubitGateMultiplicity.end()) {
          twoQubitGateMultiplicity[edge] = {1, 0};
        } else {
          twoQubitGateMultiplicity[edge].first++;
        }
      }
    }
  }

  mapUnmappedGates(twoQubitGateMultiplicity);

  node.locations = locations;
  node.qubits    = qubits;
  node.recalculateFixedCost(architecture);
  node.updateHeuristicCost(architecture, twoQubitGateMultiplicity,
                           results.config.admissibleHeuristic);

  nodes.push(node);

  while (!nodes.top().done) {
    Node current = nodes.top();
    nodes.pop();
    expandNode(consideredQubits, current, layer, twoQubitGateMultiplicity);
  }

  Node result = nodes.top();
  nodes.pop();

  // clear nodes
  while (!nodes.empty()) {
    nodes.pop();
  }

  return result;
}

void HeuristicMapper::expandNode(
    const std::unordered_set<std::uint16_t>& consideredQubits, Node& node,
    std::size_t layer, const TwoQubitMultiplicity& twoQubitGateMultiplicity) {
  std::vector<std::vector<bool>> usedSwaps;
  usedSwaps.reserve(architecture.getNqubits());
  for (int p = 0; p < architecture.getNqubits(); ++p) {
    usedSwaps.emplace_back(architecture.getNqubits());
  }

  // set up new teleportation qubits
  std::set<Edge> perms = architecture.getCouplingMap();
  architecture.getCurrentTeleportations().clear();
  architecture.getTeleportationQubits().clear();
  for (std::size_t i = 0; i < results.config.teleportationQubits; i += 2) {
    architecture.getTeleportationQubits().emplace_back(
        node.locations.at(qc.getNqubits() + i),
        node.locations.at(qc.getNqubits() + i + 1));
    Edge e;
    for (auto const& g : architecture.getCouplingMap()) {
      if (g.first == node.locations.at(qc.getNqubits() + i) &&
          g.second != node.locations.at(qc.getNqubits() + i + 1)) {
        e.first  = g.second;
        e.second = static_cast<std::uint16_t>(
            node.locations.at(qc.getNqubits() + i + 1));
        architecture.getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.second == node.locations.at(qc.getNqubits() + i) &&
          g.first != node.locations.at(qc.getNqubits() + i + 1)) {
        e.first  = g.first;
        e.second = static_cast<std::uint16_t>(
            node.locations.at(qc.getNqubits() + i + 1));
        architecture.getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.first == node.locations.at(qc.getNqubits() + i + 1) &&
          g.second != node.locations.at(qc.getNqubits() + i)) {
        e.first = g.second;
        e.second =
            static_cast<std::uint16_t>(node.locations.at(qc.getNqubits() + i));
        architecture.getCurrentTeleportations().insert(e);
        perms.insert(e);
      }
      if (g.second == node.locations.at(qc.getNqubits() + i + 1) &&
          g.first != node.locations.at(qc.getNqubits() + i)) {
        e.first = g.first;
        e.second =
            static_cast<std::uint16_t>(node.locations.at(qc.getNqubits() + i));
        architecture.getCurrentTeleportations().insert(e);
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
          expandNodeAddOneSwap(edge, node, layer, twoQubitGateMultiplicity);
        } else if (!usedSwaps.at(static_cast<std::size_t>(q1))
                        .at(static_cast<std::size_t>(q2))) {
          usedSwaps.at(static_cast<std::size_t>(q1))
              .at(static_cast<std::size_t>(q2)) = true;
          usedSwaps.at(static_cast<std::size_t>(q2))
              .at(static_cast<std::size_t>(q1)) = true;
          expandNodeAddOneSwap(edge, node, layer, twoQubitGateMultiplicity);
        }
      }
    }
  }
}

void HeuristicMapper::expandNodeAddOneSwap(
    const Edge& swap, Node& node, const std::size_t layer,
    const TwoQubitMultiplicity& twoQubitGateMultiplicity) {
  const auto& config = results.config;

  Node newNode = Node(node.qubits, node.locations, node.swaps, node.costFixed);

  if (architecture.getCouplingMap().find(swap) !=
          architecture.getCouplingMap().end() ||
      architecture.getCouplingMap().find(Edge{swap.second, swap.first}) !=
          architecture.getCouplingMap().end()) {
    newNode.applySWAP(swap, architecture);
  } else {
    newNode.applyTeleportation(swap, architecture);
  }

  newNode.updateHeuristicCost(architecture, twoQubitGateMultiplicity,
                              results.config.admissibleHeuristic);

  // calculate heuristics for the cost of the following layers
  if (config.lookahead) {
    lookahead(getNextLayer(layer), newNode);
  }

  nodes.push(newNode);
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
        for (std::uint16_t j = 0; j < architecture.getNqubits(); ++j) {
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
        for (std::uint16_t j = 0; j < architecture.getNqubits(); ++j) {
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
        auto cost = architecture.distance(
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

void HeuristicMapper::Node::applySWAP(const Edge& swap, Architecture& arch) {
  ++nswaps;
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

  if (arch.getCouplingMap().find(swap) != arch.getCouplingMap().end() ||
      arch.getCouplingMap().find(Edge{swap.second, swap.first}) !=
          arch.getCouplingMap().end()) {
    swaps.back().emplace_back(swap.first, swap.second, qc::SWAP);
  } else {
    throw QMAPException("Something wrong in applySWAP.");
  }

  if (arch.bidirectional()) {
    costFixed += COST_BIDIRECTIONAL_SWAP;
  } else {
    costFixed += COST_UNIDIRECTIONAL_SWAP;
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
  if (arch.getCouplingMap().find({swap.first, middleAnc}) !=
          arch.getCouplingMap().end() ||
      arch.getCouplingMap().find({middleAnc, swap.first}) !=
          arch.getCouplingMap().end()) {
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

void HeuristicMapper::Node::recalculateFixedCost(const Architecture& arch) {
  costFixed = 0;
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

void HeuristicMapper::Node::updateHeuristicCost(
    const Architecture&         arch,
    const TwoQubitMultiplicity& twoQubitGateMultiplicity,
    const bool                  admissibleHeuristic) {
  costHeur = 0.;
  done     = true;

  // iterating over all virtual qubit pairs, that share a gate on the
  // current layer
  for (const auto& [edge, multiplicity] : twoQubitGateMultiplicity) {
    const auto& [q1, q2] = edge;

    const auto& [straightMultiplicity, reverseMultiplicity] = multiplicity;

    // only if all qubit pairs are mapped next to each other the mapping
    // is complete
    if (arch.getCouplingMap().find(
            {static_cast<std::uint16_t>(locations.at(q1)),
             static_cast<std::uint16_t>(locations.at(q2))}) ==
            arch.getCouplingMap().end() &&
        arch.getCouplingMap().find(
            {static_cast<std::uint16_t>(locations.at(q2)),
             static_cast<std::uint16_t>(locations.at(q1))}) ==
            arch.getCouplingMap().end()) {
      done = false;
    }

    const double swapCostStraight =
        arch.distance(static_cast<std::uint16_t>(locations.at(q1)),
                      static_cast<std::uint16_t>(locations.at(q2)));
    const double swapCostReverse =
        arch.distance(static_cast<std::uint16_t>(locations.at(q2)),
                      static_cast<std::uint16_t>(locations.at(q1)));

    if (admissibleHeuristic) {
      if (straightMultiplicity > 0) {
        costHeur = std::max(costHeur, swapCostStraight);
      }
      if (reverseMultiplicity > 0) {
        costHeur = std::max(costHeur, swapCostReverse);
      }
    } else {
      costHeur += swapCostStraight * straightMultiplicity +
                  swapCostReverse * reverseMultiplicity;
    }
  }
}
