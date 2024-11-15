#include "na/nasp/CodeGenerator.hpp"

#include "Definitions.hpp"
#include "circuit_optimizer/CircuitOptimizer.hpp"
#include "datastructures/Layer.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "na/NAComputation.hpp"
#include "na/NADefinitions.hpp"
#include "na/nasp/Solver.hpp"
#include "na/operations/NAGlobalOperation.hpp"
#include "na/operations/NALocalOperation.hpp"
#include "na/operations/NAShuttlingOperation.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {
using namespace qc;

auto CodeGenerator::coordFromDiscrete(
    const NASolver::Result::Qubit q, const int64_t maxHOffset,
    const int64_t maxVOffset, const int64_t minEntanglingY,
    const int64_t maxEntanglingY, const int64_t minAtomDist,
    const int64_t noInteractionRadius, const int64_t zoneDist) -> Point {
  const auto dx = noInteractionRadius + (2LL * maxHOffset * minAtomDist);
  const auto dy = noInteractionRadius + (2LL * maxVOffset * minAtomDist);
  const auto x = q.x;
  const auto y = q.y;
  const auto h = q.h;
  const auto v = q.v;
  if (minEntanglingY == 0) {
    // no top storage zone
    if (y <= maxEntanglingY) {
      return {(x * dx) + (h * minAtomDist), (y * dy) + (v * minAtomDist)};
    }
    return {(x * dx) + (h * minAtomDist),
            zoneDist + ((y - 1) * dy) + (v * minAtomDist)};
  }
  // top storage zone
  if (y < minEntanglingY) {
    return {(x * dx) + (h * minAtomDist), (y * dy) + (v * minAtomDist)};
  }
  if (y <= maxEntanglingY) {
    return {(x * dx) + (h * minAtomDist),
            zoneDist + ((y - 1) * dy) + (v * minAtomDist)};
  }
  return {(x * dx) + (h * minAtomDist),
          (2LL * zoneDist) + ((y - 2) * dy) + (v * minAtomDist)};
}

auto CodeGenerator::generate(
    const QuantumComputation& input, const NASolver::Result& result,
    const uint16_t maxHOffset, const uint16_t maxVOffset,
    const uint16_t minEntanglingY, const uint16_t maxEntanglingY,
    const uint16_t minAtomDist, const uint16_t noInteractionRadius,
    const uint16_t zoneDist) -> NAComputation {
  auto flattened = input;
  CircuitOptimizer::flattenOperations(flattened);
  const Layer layer(flattened);
  NAComputation code;
  std::vector<std::shared_ptr<Point>> oldPositions;
  std::vector<bool> wasAOD;
  oldPositions.reserve(result.stages.front().qubits.size());
  wasAOD.reserve(result.stages.front().qubits.size());
  // initialize atoms in SLM and load required ones into AOD
  {
    std::unordered_map<std::size_t, std::set<std::size_t>> hAODLines{};
    std::unordered_map<std::size_t, std::set<std::size_t>> vAODLines{};
    for (const auto& q : result.stages.front().qubits) {
      if (q.a) {
        hAODLines[static_cast<std::size_t>(q.y)].emplace(q.r);
        vAODLines[static_cast<std::size_t>(q.x)].emplace(q.c);
      }
    }
    std::vector<std::shared_ptr<Point>> loadPositions;
    for (const auto& q : result.stages.front().qubits) {
      auto pos = std::make_shared<Point>(coordFromDiscrete(
          q, maxHOffset, maxVOffset, minEntanglingY, maxEntanglingY,
          minAtomDist, noInteractionRadius, zoneDist));
      wasAOD.emplace_back(q.a);
      if (q.a) {
        loadPositions.emplace_back(pos);
      }
      oldPositions.emplace_back(pos);
      code.emplaceInitialPosition(pos);
    }
    if (!loadPositions.empty()) {
      code.emplaceBack<NAShuttlingOperation>(LOAD, loadPositions,
                                             loadPositions);
    }
  }
  const auto ops = layer.getExecutablesOfType(H, 0);
  std::unordered_set<Qubit> affectedQubits;
  std::transform(
      ops.cbegin(), ops.cend(),
      std::inserter(affectedQubits, affectedQubits.end()),
      [](const auto& v) { return v->getOperation()->getTargets().front(); });
  if (affectedQubits.size() != flattened.getNqubits() ||
      ops.size() != affectedQubits.size()) {
    throw std::invalid_argument("Not all atoms are initialized to plus state.");
  }
  // initialize atoms in to |+> state starting in |0> state
  code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{RY, 0},
                                                       std::vector{PI_2}));
  std::for_each(ops.cbegin(), ops.cend(), [](const auto& v) { v->execute(); });
  // Reference to the executable set of the input circuit
  const auto& executableSet = layer.getExecutableSet();
  if (result.stages.front().rydberg) {
    code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{Z, 1}));
    // find and execute corresponding gates in input circuit
    for (const auto& g : result.stages.front().gates) {
      const auto& it = std::find_if(
          executableSet.begin(), executableSet.end(), [g](const auto& v) {
            if (v->getOperation()->getType() == Z &&
                v->getOperation()->getNcontrols() == 1) {
              const auto& usedQubits = v->getOperation()->getUsedQubits();
              const auto [first, second] = g.qubits;
              return std::set{first, second} == usedQubits;
            }
            return false;
          });
      if (it == executableSet.end()) {
        throw std::invalid_argument(
            "Gate in input circuit has no correspondence in solution.");
      }
      (*it)->execute();
    }
  }
  for (uint16_t t = 1; t < static_cast<uint16_t>(result.stages.size()); ++t) {
    std::vector<std::shared_ptr<Point>> newPositions;
    newPositions.reserve(oldPositions.size());
    std::vector<std::shared_ptr<Point>> startPositions;
    std::vector<std::shared_ptr<Point>> endPositions;
    std::vector<std::shared_ptr<Point>> loadStartPositions;
    std::vector<std::shared_ptr<Point>> loadEndPositions;
    std::vector<std::shared_ptr<Point>> storeStartPositions;
    std::vector<std::shared_ptr<Point>> storeEndPositions;
    for (uint16_t i = 0;
         static_cast<size_t>(i) < result.stages.at(t).qubits.size(); ++i) {
      const auto& q = result.stages.at(t).qubits.at(i);
      auto pos = std::make_shared<Point>(coordFromDiscrete(
          q, maxHOffset, maxVOffset, minEntanglingY, maxEntanglingY,
          minAtomDist, noInteractionRadius, zoneDist));
      if (wasAOD[i] && q.a) {
        startPositions.emplace_back(oldPositions[i]);
        endPositions.emplace_back(pos);
      } else if (wasAOD[i] && !q.a) {
        storeStartPositions.emplace_back(oldPositions[i]);
        storeEndPositions.emplace_back(pos);
      } else if (!wasAOD[i] && q.a) {
        loadStartPositions.emplace_back(oldPositions[i]);
        loadEndPositions.emplace_back(loadStartPositions.back());
        startPositions.emplace_back(loadEndPositions.back());
        endPositions.emplace_back(pos);
      }
      newPositions.emplace_back(pos);
      wasAOD[i] = q.a;
    }
    if (!storeEndPositions.empty()) {
      code.emplaceBack<NAShuttlingOperation>(STORE, storeStartPositions,
                                             storeEndPositions);
    }
    if (!loadEndPositions.empty()) {
      code.emplaceBack<NAShuttlingOperation>(LOAD, loadStartPositions,
                                             loadEndPositions);
    }
    if (!startPositions.empty()) {
      code.emplaceBack<NAShuttlingOperation>(MOVE, startPositions,
                                             endPositions);
    }
    if (result.stages.at(t).rydberg) {
      code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{Z, 1}));
    }
    oldPositions = std::move(newPositions);
    // find and execute corresponding gates in input circuit
    for (const auto& g : result.stages.at(t).gates) {
      const auto& it = std::find_if(
          executableSet.begin(), executableSet.end(), [g](const auto& v) {
            if (v->getOperation()->getType() == Z &&
                v->getOperation()->getNcontrols() == 1) {
              const auto& usedQubits = v->getOperation()->getUsedQubits();
              const auto [first, second] = g.qubits;
              return std::set{first, second} == usedQubits;
            }
            return false;
          });
      if (it == executableSet.end()) {
        throw std::invalid_argument(
            "Gate in input circuit has no correspondence in solution.");
      }
      (*it)->execute();
    }
  }
  if (!executableSet.empty()) {
    code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{RY, 0},
                                                         std::vector{-PI_4}));
    while (!executableSet.empty()) {
      const auto& v = (*executableSet.cbegin());
      if (v->getOperation()->getType() != H) {
        throw std::invalid_argument(
            "Not all non CZ-gates in input circuit are executed.");
      }
      const auto q = v->getOperation()->getTargets().front();
      const auto& pos = oldPositions[q];
      code.emplaceBack(std::make_unique<NALocalOperation>(
          FullOpType{RZ, 0}, std::vector{PI}, pos));
      v->execute();
    }
    code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{RY, 0},
                                                         std::vector{PI_4}));
  }
  return code;
}
} // namespace na
