#include "na/nasp/CodeGenerator.hpp"

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "datastructures/Layer.hpp"
#include "na/NAComputation.hpp"
#include "na/NADefinitions.hpp"
#include "na/operations/NAGlobalOperation.hpp"
#include "na/operations/NALocalOperation.hpp"
#include "na/operations/NAShuttlingOperation.hpp"
#include "na/nasp/Solver.hpp"
#include "ir/operations/OpType.hpp"

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
    const std::int32_t x, const std::int32_t y, const std::int32_t h,
    const std::int32_t v, const std::int32_t maxHOffset,
    const std::int32_t maxVOffset, const std::int32_t minEntanglingY,
    const std::int32_t maxEntanglingY) -> Point {
  constexpr auto minAtomDist = 1;
  constexpr auto noInteractionRadius = 10;
  constexpr auto zoneDist = 24; // incl., 2 * maxHOffset * minAtomDist
  const auto dx = static_cast<std::int64_t>(noInteractionRadius) +
                  2LL * maxHOffset * minAtomDist;
  const auto dy = static_cast<std::int64_t>(noInteractionRadius) +
                  2LL * maxVOffset * minAtomDist;
  if (minEntanglingY == 0) {
    // no top storage zone
    if (y <= maxEntanglingY) {
      return {static_cast<std::int64_t>(x) * dx +
                  static_cast<std::int64_t>(h) * minAtomDist,
              static_cast<std::int64_t>(y) * dy +
                  static_cast<std::int64_t>(v) * minAtomDist};
    }
    return {static_cast<std::int64_t>(x) * dx +
                static_cast<std::int64_t>(h) * minAtomDist,
            zoneDist + static_cast<std::int64_t>(y - 1) * dy +
                static_cast<std::int64_t>(v) * minAtomDist};
  }
  // top storage zone
  if (y < minEntanglingY) {
    return {static_cast<std::int64_t>(x) * dx +
                static_cast<std::int64_t>(h) * minAtomDist,
            static_cast<std::int64_t>(y) * dy +
                static_cast<std::int64_t>(v) * minAtomDist};
  }
  if (y <= maxEntanglingY) {
    return {static_cast<std::int64_t>(x) * dx +
                static_cast<std::int64_t>(h) * minAtomDist,
            zoneDist + static_cast<std::int64_t>(y - 1) * dy +
                static_cast<std::int64_t>(v) * minAtomDist};
  }
  return {static_cast<std::int64_t>(x) * dx +
              static_cast<std::int64_t>(h) * minAtomDist,
          2LL * zoneDist + static_cast<std::int64_t>(y - 2) * dy +
              static_cast<std::int64_t>(v) * minAtomDist};
}

auto CodeGenerator::generate(
    const QuantumComputation& input, const NASolver::Result& result,
    const std::uint16_t maxHOffset, const std::uint16_t maxVOffset,
    const std::uint16_t minEntanglingY,
    const std::uint16_t maxEntanglingY) -> NAComputation {
  Layer const layer(input);
  NAComputation code;
  std::vector<std::shared_ptr<Point>> oldPositions;
  std::vector<bool> wasAOD;
  oldPositions.reserve(result.front().numQubits());
  wasAOD.reserve(result.front().numQubits());
  // initialize atoms in SLM and load required ones into AOD
  {
    std::unordered_map<std::size_t, std::set<std::size_t>> hAODLines{};
    std::unordered_map<std::size_t, std::set<std::size_t>> vAODLines{};
    for (const auto& q : result.front().getQubits()) {
      if (q.isAOD()) {
        hAODLines[static_cast<std::size_t>(q.getY())].emplace(q.getR());
        vAODLines[static_cast<std::size_t>(q.getX())].emplace(q.getC());
      }
    }
    std::vector<std::shared_ptr<Point>> loadPositions;
    for (const auto& q : result.front().getQubits()) {
      auto pos = std::make_shared<Point>(
          coordFromDiscrete(q.getX(), q.getY(), q.getH(), q.getV(), maxHOffset,
                            maxVOffset, minEntanglingY, maxEntanglingY));
      wasAOD.emplace_back(q.isAOD());
      if (q.isAOD()) {
        loadPositions.emplace_back(pos);
      }
      oldPositions.emplace_back(pos);
      code.emplaceInitialPosition(pos);
    }
    // TODO properly load atoms into AOD with offset
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
  if (affectedQubits.size() != input.getNqubits() ||
      ops.size() != affectedQubits.size()) {
    throw std::invalid_argument("Not all atoms are initialized to plus state.");
  }
  // initialize atoms in to |+> state starting in |0> state
  code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{RY, 0},
                                                       std::vector{PI_2}));
  std::for_each(ops.cbegin(), ops.cend(), [](const auto& v) { v->execute(); });
  // Reference to the executable set of the input circuit
  const auto& executableSet = layer.getExecutableSet();
  if (result.front().isRydberg()) {
    code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{Z, 1}));
    // find and execute corresponding gates in input circuit
    for (const auto& g : result.front().getGates()) {
      const auto& it = std::find_if(
          executableSet.begin(), executableSet.end(), [g](const auto& v) {
            if (v->getOperation()->getType() == Z &&
                v->getOperation()->getNcontrols() == 1) {
              const auto& usedQubits = v->getOperation()->getUsedQubits();
              const auto [first, second] = g.getQubits();
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
  for (std::uint16_t t = 1; t < static_cast<std::uint16_t>(result.numStages());
       ++t) {
    std::vector<std::shared_ptr<Point>> newPositions;
    newPositions.reserve(oldPositions.size());
    std::vector<std::shared_ptr<Point>> startPositions;
    std::vector<std::shared_ptr<Point>> endPositions;
    std::vector<std::shared_ptr<Point>> loadStartPositions;
    std::vector<std::shared_ptr<Point>> loadEndPositions;
    std::vector<std::shared_ptr<Point>> storeStartPositions;
    std::vector<std::shared_ptr<Point>> storeEndPositions;
    for (std::uint16_t i = 0;
         i < static_cast<std::uint16_t>(result.getStage(t).numQubits()); ++i) {
      const auto& q = result.getStage(t).getQubit(i);
      auto pos = std::make_shared<Point>(
          coordFromDiscrete(q.getX(), q.getY(), q.getH(), q.getV(), maxHOffset,
                            maxVOffset, minEntanglingY, maxEntanglingY));
      if (wasAOD[i] && q.isAOD()) {
        startPositions.emplace_back(oldPositions[i]);
        endPositions.emplace_back(pos);
      } else if (wasAOD[i] && !q.isAOD()) {
        storeStartPositions.emplace_back(oldPositions[i]);
        storeEndPositions.emplace_back(pos);
      } else if (!wasAOD[i] && q.isAOD()) {
        loadStartPositions.emplace_back(oldPositions[i]);
        loadEndPositions.emplace_back(loadStartPositions.back());
        startPositions.emplace_back(loadEndPositions.back());
        endPositions.emplace_back(pos);
      }
      newPositions.emplace_back(pos);
      wasAOD[i] = q.isAOD();
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
    if (result.getStage(t).isRydberg()) {
      code.emplaceBack(std::make_unique<NAGlobalOperation>(FullOpType{Z, 1}));
    }
    oldPositions = std::move(newPositions);
    // find and execute corresponding gates in input circuit
    for (const auto& g : result.getStage(t).getGates()) {
      const auto& it = std::find_if(
          executableSet.begin(), executableSet.end(), [g](const auto& v) {
            if (v->getOperation()->getType() == Z &&
                v->getOperation()->getNcontrols() == 1) {
              const auto& usedQubits = v->getOperation()->getUsedQubits();
              const auto [first, second] = g.getQubits();
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
