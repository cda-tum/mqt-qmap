#include "na/nasp/SolverFactory.hpp"

#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "na/Architecture.hpp"
#include "na/NADefinitions.hpp"
#include "na/nasp/Solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace na {
auto SolverFactory::create(const Architecture& arch) -> NASolver {
  const auto interactionZone =
      *arch.getPropertiesOfOperation({qc::Z, 1}).zones.cbegin();
  const auto maxX =
      static_cast<uint16_t>(arch.getNColsInZone(interactionZone) - 1);
  const auto maxEntanglingY =
      static_cast<uint16_t>(arch.getNrowsInZone(interactionZone) - 1);
  const auto maxC =
      static_cast<uint16_t>(arch.getPropertiesOfShuttlingUnit(0).cols);
  const auto maxR =
      static_cast<uint16_t>(arch.getPropertiesOfShuttlingUnit(0).rows);
  const auto storageZone = arch.getInitialZones().front();
  const auto maxY =
      static_cast<uint16_t>(maxEntanglingY + arch.getNrowsInZone(storageZone));
  // the atoms are located, e.g., in the following manner:
  //   0 <-- SLM               0 <-- SLM
  // o o o <-- AOD    OR    o o o o <-- AOD
  // o o o <-- AOD          o o o o <-- AOD
  // The max number of columns and rows that can be stacked at one interaction
  // site is derived from the fact that the SLM atom must still be in the
  // Rydberg interaction radius of the lower left AOD atom. For that we first
  // assume the same value for the vertical and horizontal stacking factor
  // and solve the resulting quadratic equation. For the solution we take the
  // flooring of the positive solution. If the equation allows the horizontal
  // stacking factor to be one larger, we take that value.
  const auto minAtomDistance = arch.getMinAtomDistance();
  const auto interactionRadius = arch.getInteractionRadius();
  const auto minAtomDistanceSquared = std::pow(minAtomDistance, 2);
  const auto interactionRadiusSquared = std::pow(interactionRadius, 2);
  const auto maxVDist = static_cast<uint16_t>(
      std::min(std::floor(.2 + std::sqrt(.8 * interactionRadiusSquared /
                                             minAtomDistanceSquared -
                                         .16)),
               std::floor(0.7071067811865475244 *
                              static_cast<double>(interactionRadius) /
                              static_cast<double>(minAtomDistance) +
                          1.)));
  const auto maxHDist =
      std::pow(maxVDist * minAtomDistance, 2) +
                      std::pow(maxVDist / 2 * minAtomDistance, 2) <=
                  interactionRadiusSquared &&
              std::pow((maxVDist - 1) * minAtomDistance, 2) +
                      std::pow(maxVDist * minAtomDistance, 2) <=
                  interactionRadiusSquared
          ? static_cast<uint16_t>(maxVDist + 1)
          : maxVDist;
  const auto noInteractionRadius = arch.getNoInteractionRadius();
  const auto firstSite = arch.getSitesInZone(interactionZone)[0];
  const auto siteRight =
      arch.getNearestSiteRight(arch.getPositionOfSite(firstSite), true, true);
  const auto siteBelow =
      arch.getNearestSiteDown(arch.getPositionOfSite(firstSite), true, true);
  if (!siteRight || !siteBelow) {
    throw std::invalid_argument(
        "Unexpected architecture: There is no site to the right or below the "
        "first site in the interaction zone.");
  }
  const auto maxHOffset = static_cast<uint16_t>(
      (arch.getPositionOfSite(*siteRight) - arch.getPositionOfSite(firstSite))
          .length() -
      noInteractionRadius / 2 / minAtomDistance);
  const auto maxVOffset = static_cast<uint16_t>(
      (arch.getPositionOfSite(*siteBelow) - arch.getPositionOfSite(firstSite))
          .length() -
      noInteractionRadius / 2 / minAtomDistance);
  NASolver solver{};
  solver.init(maxX, maxY, maxC, maxR, maxHOffset, maxVOffset, maxHDist,
              maxVDist, 0, maxEntanglingY);
  return solver;
}

auto SolverFactory::getOpsForSolver(const qc::QuantumComputation& circ,
                                    const FullOpType opType, const bool quiet)
    -> std::vector<std::pair<unsigned int, unsigned int>> {
  std::vector<std::pair<unsigned int, unsigned int>> ops;
  ops.reserve(circ.size());
  for (const auto& op : circ) {
    if (op->getType() == opType.type &&
        op->getNcontrols() == opType.nControls) {
      const auto& operands = op->getUsedQubits();
      if (operands.size() != 2) {
        std::stringstream ss;
        ss << "Operation " << op->getName() << " does not have two operands.";
        throw std::invalid_argument(ss.str());
      }
      ops.emplace_back(*operands.cbegin(), *operands.rbegin());
    } else if (!quiet) {
      std::stringstream ss;
      ss << "Operation " << op->getName() << " is not of type " << opType.type
         << " or does not have " << opType.nControls << " controls.";
      throw std::invalid_argument(ss.str());
    }
  }
  return ops;
}

} // namespace na
