#include "na/zoned/Architecture.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <memory>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace na::zoned {
AOD::AOD(nlohmann::json aodSpec) {
  if (aodSpec.contains("id")) {
    if (aodSpec["id"].is_number()) {
      id = aodSpec["id"];
    } else {
      throw std::invalid_argument(
          "AOD id must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument("AOD id is missed in architecture spec");
  }
  if (aodSpec.contains("site_separation")) {
    if (aodSpec["site_separation"].is_number()) {
      siteSeparation = aodSpec["site_separation"];
    } else {
      throw std::invalid_argument(
          "AOD site separation must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "AOD site separation is missed in architecture spec");
  }
  if (aodSpec.contains("r")) {
    if (aodSpec["r"].is_number()) {
      nRows = aodSpec["r"];
    } else {
      throw std::invalid_argument(
          "AOD row number must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "AOD row number is missed in architecture spec");
  }
  if (aodSpec.contains("c")) {
    if (aodSpec["c"].is_number()) {
      nCols = aodSpec["c"];
    } else {
      throw std::invalid_argument(
          "AOD column number must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "AOD column number is missed in architecture spec");
  }
}

SLM::SLM(nlohmann::json slmSpec) {
  if (slmSpec.contains("id")) {
    if (slmSpec["id"].is_number()) {
      id = slmSpec["id"];
    } else {
      throw std::invalid_argument(
          "SLM id must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument("SLM id is missed in architecture spec");
  }
  if (slmSpec.contains("site_separation")) {
    if (slmSpec["site_separation"].is_array() &&
        slmSpec["site_separation"].size() == 2 &&
        slmSpec["site_separation"][0].is_number() &&
        slmSpec["site_separation"][1].is_number()) {
      siteSeparation = slmSpec["site_separation"];
    } else {
      throw std::invalid_argument(
          "SLM site separation must be a 2D number array in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "SLM site separation is missed in architecture spec");
  }
  if (slmSpec.contains("r")) {
    if (slmSpec["r"].is_number()) {
      nRows = slmSpec["r"];
    } else {
      throw std::invalid_argument(
          "SLM row number must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "SLM row number is missed in architecture spec");
  }
  if (slmSpec.contains("c")) {
    if (slmSpec["c"].is_number()) {
      nCols = slmSpec["c"];
    } else {
      throw std::invalid_argument(
          "SLM column number must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "SLM column number is missed in architecture spec");
  }
  if (slmSpec.contains("location")) {
    if (slmSpec["location"].is_array() && slmSpec["location"].size() == 2 &&
        slmSpec["location"][0].is_number() &&
        slmSpec["location"][1].is_number()) {
      location = slmSpec["location"];
    } else {
      throw std::invalid_argument(
          "SLM location must be a 2D number array in architecture spec");
    }
  } else {
    throw std::invalid_argument("SLM location is missed in architecture spec");
  }
}

auto SLM::operator==(const SLM& other) const -> bool {
  if (&other == this) {
    return true;
  }
  if (other.location != location) {
    return false;
  }
  if (other.nRows != nRows || other.nCols != nCols) {
    return false;
  }
  if (other.siteSeparation != siteSeparation) {
    return false;
  }
  if (other.entanglementZone_ != entanglementZone_) {
    return false;
  }
  return true;
}

Architecture::Architecture(nlohmann::json json) {
  // check if the name exists and is a string, otherwise throw an error
  // JSON Example:
  // "name": "My Super Cool Architecture"
  if (json.contains("name")) {
    if (json["name"].is_string()) {
      name = json["name"];
    } else {
      throw std::invalid_argument("Architecture name must be a string");
    }
  } else {
    throw std::invalid_argument("Architecture name is missed in architecture "
                                "spec");
  }
  // check if the operation's duration exists, otherwise print a warning
  // throw an error if the specification is invalid
  // JSON Example:
  // "operation_duration": {
  //   "rydberg": 0.36,
  //   "1q_gate": 52,
  //   "atom_transfer": 15
  // }
  if (json.contains("operation_duration")) {
    if (json["operation_duration"].is_object()) {
      operationDurations = OperationDurations{};
      if (json["operation_duration"].contains("rydberg")) {
        if (json["operation_duration"]["rydberg"].is_number()) {
          operationDurations->timeRydberg =
              json["operation_duration"]["rydberg"];
        } else {
          throw std::invalid_argument(
              "Rydberg duration must be a number in architecture spec");
        }
      } else {
        throw std::invalid_argument(
            "Operation duration must contain rydberg duration");
      }
      if (json["operation_duration"].contains("atom_transfer")) {
        if (json["operation_duration"]["atom_transfer"].is_number()) {
          operationDurations->timeAtomTransfer =
              json["operation_duration"]["atom_transfer"];
        } else {
          throw std::invalid_argument(
              "Atom transfer duration must be a number in architecture spec");
        }
      } else {
        throw std::invalid_argument(
            "Operation duration must contain atom transfer duration");
      }
      if (json["operation_duration"].contains("1q_gate")) {
        if (json["operation_duration"]["1q_gate"].is_number()) {
          operationDurations->timeSingleQubitGate =
              json["operation_duration"]["1q_gate"];
        } else {
          throw std::invalid_argument(
              "One qubit gate duration must be a number in architecture spec");
        }
      } else {
        throw std::invalid_argument(
            "Operation duration must contain 1q_gate duration");
      }
    } else {
      throw std::invalid_argument(
          "Operation duration must be a dict in architecture spec");
    }
  } else {
    std::cout << "\033[1;35m[WARN]\033[0m Operation's duration is missed in "
                 "architecture spec. Using default values.\n";
  }
  // check if the operation's fidelity exists, otherwise print a warning
  // throw an error if the specification is invalid
  // "operation_fidelity": {
  //   "two_qubit_gate": 0.995,
  //   "single_qubit_gate": 0.9997,
  //   "atom_transfer": 0.999
  // }
  if (json.contains("operation_fidelity")) {
    if (json["operation_fidelity"].is_object()) {
      operationFidelities = OperationFidelities{};
      if (json["operation_fidelity"].contains("two_qubit_gate")) {
        if (json["operation_fidelity"]["two_qubit_gate"].is_number()) {
          operationFidelities->fidelityTwoQubitGate =
              json["operation_fidelity"]["two_qubit_gate"];
        } else {
          throw std::invalid_argument("Two qubit gate fidelity must be a float "
                                      "in architecture spec");
        }
      } else {
        throw std::invalid_argument("Operation fidelity must contain two qubit "
                                    "gate fidelity");
      }
      if (json["operation_fidelity"].contains("atom_transfer")) {
        if (json["operation_fidelity"]["atom_transfer"].is_number()) {
          operationFidelities->fidelityAtomTransfer =
              json["operation_fidelity"]["atom_transfer"];
        } else {
          throw std::invalid_argument("Atom transfer fidelity must be a float "
                                      "in architecture spec");
        }
      } else {
        throw std::invalid_argument("Operation fidelity must contain atom "
                                    "transfer fidelity");
      }
      if (json["operation_fidelity"].contains("single_qubit_gate")) {
        if (json["operation_fidelity"]["single_qubit_gate"].is_number()) {
          operationFidelities->fidelitySingleQubitGate =
              json["operation_fidelity"]["single_qubit_gate"];
        } else {
          throw std::invalid_argument("One qubit gate fidelity must be a float "
                                      "in architecture spec");
        }
      } else {
        throw std::invalid_argument("Operation fidelity must contain one qubit "
                                    "gate fidelity");
      }
    } else {
      throw std::invalid_argument(
          "Operation fidelities must be a dict in architecture spec");
    }
  } else {
    std::cout << "\033[1;35m[WARN]\033[0m Operation's fidelity is missed in "
                 "architecture spec. Using default values.\n";
  }
  // check if the qubit's T1 time exists, otherwise print a warning
  // throw an error if the time is not a number
  // JSON Example:
  // "qubit_spec": {
  //   "T": 1.5e6
  // }
  if (json.contains("qubit_spec")) {
    if (json["qubit_spec"].is_object()) {
      if (json["qubit_spec"].contains("T")) {
        if (json["qubit_spec"]["T"].is_number()) {
          qubitT1 = json["qubit_spec"]["T"];
        } else {
          throw std::invalid_argument("The qubit's T1 time must be a number in "
                                      "architecture spec");
        }
      } else {
        throw std::invalid_argument("The qubit spec must contain T1 time");
      }
    } else {
      throw std::invalid_argument(
          "The qubit spec must be a dict in architecture spec");
    }
  } else {
    std::cout << "\033[1;35m[WARN]\033[0m The qubit spec is missed in "
                 "architecture spec. Using default values.\n";
  }
  // check if the architecture range exists and is valid, otherwise throw an
  // error
  // JSON Example:
  // "arch_range": [
  //   [
  //     0,
  //     0
  //   ],
  //   [
  //     60,
  //     110
  //   ]
  //
  if (json.contains("arch_range")) {
    if (json["arch_range"].is_array() && json["arch_range"].size() == 2 &&
        json["arch_range"][0].is_array() && json["arch_range"][0].size() == 2 &&
        json["arch_range"][1].is_array() && json["arch_range"][1].size() == 2 &&
        json["arch_range"][0][0].is_number() &&
        json["arch_range"][0][1].is_number() &&
        json["arch_range"][1][0].is_number() &&
        json["arch_range"][1][1].is_number()) {
      archRangeMinX = json["arch_range"][0][0];
      archRangeMinY = json["arch_range"][0][1];
      archRangeMaxX = json["arch_range"][1][0];
      archRangeMaxY = json["arch_range"][1][1];
    } else {
      throw std::invalid_argument(
          "Architecture range must be a 2x2 number array in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "Architecture range is missed in architecture spec");
  }
  // check if the rydberg range exists and is valid, otherwise throw an error
  // JSON Example:
  // "rydberg_range": [
  //   [
  //     [
  //       0,
  //       57
  //     ],
  //     [
  //       65,
  //       105
  //     ]
  //   ]
  // ]
  if (json.contains("rydberg_range")) {
    if (json["rydberg_range"].is_array() && !json["rydberg_range"].empty()) {
      for (const auto& rydbergRange : json["rydberg_range"]) {
        if (rydbergRange.is_array() && rydbergRange.size() == 2 &&
            rydbergRange[0].is_array() && rydbergRange[0].size() == 2 &&
            rydbergRange[1].is_array() && rydbergRange[1].size() == 2 &&
            rydbergRange[0][0].is_number() && rydbergRange[0][1].is_number() &&
            rydbergRange[1][0].is_number() && rydbergRange[1][1].is_number()) {
          rydbergRangeMinX.emplace_back(rydbergRange[0][0]);
          rydbergRangeMinY.emplace_back(rydbergRange[0][1]);
          rydbergRangeMaxX.emplace_back(rydbergRange[1][0]);
          rydbergRangeMaxY.emplace_back(rydbergRange[1][1]);
        } else {
          throw std::invalid_argument("Rydberg range must be a Nx2x2 number "
                                      "array in architecture spec, N > 1");
        }
      }
    } else {
      throw std::invalid_argument("Rydberg range must be a Nx2x2 number array "
                                  "in architecture spec, N > 1");
    }
  } else {
    throw std::invalid_argument("Rydberg range is missed in architecture spec");
  }
  // check if the storage zones exist and are valid, otherwise throw an error
  // JSON Example:
  // "storage_zones": [
  //   {
  //     "slms": [
  //       {
  //         "id": 0,
  //         "site_separation": [
  //           3,
  //           3
  //         ],
  //         "r": 20,
  //         "c": 20,
  //         "location": [
  //           0,
  //           0
  //         ]
  //       }
  //     ],
  //     "offset": [
  //       0,
  //       0
  //     ],
  //     "dimension": [
  //       60,
  //       60
  //     ]
  //   }
  // ]
  if (json.contains("storage_zones")) {
    if (json["storage_zones"].is_array()) {
      if (!json["storage_zones"].empty()) {
        for (const auto& zone : json["storage_zones"]) {
          if (zone.contains("slms") && zone["slms"].is_array()) {
            for (const auto& slmSpec : zone["slms"]) {
              storageZones.emplace_back(std::make_unique<SLM>(slmSpec));
            }
          } else {
            throw std::invalid_argument(
                "Storage zone configuration must contain an array of slms");
          }
        }
      } else {
        throw std::invalid_argument(
            "Storage zone configuration must contain at least one zone in "
            "architecture spec");
      }
    } else {
      throw std::invalid_argument(
          "Storage zone configuration must be an array in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "Storage zone configuration is missed in architecture spec");
  }
  // check if the entanglement zones exist and are valid, otherwise throw an
  // error
  // JSON Example:
  // "entanglement_zones": [
  //   {
  //     "zone_id": 0,
  //     "slms": [
  //       {
  //         "id": 1,
  //         "site_separation": [
  //           12,
  //           10
  //         ],
  //         "r": 4,
  //         "c": 4,
  //         "location": [
  //           5,
  //           70
  //         ]
  //       },
  //       {
  //         "id": 2,
  //         "site_separation": [
  //           12,
  //           10
  //         ],
  //         "r": 4,
  //         "c": 4,
  //         "location": [
  //           7,
  //           70
  //         ]
  //       }
  //     ],
  //     "offset": [
  //       5,
  //       70
  //      ],
  //     "dimension": [
  //       50,
  //       40
  //     ]
  //   }
  // ]
  if (json.contains("entanglement_zones")) {
    if (json["entanglement_zones"].is_array()) {
      if (!json["entanglement_zones"].empty()) {
        for (const auto& zone : json["entanglement_zones"]) {
          if (zone.contains("slms") && zone["slms"].is_array()) {
            if (zone["slms"].size() != 2) {
              throw std::invalid_argument("entanglement zone must contain two "
                                          "slms in architecture spec");
            }
            auto& slmPair = *entanglementZones.emplace_back(
                std::make_unique<std::array<SLM, 2>>(std::array<SLM, 2>{
                    SLM(zone["slms"].front()), SLM(zone["slms"].back())}));
            slmPair.front().entanglementId_ = zone["zone_id"];
            slmPair.front().entanglementZone_ = &slmPair;
            slmPair.back().entanglementId_ = zone["zone_id"];
            slmPair.back().entanglementZone_ = &slmPair;
          } else {
            throw std::invalid_argument("Entanglement zone configuration must "
                                        "contain an array of slms");
          }
        }
      } else {
        throw std::invalid_argument("Entanglement zone configuration must "
                                    "contain at least one zone in architecture "
                                    "spec");
      }
    } else {
      throw std::invalid_argument(
          "Entanglement zone configuration must be an array in architecture "
          "spec");
    }
  } else {
    throw std::invalid_argument(
        "Entanglement zone configuration is missed in architecture spec");
  }
  // check if the AODs exist and are valid, otherwise throw an error
  // JSON Example:
  // "aods": [
  //   {
  //     "id": 0,
  //     "site_separation": 2,
  //     "r": 20,
  //     "c": 20
  //   }
  // ]
  if (json.contains("aods")) {
    if (json["aods"].is_array()) {
      for (const auto& aodSpec : json["aods"]) {
        aods.emplace_back(std::make_unique<AOD>(aodSpec));
      }
    } else {
      throw std::invalid_argument(
          "AOD configuration must be an array in architecture spec");
    }
  } else {
    throw std::invalid_argument("AOD is missed in architecture spec");
  }
  // preprocess the created architecture, i.e., calculate the nearest sites for
  // entanglement and storage zones
  preprocessing();
}

auto Architecture::exportNAVizMachine() const -> std::string {
  std::stringstream ss;
  ss << "name: \"" << name << "\"\n";
  ss << "movement {\n    max_speed: 30\n}\n";
  ss << "time {\n    load: 1\n    store: 1\n    ry: 1\n    rz: 1\n    cz: 1\n  "
        "  unit: \"us\"\n}\n";
  ss << "distance {\n    interaction: 10\n    unit: \"um\"\n}\n";
  for (size_t i = 0; i < entanglementZones.size(); ++i) {
    const auto minX = rydbergRangeMinX[i];
    const auto minY = rydbergRangeMinY[i];
    const auto maxX = rydbergRangeMaxX[i];
    const auto maxY = rydbergRangeMaxY[i];
    ss << "zone zone_cz" << *entanglementZones[i]->front().entanglementId_
       << " {\n    from: (" << minX << ", " << minY << ")\n    to: (" << maxX
       << ", " << maxY << ")\n}\n";
  }
  // Generate traps for storage slms
  for (const auto& slm : storageZones) {
    for (std::size_t row = 0; row < slm->nRows; ++row) {
      for (std::size_t col = 0; col < slm->nCols; ++col) {
        const auto& [x, y] = exactSlmLocation(*slm, row, col);
        ss << "trap trap" << slm->id << "_" << row << "_" << col << " {\n";
        ss << "    position: (" << x << ", " << y << ")\n";
        ss << "}\n";
      }
    }
  }
  // do the same for entanglement slms
  for (const auto& zone : entanglementZones) {
    for (const auto& slm : *zone) {
      for (std::size_t row = 0; row < slm.nRows; ++row) {
        for (std::size_t col = 0; col < slm.nCols; ++col) {
          const auto& [x, y] = exactSlmLocation(slm, row, col);
          ss << "trap trap" << slm.id << "_" << row << "_" << col << " {\n";
          ss << "    position: (" << x << ", " << y << ")\n";
          ss << "}\n";
        }
      }
    }
  }
  return ss.str();
}

auto Architecture::isValidSlmPosition(const SLM& slm, const std::size_t r,
                                      const std::size_t c) const -> bool {
  return r < slm.nRows && c < slm.nCols;
}

auto Architecture::exactSlmLocation(const SLM& slm, const std::size_t r,
                                    const std::size_t c) const
    -> std::pair<std::size_t, std::size_t> {
  assert(isValidSlmPosition(slm, r, c));
  return {(slm.siteSeparation.first * c) + slm.location.first,
          (slm.siteSeparation.second * r) + slm.location.second};
}

auto Architecture::findNearestStorageSLM(const size_t x, const size_t y) const
    -> const SLM& {
  double minimumDistance = std::numeric_limits<double>::max();
  const SLM* nearestStorageSLM = nullptr;
  for (const auto& storageSLM : storageZones) {
    std::size_t minimalXDistance = 0;
    if (x < storageSLM->location.first) {
      minimalXDistance = storageSLM->location.first - x;
    } else if (const auto maxX =
                   storageSLM->location.first +
                   ((storageSLM->nCols - 1) * storageSLM->siteSeparation.first);
               x > maxX) {
      minimalXDistance = x - maxX;
    }
    std::size_t minimalYDistance = 0;
    if (y < storageSLM->location.second) {
      minimalYDistance = storageSLM->location.second - y;
    } else if (const auto maxY = storageSLM->location.second +
                                 ((storageSLM->nRows - 1) *
                                  storageSLM->siteSeparation.second);
               y > maxY) {
      minimalYDistance = y - maxY;
    }
    const auto minimalDistance = std::sqrt(std::pow(minimalXDistance, 2) +
                                           std::pow(minimalYDistance, 2));
    if (minimalDistance < minimumDistance) {
      minimumDistance = minimalDistance;
      nearestStorageSLM = storageSLM.get();
    }
  }
  assert(nearestStorageSLM != nullptr);
  return *nearestStorageSLM;
}

auto Architecture::findNearestEntanglementSLM(const size_t x, const size_t y,
                                              const size_t otherX,
                                              const size_t otherY) const
    -> const SLM& {
  // In the loop, we will calculate a lower bound of the distance
  // between the entanglement site and a storage SLM. Any site in
  // the storage SLM will have at least this distance to the
  // entanglement site. This distance will be the variable @c
  // minimalDistance. Among all storage SLMs, we will find the one
  // that has the minimum distance to the entanglement site, the @c
  // minimumDistance.
  double minimumDistance = std::numeric_limits<double>::max();
  const SLM* nearestEntanglementSLM = nullptr;
  for (const auto& entangleSLMs : entanglementZones) {
    for (const auto& entangleSLM : *entangleSLMs) {
      std::size_t minimalXDistance = (x > otherX ? x - otherX : otherX - x);
      if (x < entangleSLM.location.first &&
          otherX < entangleSLM.location.first) {
        minimalXDistance +=
            2 * (entangleSLM.location.first - std::max(x, otherX));
      } else if (const auto maxX = entangleSLM.location.first +
                                   ((entangleSLM.nCols - 1) *
                                    entangleSLM.siteSeparation.first);
                 x > maxX && otherX > maxX) {
        minimalXDistance += 2 * (std::min(x, otherX) - maxX);
      }
      std::size_t minimalYDistance = (y > otherY ? y - otherY : otherY - y);
      if (y < entangleSLM.location.second) {
        minimalYDistance +=
            2 * (entangleSLM.location.second - std::max(y, otherY));
      } else if (const auto maxY = entangleSLM.location.second +
                                   ((entangleSLM.nRows - 1) *
                                    entangleSLM.siteSeparation.second);
                 y > maxY) {
        minimalYDistance += 2 * (std::min(y, otherY) - maxY);
      }
      const auto minimalDistance =
          std::hypot(minimalXDistance, minimalYDistance);
      if (minimalDistance < minimumDistance) {
        minimumDistance = minimalDistance;
        nearestEntanglementSLM = &entangleSLM;
      }
    }
  }
  assert(nearestEntanglementSLM != nullptr);
  return *nearestEntanglementSLM;
}

auto Architecture::preprocessing() -> void {
  //===--------------------------------------------------------------------===//
  // calculate the nearest storage site for each entanglement site
  //===--------------------------------------------------------------------===//
  entanglementToNearestStorageSite.clear();
  for (const auto& slms : entanglementZones) {
    for (const auto& slm : *slms) {
      entanglementToNearestStorageSite.emplace(
          slm,
          std::vector<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                             std::size_t, std::size_t>>>{});
      entanglementToNearestStorageSite.at(slm).reserve(slm.nRows);
      for (std::size_t row = 0; row < slm.nRows; ++row) {
        entanglementToNearestStorageSite.at(slm).emplace_back();
        entanglementToNearestStorageSite.at(slm).back().reserve(slm.nCols);
        for (std::size_t col = 0; col < slm.nCols; ++col) {
          const auto& [x, y] = exactSlmLocation(slm, row, col);
          //===------------------------------------------------------------===//
          // In the first step, find the nearest storage SLM (not the specific
          // site in the storage SLM yet)
          //===------------------------------------------------------------===//
          const auto& nearestStorageSLM = findNearestStorageSLM(x, y);
          //===------------------------------------------------------------===//
          // In the second step, find the specific site in the determined
          // storage SLM
          //===------------------------------------------------------------===//
          std::size_t storageCol = 0;
          if (const auto maxX = nearestStorageSLM.location.first +
                                ((nearestStorageSLM.nCols - 1) *
                                 nearestStorageSLM.siteSeparation.first);
              x > maxX) {
            storageCol = nearestStorageSLM.nCols - 1;
          } else if (x >= nearestStorageSLM.location.first) {
            storageCol = (x - nearestStorageSLM.location.first +
                          (nearestStorageSLM.siteSeparation.first / 2)) /
                         nearestStorageSLM.siteSeparation.first;
          }
          std::size_t storageRow = 0;
          if (const auto maxY = nearestStorageSLM.location.second +
                                ((nearestStorageSLM.nRows - 1) *
                                 nearestStorageSLM.siteSeparation.second);
              y > maxY) {
            storageRow = nearestStorageSLM.nRows - 1;
          } else if (y >= nearestStorageSLM.location.second) {
            storageRow = (y - nearestStorageSLM.location.second +
                          (nearestStorageSLM.siteSeparation.second / 2)) /
                         nearestStorageSLM.siteSeparation.second;
          }
          entanglementToNearestStorageSite.at(slm).back().emplace_back(
              nearestStorageSLM, storageRow, storageCol);
        }
      }
    }
  }
  //===--------------------------------------------------------------------===//
  // calculate the nearest entanglement site for each storage site
  //===--------------------------------------------------------------------===//
  storageToNearestEntanglementSite.clear();
  for (const auto& slm : storageZones) {
    storageToNearestEntanglementSite.emplace(
        *slm,
        std::vector<std::vector<std::unordered_map<
            std::reference_wrapper<const SLM>,
            std::vector<std::vector<
                std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>,
            std::hash<SLM>, std::equal_to<SLM>>>>{});
    storageToNearestEntanglementSite.at(*slm).reserve(slm->nRows);
    for (std::size_t row = 0; row < slm->nRows; ++row) {
      storageToNearestEntanglementSite.at(*slm).emplace_back();
      storageToNearestEntanglementSite.at(*slm).back().reserve(slm->nCols);
      for (std::size_t col = 0; col < slm->nCols; ++col) {
        const auto& [x, y] = exactSlmLocation(*slm, row, col);
        auto& storageToNearestEntanglementSiteForThisSite =
            storageToNearestEntanglementSite.at(*slm).back().emplace_back();
        for (const auto& otherSlm : storageZones) {
          if (otherSlm < slm) {
            continue;
          }
          storageToNearestEntanglementSiteForThisSite.emplace(
              *otherSlm,
              std::vector<
                  std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                         std::size_t, std::size_t>>>{});
          storageToNearestEntanglementSiteForThisSite.at(*otherSlm).reserve(
              slm == otherSlm ? otherSlm->nRows - row : otherSlm->nRows);
          for (std::size_t otherRow = (slm == otherSlm ? row : 0);
               otherRow < otherSlm->nRows; ++otherRow) {
            storageToNearestEntanglementSiteForThisSite[*otherSlm]
                .emplace_back();
            storageToNearestEntanglementSiteForThisSite[*otherSlm]
                .back()
                .reserve(slm == otherSlm && row == otherRow
                             ? otherSlm->nCols - col
                             : otherSlm->nCols);
            for (std::size_t otherCol =
                     (slm == otherSlm && row == otherRow ? col : 0);
                 otherCol < otherSlm->nCols; ++otherCol) {
              const auto& [otherX, otherY] =
                  exactSlmLocation(*otherSlm, otherRow, otherCol);
              //===------------------------------------------------------------===//
              // In the first step, find the nearest storage SLM (not the
              // specific site in the storage SLM yet)
              //===------------------------------------------------------------===//
              const auto& nearestEntanglementSLM =
                  findNearestEntanglementSLM(x, y, otherX, otherY);
              //===------------------------------------------------------------===//
              // In the second step, find the specific site in the determined
              // storage SLM
              //===------------------------------------------------------------===//
              std::size_t entangleCol = 0;
              if (const auto maxX =
                      nearestEntanglementSLM.location.first +
                      ((nearestEntanglementSLM.nCols - 1) *
                       nearestEntanglementSLM.siteSeparation.first);
                  x > maxX) {
                entangleCol = nearestEntanglementSLM.nCols - 1;
              } else if (x >= nearestEntanglementSLM.location.first) {
                entangleCol =
                    (x - nearestEntanglementSLM.location.first +
                     (nearestEntanglementSLM.siteSeparation.first / 2)) /
                    nearestEntanglementSLM.siteSeparation.first;
              }
              std::size_t entangleRow = 0;
              if (const auto maxY =
                      nearestEntanglementSLM.location.second +
                      ((nearestEntanglementSLM.nRows - 1) *
                       nearestEntanglementSLM.siteSeparation.second);
                  y > maxY) {
                entangleRow = nearestEntanglementSLM.nRows - 1;
              } else if (y >= nearestEntanglementSLM.location.second) {
                entangleRow =
                    (y - nearestEntanglementSLM.location.second +
                     (nearestEntanglementSLM.siteSeparation.second / 2)) /
                    nearestEntanglementSLM.siteSeparation.second;
              }
              storageToNearestEntanglementSiteForThisSite.at(*otherSlm)
                  .back()
                  .emplace_back(nearestEntanglementSLM, entangleRow,
                                entangleCol);
            }
          }
        }
      }
    }
  }
}

auto Architecture::distance(const SLM& idx1, const std::size_t r1,
                            const std::size_t c1, const SLM& idx2,
                            const std::size_t r2, const std::size_t c2) const
    -> double {
  const auto& [x1, y1] = exactSlmLocation(idx1, r1, c1);
  const auto& [x2, y2] = exactSlmLocation(idx2, r2, c2);
  return std::hypot(static_cast<double>(x1) - static_cast<double>(x2),
                    static_cast<double>(y1) - static_cast<double>(y2));
}

auto Architecture::nearestStorageSite(const SLM& slm, const std::size_t r,
                                      const std::size_t c) const -> const
    std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t>& {
  return entanglementToNearestStorageSite.at(std::cref(slm))[r][c];
}

auto Architecture::nearestEntanglementSite(
    const SLM& idx1, const std::size_t r1, const std::size_t c1,
    const SLM& idx2, const std::size_t r2, const std::size_t c2) const -> const
    std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t>& {
  assert(idx1.isStorage() && idx2.isStorage());
  if (&idx1 > &idx2 || (&idx1 == &idx2 && r1 > r2) ||
      (&idx1 == &idx2 && r1 == r2 && c1 > c2)) {
    return nearestEntanglementSite(idx2, r2, c2, idx1, r1, c1);
  }
  assert(&idx1 < &idx2 || (&idx1 == &idx2 && r1 < r2) ||
         (&idx1 == &idx2 && r1 == r2 && c1 <= c2));
  return storageToNearestEntanglementSite.at(std::cref(idx1))[r1][c1].at(
      std::cref(idx2))[&idx1 == &idx2 ? r2 - r1 : r2]
                      [&idx1 == &idx2 && r1 == r2 ? c2 - c1 : c2];
}

auto Architecture::nearestEntanglementSiteDistance(
    const SLM& slm1, const std::size_t r1, const std::size_t c1,
    const SLM& slm2, const std::size_t r2, const std::size_t c2) const
    -> double {
  const auto& [x1, y1] = exactSlmLocation(slm1, r1, c1);
  const auto& [x2, y2] = exactSlmLocation(slm2, r2, c2);
  const auto& [entangleSlm, entangleRow, entangleCol] =
      nearestEntanglementSite(slm1, r1, c1, slm2, r2, c2);
  const auto& [entangleX, entangleY] =
      exactSlmLocation(entangleSlm, entangleRow, entangleCol);
  auto dis = std::numeric_limits<double>::max();
  if (r1 == r2 and &slm1 == &slm2) {
    dis = std::min(
        std::max(std::hypot(
                     static_cast<double>(x1) - static_cast<double>(entangleX),
                     static_cast<double>(y1) - static_cast<double>(entangleY)),
                 std::hypot(
                     static_cast<double>(x2) - static_cast<double>(entangleX),
                     static_cast<double>(y2) - static_cast<double>(entangleY))),
        dis);
  } else {
    dis = std::min(
        std::hypot(static_cast<double>(x1) - static_cast<double>(entangleX),
                   static_cast<double>(y1) - static_cast<double>(entangleY)) +
            std::hypot(static_cast<double>(x2) - static_cast<double>(entangleX),
                       static_cast<double>(y2) -
                           static_cast<double>(entangleY)),
        dis);
  }
  return dis;
}

auto Architecture::otherEntanglementSite(const SLM& slm, std::size_t r,
                                         std::size_t c) const
    -> std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t> {
  assert(slm.entanglementZone_);
  const auto& otherSlm = &slm.entanglementZone_->front() == &slm
                             ? slm.entanglementZone_->back()
                             : slm.entanglementZone_->front();
  assert(slm.nCols == otherSlm.nCols);
  assert(slm.nRows == otherSlm.nRows);
  return {otherSlm, r, c};
}

} // namespace na::zoned
