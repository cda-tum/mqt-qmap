#include "na/azac/Architecture.hpp"

#include "na/azac/Utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <istream>
#include <limits>
#include <memory>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace na {
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

SLM::SLM(nlohmann::json slmSpec,
         const decltype(entanglementZone) entanglementZone,
         const std::size_t entanglementId)
    : SLM(slmSpec) {
  if (entanglementZone == nullptr) {
    throw std::invalid_argument("If set, entanglementZone must not be nullptr");
  }
  this->entanglementZone = entanglementZone;
  this->entanglementId = entanglementId;
}
auto Architecture::load(const nlohmann::json&& architectureSpec) -> void {
  ;
  if (architectureSpec.contains("name")) {
    if (architectureSpec["name"].is_string()) {
      name = architectureSpec["name"];
    } else {
      throw std::invalid_argument("architecture name must be a string");
    }
  } else {
    throw std::invalid_argument("architecture name is missed in architecture "
                                "spec");
  }
  if (architectureSpec.contains("operation_duration")) {
    if (architectureSpec["operation_duration"].is_object()) {
      if (architectureSpec["operation_duration"].contains("rydberg")) {
        if (architectureSpec["operation_duration"]["rydberg"].is_number()) {
          timeRydberg = architectureSpec["operation_duration"]["rydberg"];
        } else {
          throw std::invalid_argument(
              "rydberg duration must be a number in architecture spec");
        }
      } else {
        throw std::invalid_argument(
            "operation duration must contain rydberg duration");
      }
      if (architectureSpec["operation_duration"].contains("atom_transfer")) {
        if (architectureSpec["operation_duration"]["atom_transfer"]
                .is_number()) {
          timeAtomTransfer =
              architectureSpec["operation_duration"]["atom_transfer"];
        } else {
          throw std::invalid_argument(
              "atom transfer duration must be a number in architecture spec");
        }
      } else {
        throw std::invalid_argument(
            "operation duration must contain atom transfer duration");
      }
      if (architectureSpec["operation_duration"].contains("1qGate")) {
        time1QGate = architectureSpec["operation_duration"]["1qGate"];
      } else {
        throw std::invalid_argument(
            "operation duration must contain 1qGate duration");
      }
    } else {
      throw std::invalid_argument(
          "operation duration must be an dict in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "operation duration is missed in architecture spec");
  }
  if (architectureSpec.contains("arch_range")) {
    if (architectureSpec["arch_range"].is_array() &&
        architectureSpec["arch_range"].size() == 2 &&
        architectureSpec["arch_range"][0].is_array() &&
        architectureSpec["arch_range"][0].size() == 2 &&
        architectureSpec["arch_range"][1].is_array() &&
        architectureSpec["arch_range"][1].size() == 2 &&
        architectureSpec["arch_range"][0][0].is_number() &&
        architectureSpec["arch_range"][0][1].is_number() &&
        architectureSpec["arch_range"][1][0].is_number() &&
        architectureSpec["arch_range"][1][1].is_number()) {
      archRangeMinX = architectureSpec["arch_range"][0][0];
      archRangeMaxX = architectureSpec["arch_range"][0][1];
      archRangeMinY = architectureSpec["arch_range"][1][0];
      archRangeMaxY = architectureSpec["arch_range"][1][1];
    } else {
      throw std::invalid_argument(
          "architecture range must be a 2x2 number array in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "architecture range is missed in architecture spec");
  }
  if (architectureSpec.contains("rydberg_range")) {
    if (architectureSpec["rydberg_range"].is_array() &&
        !architectureSpec["rydberg_range"].empty()) {
      for (const auto& rydbergRange : architectureSpec["rydberg_range"]) {
        if (rydbergRange.is_array() && rydbergRange.size() == 2 &&
            rydbergRange[0].is_array() && rydbergRange[0].size() == 2 &&
            rydbergRange[1].is_array() && rydbergRange[1].size() == 2 &&
            rydbergRange[0][0].is_number() && rydbergRange[0][1].is_number() &&
            rydbergRange[1][0].is_number() && rydbergRange[1][1].is_number()) {
          rydbergRangeMinX.emplace_back(rydbergRange[0][0]);
          rydbergRangeMaxX.emplace_back(rydbergRange[0][1]);
          rydbergRangeMinY.emplace_back(rydbergRange[1][0]);
          rydbergRangeMaxY.emplace_back(rydbergRange[1][1]);
        } else {
          throw std::invalid_argument("rydberg range must be a Nx2x2 number "
                                      "array in architecture spec, N > 1");
        }
      }
    } else {
      throw std::invalid_argument("rydberg range must be a Nx2x2 number array "
                                  "in architecture spec, N > 1");
    }
  } else {
    throw std::invalid_argument(
        "architecture range is missed in architecture spec");
  }
  if (architectureSpec.contains("storage_zones")) {
    if (architectureSpec["storage_zones"].is_array()) {
      for (const auto& zone : architectureSpec["storage_zones"]) {
        if (zone.contains("slms") && zone["slms"].is_array()) {
          for (const auto& slmSpec : zone["slms"]) {
            storageZones.emplace_back(std::make_unique<SLM>(slmSpec));
          }
        } else {
          throw std::invalid_argument(
              "storage zone configuration must contain an array of slms");
        }
      }
    } else {
      throw std::invalid_argument(
          "storage zone configuration must be an array in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "storage zone configuration is missed in architecture spec");
  }
  if (architectureSpec.contains("entanglement_zones")) {
    if (architectureSpec["entanglement_zones"].is_array()) {
      // corresponding slms that form an entanglement group/zone are matched by
      // their y-coordinate, i.e., slms with the same y-coordinate are in the
      // same group
      std::unordered_map<std::size_t, std::vector<std::unique_ptr<SLM>>*>
          ySlm{};
      for (const auto& zone : architectureSpec["entanglement_zones"]) {
        if (zone.contains("slms") && zone["slms"].is_array()) {
          for (const auto& slmSpec : zone["slms"]) {
            if (!slmSpec.contains("location") &&
                !slmSpec["location"].is_array() &&
                slmSpec["location"].size() != 2) {
              throw std::invalid_argument("location is missed in slm spec or "
                                          "it is not a 2-element array");
            }
            const std::size_t y = slmSpec["location"][1];
            auto it = ySlm.find(y);
            if (it == ySlm.end()) {
              ySlm.emplace(y, &entanglementZones.emplace_back());
            }
            ySlm[y]->emplace_back(
                std::make_unique<SLM>(slmSpec, ySlm[y], zone["zone_id"]));
          }
        } else {
          throw std::invalid_argument(
              "entanglement zone configuration must contain an array of slms");
        }
      }
    } else {
      throw std::invalid_argument(
          "entanglement zone configuration must be an array in architecture "
          "spec");
    }
  } else {
    throw std::invalid_argument(
        "entanglement zone configuration is missed in architecture spec");
  }
  if (architectureSpec.contains("aods")) {
    if (architectureSpec["aods"].is_array()) {
      for (const auto& aodSpec : architectureSpec["aods"]) {
        aods.emplace_back(std::make_unique<AOD>(aodSpec));
      }
    } else {
      throw std::invalid_argument(
          "AOD configuration must be an array in architecture spec");
    }
  } else {
    throw std::invalid_argument("AOD is missed in architecture spec");
  }
}
auto Architecture::exportNAVizMachine() const -> std::string {
  std::stringstream ss;
  ss << "name: \"" << name << "\"\n";
  ss << "movement {\n    max_speed: 30\n}\n";
  ss << "time {\n    load: 1\n    store: 1\n    ry: 1\n    rz: 1\n    cz: 1\n  "
        "  unit: \"us\"\n}\n";
  ss << "distance {\n    interaction: 10\n    unit: \"um\"\n}\n";
  if (entanglementZones.size() != 1) {
    throw std::runtime_error(
        "Right now, only one entanglement zone is supported");
  }
  const auto& slm1 = entanglementZones.front().front();
  const auto& slm2 = entanglementZones.front().back();
  const auto minX = std::min(slm1->location.first, slm2->location.first);
  const auto minY = std::min(slm1->location.second, slm2->location.second);
  const auto maxX = std::max(
      slm1->location.first + (slm1->nCols - 1) * slm1->siteSeparation.first,
      slm2->location.first + (slm2->nCols - 1) * slm2->siteSeparation.first);
  const auto maxY = std::max(
      slm1->location.second + (slm1->nRows - 1) * slm1->siteSeparation.second,
      slm2->location.second + (slm2->nRows - 1) * slm2->siteSeparation.second);
  ss << "zone zone_cz0 {\n    from: (" << minX - 10 << ", " << minY - 10
     << ")\n    to: (" << maxX + 10 << ", " << maxY + 10 << ")\n}\n";
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
    for (const auto& slm : zone) {
      for (std::size_t row = 0; row < slm->nRows; ++row) {
        for (std::size_t col = 0; col < slm->nCols; ++col) {
          const auto& [x, y] = exactSlmLocation(*slm, row, col);
          ss << "trap trap" << slm->id << "_" << row << "_" << col << " {\n";
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
auto Architecture::preprocessing() -> void {
  //===--------------------------------------------------------------------===//
  // calculate the nearest storage site for each entanglement site
  //===--------------------------------------------------------------------===//
  entanglementToNearestStorageSite.clear();
  for (const auto& slms : entanglementZones) {
    for (const auto& slm : slms) {
      entanglementToNearestStorageSite.emplace(
          slm.get(), std::vector<std::vector<std::vector<
                         std::tuple<const SLM*, std::size_t, std::size_t>>>>{});
      entanglementToNearestStorageSite[slm.get()].reserve(slm->nRows);
      for (std::size_t row = 0; row < slm->nRows; ++row) {
        entanglementToNearestStorageSite[slm.get()].emplace_back();
        entanglementToNearestStorageSite[slm.get()].back().reserve(slm->nCols);
        for (std::size_t col = 0; col < slm->nCols; ++col) {
          const auto& entanglementSiteLocation =
              exactSlmLocation(*slm, row, col);
          //===------------------------------------------------------------===//
          // get all storage sites and put them in a vector with their
          // distance to the entanglement site
          //===------------------------------------------------------------===//
          std::vector<std::pair<
              std::tuple<const SLM*, std::size_t, std::size_t>, double>>
              allStorageSitesWithTheirDistance;
          for (const auto& storageSlm : storageZones) {
            for (std::size_t storageRow = 0; storageRow < storageSlm->nRows;
                 ++storageRow) {
              for (std::size_t storageCol = 0; storageCol < storageSlm->nCols;
                   ++storageCol) {
                const auto distance = na::distance(
                    entanglementSiteLocation,
                    exactSlmLocation(*storageSlm, storageRow, storageCol));
                allStorageSitesWithTheirDistance.emplace_back(
                    std::make_tuple(storageSlm.get(), storageRow, storageCol),
                    distance);
              }
            }
          }
          //===------------------------------------------------------------===//
          // sort the storage sites by their distance to the entanglement site
          //===------------------------------------------------------------===//
          std::sort(
              allStorageSitesWithTheirDistance.begin(),
              allStorageSitesWithTheirDistance.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
          //===------------------------------------------------------------===//
          // put the sorted storage sites in the map
          //===------------------------------------------------------------===//
          entanglementToNearestStorageSite[slm.get()].back().emplace_back();
          entanglementToNearestStorageSite[slm.get()].back().back().reserve(
              allStorageSitesWithTheirDistance.size());
          std::transform(
              allStorageSitesWithTheirDistance.begin(),
              allStorageSitesWithTheirDistance.end(),
              std::back_inserter(
                  entanglementToNearestStorageSite[slm.get()].back().back()),
              [](const auto& t) { return t.first; });
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
        slm.get(),
        std::vector<std::vector<
            std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>>{});
    storageToNearestEntanglementSite[slm.get()].reserve(slm->nRows);
    for (std::size_t row = 0; row < slm->nRows; ++row) {
      storageToNearestEntanglementSite[slm.get()].emplace_back();
      storageToNearestEntanglementSite[slm.get()].back().reserve(slm->nCols);
      for (std::size_t col = 0; col < slm->nCols; ++col) {
        const auto& storageSiteLocation = exactSlmLocation(*slm, row, col);
        //===--------------------------------------------------------------===//
        // get all entanglement sites and put them in a vector with their
        // distance to the storage site
        //===--------------------------------------------------------------===//
        std::vector<
            std::pair<std::tuple<const SLM*, std::size_t, std::size_t>, double>>
            allEntanglementSitesWithTheirDistance;
        for (const auto& slms : entanglementZones) {
          for (const auto& entanglementSlm : slms) {
            for (std::size_t entanglementRow = 0;
                 entanglementRow < entanglementSlm->nRows; ++entanglementRow) {
              for (std::size_t entanglementCol = 0;
                   entanglementCol < entanglementSlm->nCols;
                   ++entanglementCol) {
                const auto distance = na::distance(
                    storageSiteLocation,
                    exactSlmLocation(*entanglementSlm, entanglementRow,
                                     entanglementCol));
                allEntanglementSitesWithTheirDistance.emplace_back(
                    std::make_tuple(entanglementSlm.get(), entanglementRow,
                                    entanglementCol),
                    distance);
              }
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
  return na::distance(exactSlmLocation(idx1, r1, c1),
                      exactSlmLocation(idx2, r2, c2));
}
auto Architecture::nearestStorageSite(const SLM& slm, const std::size_t r,
                                      const std::size_t c) const
    -> const std::tuple<const SLM*, std::size_t, std::size_t>& {
  return entanglementToNearestStorageSite.at(&slm)[r][c].front();
}
auto Architecture::nearestStorageSitesAsc(const SLM& slm, std::size_t r,
                                          std::size_t c) const
    -> const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>& {
  return entanglementToNearestStorageSite.at(&slm)[r][c];
}
auto Architecture::nearestEntanglementSite(const SLM* const idx,
                                           const std::size_t r,
                                           const std::size_t c) const
    -> const std::tuple<const SLM*, std::size_t, std::size_t>& {
  return storageToNearestEntanglementSite.at(idx)[r][c].front();
}
auto Architecture::nearestEntanglementSiteDistance(const SLM* const idx,
                                                   const std::size_t r,
                                                   const std::size_t c) const
    -> double {
  return distance({idx, r, c},
                  storageToNearestEntanglementSite.at(idx)[r][c].front());
}
auto Architecture::nearestEntanglementSite(
    const SLM* const idx1, const std::size_t r1, const std::size_t c1,
    const SLM* const idx2, const std::size_t r2, const std::size_t c2) const
    -> std::tuple<const SLM*, std::size_t, std::size_t> {
  const auto& site1 = storageToNearestEntanglementSite.at(idx1)[r1][c1].front();
  const auto& site2 = storageToNearestEntanglementSite.at(idx2)[r2][c2].front();
  // the nearest entanglement zone for both qubits is the same
  if (site1 == site2) {
    return site1;
  }
  if (std::get<0>(site1) == std::get<0>(site2)) {
    const auto middleSite_c = (std::get<2>(site1) + std::get<2>(site2)) / 2;
    return {std::get<0>(site1), std::get<1>(site1), middleSite_c};
  }
  throw std::invalid_argument(
      "[ERROR] The nearest entanglement site is not in the same entanglement "
      "zone. This feature is not supported yet.");
}
auto Architecture::nearestEntanglementSiteDistance(
    const SLM* const slm1, const std::size_t r1, const std::size_t c1,
    const SLM* const slm2, const std::size_t r2, const std::size_t c2) const
    -> double {
  const auto& storageSite1 = exactSlmLocation(slm1, r1, c1);
  const auto& storageSite2 = exactSlmLocation(slm2, r2, c2);
  const auto& entanglementSite =
      exactSlmLocation(nearestEntanglementSite(slm1, r1, c1, slm2, r2, c2));
  auto dis = std::numeric_limits<double>::max();
  if (r1 == r2 and slm1 == slm2) {
    dis = std::min(std::max(na::distance(storageSite1, entanglementSite),
                            na::distance(storageSite2, entanglementSite)),
                   dis);
  } else {
    dis = std::min(na::distance(storageSite1, entanglementSite) +
                       na::distance(storageSite2, entanglementSite),
                   dis);
  }
  return dis;
}
auto Architecture::movementDuration(const std::size_t x1, const std::size_t y1,
                                    const std::size_t x2, const std::size_t y2)
    -> double {
  // todo: add reference for constant
  constexpr double a = 0.00275;
  const auto d = na::distance(std::pair{x1, y1}, std::pair{x2, y2});
  const auto t = std::sqrt(d / a);
  return t;
}

} // namespace na
