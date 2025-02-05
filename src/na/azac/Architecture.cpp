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
  if (aodSpec.contains("site_seperation")) {
    if (aodSpec["site_seperation"].is_number()) {
      site_separation = aodSpec["site_seperation"];
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
      n_r = aodSpec["r"];
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
      n_c = aodSpec["c"];
    } else {
      throw std::invalid_argument(
          "AOD column number must be a number in architecture spec");
    }
  } else {
    throw std::invalid_argument(
        "AOD column number is missed in architecture spec");
  }
}

SLM::SLM(nlohmann::json slmSpec, const decltype(entanglement_id) entanglementId)
    : entanglement_id(entanglementId) {
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
  if (slmSpec.contains("site_seperation")) {
    if (slmSpec["site_seperation"].is_array() &&
        slmSpec["site_seperation"].size() == 2 &&
        slmSpec["site_seperation"][0].is_number() &&
        slmSpec["site_seperation"][1].is_number()) {
      site_separation = slmSpec["site_seperation"];
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
      n_r = slmSpec["r"];
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
      n_c = slmSpec["c"];
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
auto Architecture::load(const nlohmann::json&& architectureSpec) -> void {;
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
          time_rydberg = architectureSpec["operation_duration"]["rydberg"];
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
          time_atom_transfer =
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
        time_1qGate = architectureSpec["operation_duration"]["1qGate"];
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
      arch_range_min_x = architectureSpec["arch_range"][0][0];
      arch_range_max_x = architectureSpec["arch_range"][0][1];
      arch_range_min_y = architectureSpec["arch_range"][1][0];
      arch_range_max_y = architectureSpec["arch_range"][1][1];
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
      for (const auto& rydberg_range : architectureSpec["rydberg_range"]) {
        if (rydberg_range.is_array() && rydberg_range.size() == 2 &&
            rydberg_range[0].is_array() && rydberg_range[0].size() == 2 &&
            rydberg_range[1].is_array() && rydberg_range[1].size() == 2 &&
            rydberg_range[0][0].is_number() &&
            rydberg_range[0][1].is_number() &&
            rydberg_range[1][0].is_number() &&
            rydberg_range[1][1].is_number()) {
          rydberg_range_min_x.emplace_back(rydberg_range[0][0]);
          rydberg_range_max_x.emplace_back(rydberg_range[0][1]);
          rydberg_range_min_y.emplace_back(rydberg_range[1][0]);
          rydberg_range_max_y.emplace_back(rydberg_range[1][1]);
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
          for (const auto& slm_spec : zone["slms"]) {
            storage_zone.emplace_back(std::make_unique<SLM>(slm_spec));
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
              ySlm.emplace(y, &entanglement_zone.emplace_back());
            }
            ySlm[y]->emplace_back(std::make_unique<SLM>(slmSpec, ySlm[y]));
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
        dict_AOD.emplace_back(std::make_unique<AOD>(aodSpec));
      }
    } else {
      throw std::invalid_argument(
          "AOD configuration must be an array in architecture spec");
    }
  } else {
    throw std::invalid_argument("AOD is missed in architecture spec");
  }
}
auto Architecture::is_valid_SLM_position(const SLM& slm, const std::size_t r,
                                         const std::size_t c) const -> bool {
  return r < slm.n_r && c < slm.n_c;
}
auto Architecture::exact_SLM_location(const SLM& slm, const std::size_t r,
                                      const std::size_t c) const
    -> std::pair<std::size_t, std::size_t> {
  assert(is_valid_SLM_position(slm, r, c));
  return {(slm.site_separation.first * c) + slm.location.first,
          (slm.site_separation.second * r) + slm.location.second};
}
auto Architecture::preprocessing() -> void {
  //===--------------------------------------------------------------------===//
  // calculate the nearest storage site for each entanglement site
  //===--------------------------------------------------------------------===//
  entanglementToNearestStorageSite.clear();
  for (const auto& slms : entanglement_zone) {
    for (const auto& slm : slms) {
      entanglementToNearestStorageSite.emplace(
          slm.get(),
          std::vector<std::vector<
              std::tuple<const SLM* const, std::size_t, std::size_t>>>{});
      entanglementToNearestStorageSite[slm.get()].reserve(slm->n_r);
      for (std::size_t row = 0; row < slm->n_r; ++row) {
        entanglementToNearestStorageSite[slm.get()].emplace_back();
        entanglementToNearestStorageSite[slm.get()].back().reserve(slm->n_c);
        for (std::size_t col = 0; col < slm->n_c; ++col) {
          const auto& [x, y] = exact_SLM_location(*slm, row, col);
          //===------------------------------------------------------------===//
          // In the first step, find the nearest storage SLM (not the specific
          // site in the storage SLM yet)
          //===------------------------------------------------------------===//
          // In the loop, we will calculate a lower bound of the distance
          // between the entanglement site and a storage SLM. Any site in the
          // storage SLM will have at least this distance to the entanglement
          // site. This distance will be the variable @c minimalDistance.
          // Among all storage SLMs, we will find the one that has the minimum
          // distance to the entanglement site, the @c minimumDistance.
          double minimumDistance = std::numeric_limits<double>::max();
          const SLM* nearestStorageSLM = nullptr;
          for (const auto& storageSLM : storage_zone) {
            std::size_t minimalXDistance = 0;
            if (x < storageSLM->location.first) {
              minimalXDistance = storageSLM->location.first - x;
            } else if (const auto maxX = storageSLM->location.first +
                                         ((storageSLM->n_c - 1) *
                                          storageSLM->site_separation.first);
                       x > maxX) {
              minimalXDistance = x - maxX;
            }
            std::size_t minimalYDistance = 0;
            if (y < storageSLM->location.second) {
              minimalYDistance = storageSLM->location.second - y;
            } else if (const auto maxY = storageSLM->location.second +
                                         ((storageSLM->n_r - 1) *
                                          storageSLM->site_separation.second);
                       y > maxY) {
              minimalYDistance = y - maxY;
            }
            const auto minimalDistance = std::sqrt(
                std::pow(minimalXDistance, 2) + std::pow(minimalYDistance, 2));
            if (minimalDistance < minimumDistance) {
              minimumDistance = minimalDistance;
              nearestStorageSLM = storageSLM.get();
            }
          }
          //===------------------------------------------------------------===//
          // In the second step, find the specific site in the determined
          // storage SLM
          //===------------------------------------------------------------===//
          std::size_t storageCol = 0;
          if (const auto maxX = nearestStorageSLM->location.first +
                                ((nearestStorageSLM->n_c - 1) *
                                 nearestStorageSLM->site_separation.first);
              x > maxX) {
            storageCol = nearestStorageSLM->n_c - 1;
          } else if (x >= nearestStorageSLM->location.first) {
            storageCol = (x - nearestStorageSLM->location.first +
                          (nearestStorageSLM->site_separation.first / 2)) /
                         nearestStorageSLM->site_separation.first;
          }
          std::size_t storageRow = 0;
          if (const auto maxY = nearestStorageSLM->location.second +
                                ((nearestStorageSLM->n_r - 1) *
                                 nearestStorageSLM->site_separation.second);
              y > maxY) {
            storageRow = nearestStorageSLM->n_r - 1;
          } else if (y >= nearestStorageSLM->location.second) {
            storageRow = (y - nearestStorageSLM->location.second +
                          (nearestStorageSLM->site_separation.second / 2)) /
                         nearestStorageSLM->site_separation.second;
          }
          entanglementToNearestStorageSite[slm.get()].back().emplace_back(
              nearestStorageSLM, storageRow, storageCol);
        }
      }
    }
  }
  //===--------------------------------------------------------------------===//
  // calculate the nearest entanglement site for each storage site
  //===--------------------------------------------------------------------===//
  storageToNearestEntanglementSite.clear();
  storageToNearestEntanglementSiteDistance.clear();
  for (const auto& slm : storage_zone) {
    storageToNearestEntanglementSite.emplace(
        slm.get(),
        std::vector<std::vector<
            std::tuple<const SLM* const, std::size_t, std::size_t>>>{});
    storageToNearestEntanglementSite[slm.get()].reserve(slm->n_r);
    storageToNearestEntanglementSiteDistance.emplace(
        slm.get(), std::vector<std::vector<double>>{});
    storageToNearestEntanglementSiteDistance[slm.get()].reserve(slm->n_r);
    for (std::size_t row = 0; row < slm->n_r; ++row) {
      storageToNearestEntanglementSite[slm.get()].emplace_back();
      storageToNearestEntanglementSite[slm.get()].back().reserve(slm->n_c);
      storageToNearestEntanglementSiteDistance[slm.get()].emplace_back();
      storageToNearestEntanglementSiteDistance[slm.get()].back().reserve(
          slm->n_c);
      for (std::size_t col = 0; col < slm->n_c; ++col) {
        const auto& [x, y] = exact_SLM_location(*slm, row, col);
        //===--------------------------------------------------------------===//
        // in the first step, find the nearest storage SLM (not the specific
        // site in the storage SLM yet)
        //===--------------------------------------------------------------===//
        // in the loop, we will calculate a lower bound of the distance
        // between the storage site and an entanglement SLM. Any site in the
        // entanglement SLM will have at least this distance to the storage
        // site. This distance will be the variable @c minimalDistance.
        // Among all entanglement SLMs, we will find the one that has the
        // minimum distance to the storage site, the @c minimumDistance.
        double minimumDistance = std::numeric_limits<double>::max();
        const SLM* nearestEntanglementSLM = nullptr;
        for (const auto& entanglementSLMs : entanglement_zone) {
          for (const auto& entanglementSLM : entanglementSLMs) {
            std::size_t minimalXDistance = 0;
            if (x < entanglementSLM->location.first) {
              minimalXDistance = entanglementSLM->location.first - x;
            } else if (const auto maxX =
                           entanglementSLM->location.first +
                           ((entanglementSLM->n_c - 1) *
                            entanglementSLM->site_separation.first);
                       x > maxX) {
              minimalXDistance = x - maxX;
            }
            std::size_t minimalYDistance = 0;
            if (y < entanglementSLM->location.second) {
              minimalYDistance = entanglementSLM->location.second - y;
            } else if (const auto maxY =
                           entanglementSLM->location.second +
                           ((entanglementSLM->n_r - 1) *
                            entanglementSLM->site_separation.second);
                       y > maxY) {
              minimalYDistance = y - maxY;
            }
            const auto minimalDistance = std::sqrt(
                std::pow(minimalXDistance, 2) + std::pow(minimalYDistance, 2));
            if (minimalDistance < minimumDistance) {
              minimumDistance = minimalDistance;
              nearestEntanglementSLM = entanglementSLM.get();
            }
          }
          //===------------------------------------------------------------===//
          // In the second step, find the specific site in the determined
          // storage SLM
          //===------------------------------------------------------------===//
          std::size_t storageCol = 0;
          if (const auto maxX = nearestEntanglementSLM->location.first +
                                ((nearestEntanglementSLM->n_c - 1) *
                                 nearestEntanglementSLM->site_separation.first);
              x > maxX) {
            storageCol = nearestEntanglementSLM->n_c - 1;
          } else if (x >= nearestEntanglementSLM->location.first) {
            storageCol = (x - nearestEntanglementSLM->location.first +
                          (nearestEntanglementSLM->site_separation.first / 2)) /
                         nearestEntanglementSLM->site_separation.first;
          }
          std::size_t storageRow = 0;
          if (const auto maxY =
                  nearestEntanglementSLM->location.second +
                  ((nearestEntanglementSLM->n_r - 1) *
                   nearestEntanglementSLM->site_separation.second);
              y > maxY) {
            storageRow = nearestEntanglementSLM->n_r - 1;
          } else if (y >= nearestEntanglementSLM->location.second) {
            storageRow =
                (y - nearestEntanglementSLM->location.second +
                 (nearestEntanglementSLM->site_separation.second / 2)) /
                nearestEntanglementSLM->site_separation.second;
          }
          storageToNearestEntanglementSite[slm.get()].back().emplace_back(
              nearestEntanglementSLM, storageRow, storageCol);
          storageToNearestEntanglementSiteDistance[slm.get()]
              .back()
              .emplace_back(distance(slm.get(), row, col,
                                     nearestEntanglementSLM, storageRow,
                                     storageCol));
        }
      }
    }
  }
}
auto Architecture::distance(const SLM& idx1, const std::size_t r1,
                            const std::size_t c1, const SLM& idx2,
                            const std::size_t r2, const std::size_t c2) const
    -> double {
  return na::distance(exact_SLM_location(idx1, r1, c1),
                      exact_SLM_location(idx2, r2, c2));
}
auto Architecture::nearest_storage_site(const SLM& slm, const std::size_t r,
                                        const std::size_t c) const
    -> std::tuple<const SLM* const, std::size_t, std::size_t> {
  return entanglementToNearestStorageSite.at(&slm)[r][c];
}
auto Architecture::nearest_entanglement_site(const SLM* const idx,
                                             const std::size_t r,
                                             const std::size_t c) const
    -> std::tuple<const SLM* const, std::size_t, std::size_t> {
  return storageToNearestEntanglementSite.at(idx)[r][c];
}
auto Architecture::nearest_entanglement_site_distance(const SLM* const idx,
                                                      const std::size_t r,
                                                      const std::size_t c) const
    -> double {
  return storageToNearestEntanglementSiteDistance.at(idx)[r][c];
}
auto Architecture::nearest_entanglement_site(
    const SLM* const idx1, const std::size_t r1, const std::size_t c1,
    const SLM* const idx2, const std::size_t r2, const std::size_t c2) const
    -> std::tuple<const SLM* const, std::size_t, std::size_t> {
  const auto& site1 = storageToNearestEntanglementSite.at(idx1)[r1][c1];
  const auto& site2 = storageToNearestEntanglementSite.at(idx2)[r2][c2];
  // the nearest entanglement zone for both qubits is the same
  if (site1 == site2) {
    return site1;
  }
  if (std::get<0>(site1) == std::get<0>(site2)) {
    const auto middle_site_c = (std::get<2>(site1) + std::get<2>(site2)) / 2;
    return {std::get<0>(site1), std::get<1>(site1), middle_site_c};
  }
  throw std::invalid_argument(
      "[ERROR] The nearest entanglement site is not in the same entanglement "
      "zone. This feature is not supported yet.");
}
auto Architecture::nearest_entanglement_site_distance(
    const SLM* const slm1, const std::size_t r1, const std::size_t c1,
    const SLM* const slm2, const std::size_t r2, const std::size_t c2) const
    -> double {
  const auto& storageSite1 = exact_SLM_location(slm1, r1, c1);
  const auto& storageSite2 = exact_SLM_location(slm2, r2, c2);
  const auto& entanglementSite =
      exact_SLM_location(nearest_entanglement_site(slm1, r1, c1, slm2, r2, c2));
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
auto Architecture::movement_duration(const std::size_t x1, const std::size_t y1,
                                     const std::size_t x2, const std::size_t y2)
    -> double {
  // todo: add reference for constant
  constexpr double a = 0.00275;
  const auto d = na::distance(std::pair{x1, y1}, std::pair{x2, y2});
  const auto t = std::sqrt(d / a);
  return t;
}

} // namespace na