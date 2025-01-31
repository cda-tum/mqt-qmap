#include "na/azac/Architecture.hpp"

#include "na/azac/Util.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
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
    id = aodSpec["id"];
  } else {
    throw std::invalid_argument("AOD id is missed in architecture spec");
  }
  if (aodSpec.contains("site_seperation")) {
    site_separation = aodSpec["site_seperation"];
  } else {
    throw std::invalid_argument(
        "AOD site separation is missed in architecture spec");
  }
  if (aodSpec.contains("r")) {
    n_r = aodSpec["r"];
  } else {
    throw std::invalid_argument(
        "AOD row number is missed in architecture spec");
  }
  if (aodSpec.contains("c")) {
    n_c = aodSpec["c"];
  } else {
    throw std::invalid_argument(
        "AOD column number is missed in architecture spec");
  }
}

SLM::SLM(nlohmann::json slmSpec, const decltype(entanglement_id) entanglementId)
    : entanglement_id(entanglementId) {
  if (slmSpec.contains("id")) {
    id = slmSpec["id"];
  } else {
    throw std::invalid_argument("SLM id is missed in architecture spec");
  }
  if (slmSpec.contains("site_seperation")) {
    site_separation = slmSpec["site_seperation"];
  } else {
    throw std::invalid_argument(
        "SLM site separation is missed in architecture spec");
  }
  if (slmSpec.contains("r")) {
    n_r = slmSpec["r"];
  } else {
    throw std::invalid_argument(
        "SLM row number is missed in architecture spec");
  }
  if (slmSpec.contains("c")) {
    n_c = slmSpec["c"];
  } else {
    throw std::invalid_argument(
        "SLM column number is missed in architecture spec");
  }
  if (slmSpec.contains("location")) {
    location = slmSpec["location"];
  } else {
    throw std::invalid_argument("SLM location is missed in architecture spec");
  }
}
auto Architecture::load(std::ifstream&& ifs) -> void {
  nlohmann::json architecture_spec{};
  std::move(ifs) >> architecture_spec;
  if (architecture_spec.contains("name")) {
    name = architecture_spec["name"];
  }
  if (architecture_spec.contains("operation_duration")) {
    if (architecture_spec["operation_duration"].contains("rydberg")) {
      operation_duration["rydberg"] =
          architecture_spec["operation_duration"]["rydberg"];
      time_rydberg = architecture_spec["operation_duration"]["rydberg"];
    }
    if (architecture_spec["operation_duration"].contains("atom_transfer")) {
      operation_duration["atom_transfer"] =
          architecture_spec["operation_duration"]["atom_transfer"];
      time_atom_transfer =
          architecture_spec["operation_duration"]["atom_transfer"];
    }
    if (architecture_spec["operation_duration"].contains("1qGate")) {
      operation_duration["1qGate"] =
          architecture_spec["operation_duration"]["1qGate"];
      time_1qGate = architecture_spec["operation_duration"]["1qGate"];
    }
  }
  if (architecture_spec.contains("arch_range")) {
    arch_range_min_x = architecture_spec["arch_range"][0][0];
    arch_range_max_x = architecture_spec["arch_range"][0][1];
    arch_range_min_y = architecture_spec["arch_range"][1][0];
    arch_range_max_y = architecture_spec["arch_range"][1][1];
  } else {
    throw std::invalid_argument(
        "architecture range is missed in architecture spec");
  }
  if (architecture_spec.contains("rydberg_range")) {
    rydberg_range_min_x = architecture_spec["rydberg_range"][0][0];
    rydberg_range_max_x = architecture_spec["rydberg_range"][0][1];
    rydberg_range_min_y = architecture_spec["rydberg_range"][1][0];
    rydberg_range_max_y = architecture_spec["rydberg_range"][1][1];
  } else {
    throw std::invalid_argument(
        "architecture range is missed in architecture spec");
  }
  if (architecture_spec.contains("storage_zones")) {
    for (const auto& zone : architecture_spec["storage_zones"]) {
      for (const auto& slm_spec : zone["slms"]) {
        storage_zone.emplace_back(std::make_unique<SLM>(slm_spec));
      }
    }
  } else {
    throw std::invalid_argument(
        "storage zone configuration is missed in architecture spec");
  }
  if (architecture_spec.contains("entanglement_zones")) {
    std::unordered_map<std::size_t, std::vector<std::unique_ptr<SLM>>&> ySlm{};
    for (const auto& zone : architecture_spec["entanglement_zones"]) {
      for (const auto& slmSpec : zone["slms"]) {
        const std::size_t y = slmSpec["location"][1];
        auto it = ySlm.find(y);
        if (it == ySlm.end()) {
          ySlm[y] = entanglement_zone.emplace_back();
        }
        ySlm[y].emplace_back(std::make_unique<SLM>(slmSpec, &ySlm[y]));
      }
    }
  } else {
    throw std::invalid_argument(
        "entanglement zone configuration is missed in architecture spec");
  }
  if (architecture_spec.contains("aods")) {
    for (const auto& aodSpec : architecture_spec["aods"]) {
      dict_AOD.emplace_back(std::make_unique<AOD>(aodSpec));
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
  // for every pair (d, slm) records the vertical space below the slm to
  // the next one, where d is half of the total distance
  std::vector<std::pair<std::size_t, const SLM* const>>
      entanglement_site_row_space;
  // same as above for the horizontal space
  std::unordered_map<const SLM* const, std::vector<std::size_t>>
      entanglement_site_col_space;
  //===------------------------------------------------------------------===//
  // split the vertical area for SLM sites
  //===------------------------------------------------------------------===//
  // In the entanglement zone, SLM are offset in the x direction
  // Following we collect the first SLM for each pair of SLM in the
  // entanglement zone
  std::vector<const SLM*> y_site;
  for (const auto& slms : entanglement_zone) {
    y_site.emplace_back(slms.front().get());
  }
  // we sort the slm in y_site by their y location in ascending order
  std::sort(y_site.begin(), y_site.end(),
            [](const SLM* const a, const SLM* const b) {
              return a->location.second < b->location.second;
            });
  for (std::size_t i = 0; i < y_site.size() - 1; ++i) {
    const auto low_y =
        y_site[i]->location.second +
        (y_site[i]->site_separation.second * (y_site[i]->n_r - 1));
    const auto high_y = y_site[i + 1]->location.second;
    entanglement_site_row_space.emplace_back((high_y + low_y) / 2, y_site[i]);
  }
  // for the slm with the highest y location, insert a fake entry
  // of the maximum y value
  entanglement_site_row_space.emplace_back(
      std::numeric_limits<std::size_t>::max(), y_site.back());
  //===------------------------------------------------------------------===//
  // split the horizontal area for SLM sites
  //===------------------------------------------------------------------===//
  for (const auto& slms : entanglement_zone) {
    entanglement_site_col_space.emplace(slms.front().get(),
                                        std::vector<std::size_t>{});
    const auto x = slms.front()->location.first +
                   (slms.front()->site_separation.first / 2);
    for (std::size_t c = 0; c < slms.front()->n_c - 1; ++c) {
      entanglement_site_col_space[slms.front().get()].emplace_back(
          x + (c * slms.front()->site_separation.first));
    }
    entanglement_site_col_space[slms.front().get()].emplace_back(
        std::numeric_limits<std::size_t>::max());
  }
  //===------------------------------------------------------------------===//
  // compute the nearest Rydberg site for each storage site
  //===------------------------------------------------------------------===//
  storage_site_nearest_Rydberg_site.clear();
  storage_site_nearest_Rydberg_site_dis.clear();
  Rydberg_site_nearest_storage_site.clear();
  // initialize the vectors in Rydberg_site_nearest_storage_site such that
  // in consecutive accesses we can assume the existence of the respective
  // vectors and do not insert the element one by one
  for (const auto& slms : entanglement_zone) {
    Rydberg_site_nearest_storage_site.emplace(
        slms.front().get(),
        std::vector(
            2,
            std::vector(slms.front()->n_c,
                        std::tuple<const SLM* const, std::size_t, std::size_t>{
                            nullptr, 0, 0})));
  }
  // for all slm in the storage zones
  for (const auto& slm : storage_zone) {
    // initialize the vectors in storage_site_nearest_Rydberg_site such that
    // in consecutive accesses we can assume the existence of the respective
    // vectors and do not insert the element one by one
    storage_site_nearest_Rydberg_site.emplace(
        slm.get(),
        std::vector(
            slm->n_r,
            std::vector(slm->n_c,
                        std::tuple<const SLM* const, std::size_t, std::size_t>{
                            nullptr, 0, 0})));
    storage_site_nearest_Rydberg_site_dis.emplace(
        slm.get(), std::vector(slm->n_r, std::vector(slm->n_c, 0.0)));
    // in the following we find the nearest entanglement zone
    // for the current storage slm
    auto [loc_x, loc_y] = slm->location;
    const auto* nearest_slm = entanglement_site_row_space.back().second;
    const auto nearest_slm_half_r = nearest_slm->n_r / 2;
    const SLM* next_nearest_slm = nullptr;
    auto y_lim = entanglement_site_row_space.back().first;
    const auto row_y_l = nearest_slm->location.second;
    const auto row_y = row_y_l + ((nearest_slm->n_r - 1) *
                                  nearest_slm->site_separation.second);
    std::size_t row = (loc_y > row_y_l ? loc_y - row_y_l : row_y_l - loc_y) <
                              (loc_y > row_y ? loc_y - row_y : row_y - loc_y)
                          ? 0
                          : nearest_slm->n_r - 1;
    bool has_increase_y = false;
    // find the entanglement slm for the row
    for (std::size_t i = 0; i < entanglement_site_row_space.size() - 1; ++i) {
      if (loc_y < entanglement_site_row_space[i].first) {
        nearest_slm = entanglement_site_row_space[i].second;
        next_nearest_slm = entanglement_site_row_space[i + 1].second;
        y_lim = entanglement_site_row_space[i].first;
        row = nearest_slm->n_r - 1;
        has_increase_y = true;
        break;
      }
    }
    const auto init_x = loc_x;
    auto init_x_lim = entanglement_site_col_space[nearest_slm].back();
    auto init_col = nearest_slm->n_c - 1;
    for (std::size_t i = 0; i < entanglement_site_col_space[nearest_slm].size();
         ++i) {
      if (loc_x < entanglement_site_col_space[nearest_slm][i]) {
        init_x_lim = entanglement_site_col_space[nearest_slm][i];
        init_col = i;
        break;
      }
    }
    for (std::size_t r = 0; r < slm->n_r; ++r) {
      auto x_lim = init_x_lim;
      auto col = init_col;
      auto x = init_x;
      for (std::size_t c = 0; c < slm->n_c; ++c) {
        storage_site_nearest_Rydberg_site[*slm][r][c] = {nearest_slm, row, col};
        storage_site_nearest_Rydberg_site_dis[*slm][r][c] =
            distance(*slm, r, c, *nearest_slm, row, col);
        const std::size_t r_idx = r < nearest_slm_half_r ? 0 : 1;
        if (std::get<0>(
                Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col]) ==
            nullptr) {
          Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col] = {
              slm.get(), r, c};
        } else {
          const auto& [prev_idx, prev_r, prev_c] =
              Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col];
          const auto prev_dis =
              storage_site_nearest_Rydberg_site_dis[prev_idx][prev_r][prev_c];
          if (prev_dis >
              storage_site_nearest_Rydberg_site_dis[slm.get()][r][c]) {
            Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col] = {
                slm.get(), r, c};
          }
        }

        x += slm->site_separation_x;
        if (x > x_lim && col + 1 < nearest_slm->n_c) {
          col += 1;
          x_lim = entanglement_site_col_space[nearest_slm][col];
        }
      }
      loc_y += slm->site_separation_y;
      if (has_increase_y && loc_y > y_lim && next_nearest_slm != nullptr) {
        has_increase_y = false;
        nearest_slm = next_nearest_slm;
        row = 0;
      }
    }
  }

  for (const auto& [key, _] : Rydberg_site_nearest_storage_site) {
    // check if row_0 is all -1:
    std::vector<std::optional<std::size_t>> idx_first_none_empty{std::nullopt,
                                                                 std::nullopt};
    std::vector<std::optional<std::size_t>> idx_last_none_empty{std::nullopt,
                                                                std::nullopt};
    for (std::size_t i = 0; i < Rydberg_site_nearest_storage_site[key].size();
         ++i) {
      const auto& row = Rydberg_site_nearest_storage_site[key][i];
      for (std::size_t j = 0; j < row.size(); ++j) {
        const auto& [site, r, c] = row[j];
        if (site != nullptr) {
          if (!idx_first_none_empty[i]) {
            idx_first_none_empty[i] = j;
          }
          idx_last_none_empty[i] = j;
        }
      }
    }
    // assume at least one row is non empty
    if (idx_first_none_empty[0] == -1 and idx_last_none_empty[0] == -1 and
        idx_first_none_empty[1] == -1 and idx_last_none_empty[1] == -1) {
      throw std::runtime_error("[ERROR] Both rows are empty");
    }
    for (std::size_t i = 0; i < Rydberg_site_nearest_storage_site[key].size();
         ++i) {
      if (idx_first_none_empty[i] && idx_first_none_empty[i] > 0) {
        for (std::size_t j = *idx_first_none_empty[i] - 1; j >= 0; --j) {
          Rydberg_site_nearest_storage_site[key][i][j] =
              Rydberg_site_nearest_storage_site[key][i]
                                               [*idx_first_none_empty[i]];
        }
      }
      if (idx_last_none_empty[i]) {
        for (std::size_t j = *idx_last_none_empty[i] + 1;
             j < Rydberg_site_nearest_storage_site[key][i].size(); ++j) {
          Rydberg_site_nearest_storage_site[key][i][j] =
              Rydberg_site_nearest_storage_site[key][i]
                                               [*idx_last_none_empty[i]];
        }
      }
    }
    if (!idx_first_none_empty[0] && !idx_last_none_empty[0]) {
      Rydberg_site_nearest_storage_site[key][0] =
          Rydberg_site_nearest_storage_site[key][1];
    } else if (!idx_first_none_empty[1] && !idx_last_none_empty[1]) {
      Rydberg_site_nearest_storage_site[key][1] =
          Rydberg_site_nearest_storage_site[key][0];
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
  // get the unique first slm in the group of entanglement slms
  const auto& uniqueSlm = *slm.entanglement_id->front();
  return Rydberg_site_nearest_storage_site.at(
      &uniqueSlm)[r < uniqueSlm.n_r / 2 ? 0 : 1][c];
}
auto Architecture::nearest_entanglement_site(const SLM* const idx,
                                             const std::size_t r,
                                             const std::size_t c) const
    -> std::tuple<const SLM* const, std::size_t, std::size_t> {
  return storage_site_nearest_Rydberg_site.at(idx)[r][c];
}
auto Architecture::nearest_entanglement_site_distance(const SLM* const idx,
                                                      const std::size_t r,
                                                      const std::size_t c) const
    -> double {
  return storage_site_nearest_Rydberg_site_dis.at(idx)[r][c];
}
auto Architecture::nearest_entanglement_site(
    const SLM* const idx1, const std::size_t r1, const std::size_t c1,
    const SLM* const idx2, const std::size_t r2, const std::size_t c2) const
    -> std::tuple<const SLM* const, std::size_t, std::size_t> {
  const auto& site1 = storage_site_nearest_Rydberg_site.at(idx1)[r1][c1];
  const auto& site2 = storage_site_nearest_Rydberg_site.at(idx2)[r2][c2];
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