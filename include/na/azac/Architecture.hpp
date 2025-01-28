#pragma once
#include <filesystem>
#include <fstream>
#include <list>
#include <memory>
#include <nlohmann/json.hpp>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {
/// An 2D-Array of AOD traps
struct AOD {
  std::size_t id = 0;
  std::size_t site_seperation_x = 0;
  std::size_t site_seperation_y = 0;
  std::size_t n_r = 0;
  std::size_t n_c = 0;

  explicit AOD(nlohmann::json aod_spec) {
    if (aod_spec.contains("id")) {
      id = aod_spec["id"];
    } else {
      throw std::invalid_argument("AOD id is missed in architecture spec");
    }
    if (aod_spec.contains("site_seperation")) {
      site_seperation_x = aod_spec["site_seperation"][0];
      site_seperation_y = aod_spec["site_seperation"][1];
    } else {
      throw std::invalid_argument(
          "AOD site seperation is missed in architecture spec");
    }
    if (aod_spec.contains("r")) {
      n_r = aod_spec["r"];
    } else {
      throw std::invalid_argument(
          "AOD row number is missed in architecture spec");
    }
    if (aod_spec.contains("c")) {
      n_c = aod_spec["c"];
    } else {
      throw std::invalid_argument(
          "AOD column number is missed in architecture spec");
    }
  }
};

/// An 2D-array of SLM traps
struct SLM {
  std::size_t id = 0;
  std::size_t site_seperation_x = 0;
  std::size_t site_seperation_y = 0;
  std::size_t n_r = 0;
  std::size_t n_c = 0;
  std::size_t location_x = 0; ///< x-coordinate of the left uppermost SLM
  std::size_t location_y = 0; ///< y-coordinate of the left uppermost SLM
  std::optional<std::vector<std::unique_ptr<SLM>>&> entanglement_id = std::nullopt;

  explicit SLM(nlohmann::json slm_spec, decltype(entanglement_id) entanglement_id = std::nullopt)
      : entanglement_id(entanglement_id) {
    if (slm_spec.contains("id")) {
      id = slm_spec["id"];
    } else {
      throw std::invalid_argument("SLM id is missed in architecture spec");
    }
    if (slm_spec.contains("site_seperation")) {
      site_seperation_x = slm_spec["site_seperation"][0];
      site_seperation_y = slm_spec["site_seperation"][1];
    } else {
      throw std::invalid_argument(
          "SLM site seperation is missed in architecture spec");
    }
    if (slm_spec.contains("r")) {
      n_r = slm_spec["r"];
    } else {
      throw std::invalid_argument(
          "SLM row number is missed in architecture spec");
    }
    if (slm_spec.contains("c")) {
      n_c = slm_spec["c"];
    } else {
      throw std::invalid_argument(
          "SLM column number is missed in architecture spec");
    }
    if (slm_spec.contains("location")) {
      location_x = slm_spec["location"][0];
      location_y = slm_spec["location"][1];
    } else {
      throw std::invalid_argument(
          "SLM location is missed in architecture spec");
    }
  }
};

/// Class to define zone architecture
struct Architecture {
  std::string name;
  std::unordered_map<std::string, double> operation_duration;
  std::vector<std::unique_ptr<SLM>> storage_zone;
  std::vector<std::vector<std::unique_ptr<SLM>>> entanglement_zone;
  std::vector<std::unique_ptr<AOD>> dict_AOD;
  double time_atom_transfer = 15; ///< µs
  double time_rydberg = 0.36;     ///< µs
  double time_1qGate = 0.625;     ///< µs
  std::size_t arch_range_min_x = 0;
  std::size_t arch_range_max_x = 0;
  std::size_t arch_range_min_y = 0;
  std::size_t arch_range_max_y = 0;
  std::size_t rydberg_range_min_x = 0;
  std::size_t rydberg_range_max_x = 0;
  std::size_t rydberg_range_min_y = 0;
  std::size_t rydberg_range_max_y = 0;
  std::unordered_map<const SLM*, std::vector<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>> Rydberg_site_nearest_storage_site{};
  std::unordered_map<const SLM*, std::vector<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>> storage_site_nearest_Rydberg_site{};
  std::unordered_map<const SLM*, std::vector<std::vector<double>>> storage_site_nearest_Rydberg_site_dis{};

  Architecture() = default;
  explicit Architecture(const std::string& filename)
      : Architecture(std::filesystem::path(filename)) {}
  explicit Architecture(const std::filesystem::path& filepath)
      : Architecture(std::ifstream(filepath)) {}
  explicit Architecture(std::ifstream& ifs) : Architecture(std::move(ifs)) {}
  explicit Architecture(std::ifstream&& ifs) { load(std::move(ifs)); }
  auto load(const std::string& filename) -> void {
    load(std::filesystem::path(filename));
  }
  auto load(const std::filesystem::path& filepath) -> void {
    load(std::ifstream(filepath));
  }
  auto load(std::ifstream& ifs) -> void { load(std::move(ifs)); }
  auto load(std::ifstream&& ifs) -> void {
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
      std::unordered_map<std::size_t, std::vector<std::unique_ptr<SLM>>&>
          y_slm{};
      for (const auto& zone : architecture_spec["entanglement_zones"]) {
        for (const auto& slm_spec : zone["slms"]) {
          const std::size_t y = slm_spec["location"][1];
          auto it = y_slm.find(y);
          if (it == y_slm.end()) {
            y_slm[y] = entanglement_zone.emplace_back();
          }
          y_slm[y].emplace_back(std::make_unique<SLM>(slm_spec, y_slm[y]));
        }
      }
    } else {
      throw std::invalid_argument(
          "entanglement zone configuration is missed in architecture spec");
    }
    if (architecture_spec.contains("aods")) {
      for (const auto& aod_spec : architecture_spec["aods"]) {
        dict_AOD.emplace_back(std::make_unique<AOD>(aod_spec));
      }
    } else {
      throw std::invalid_argument("AOD is missed in architecture spec");
    }
  }

  auto is_valid_SLM_position(const SLM* slm, const std::size_t r,
                             const std::size_t c) -> bool {
    return is_valid_SLM_position(*slm, r, c);
  }

  auto is_valid_SLM_position(const SLM& slm, const std::size_t r,
                             const std::size_t c) -> bool {
    return r < slm.n_r && c < slm.n_c;
  }

  auto exact_SLM_location(const SLM* slm, const std::size_t r,
                          const std::size_t c)
      -> std::pair<std::size_t, std::size_t> {
    return exact_SLM_location(*slm, r, c);
  }

  auto exact_SLM_location(const SLM& slm, const std::size_t r,
                          const std::size_t c)
      -> std::pair<std::size_t, std::size_t> {
    assert(is_valid_SLM_position(slm, r, c));
    return {(slm.site_seperation_x * c) + slm.location_x,
            (slm.site_seperation_y * r) + slm.location_y};
  }

  /// Compute the site region for entanglement zone and the nearest Rydberg site
  /// for each storage site. We assume, we only have one storage zone or one
  /// entanglement zone per row
  auto preprocessing() -> void {
    std::vector<std::pair<std::size_t, const SLM*>> entanglement_site_row_space; // 2d array. [[y, idx]] => if y' < y, the nearest row to y is the row in zone idx
    std::unordered_map<const SLM*, std::vector<std::size_t>> entanglement_site_col_space;
    std::vector<std::pair<std::size_t, std::size_t>> y_site;
    for (std::size_t i = 0; i < entanglement_zone.size(); ++i) {
      const auto& slm = *entanglement_zone[i].front();
      y_site.emplace_back(slm.location_y, i);
    }
    std::sort(y_site.begin(), y_site.end(),
              [](const std::pair<std::size_t, std::size_t>& a,
                 const std::pair<std::size_t, std::size_t>& b) {
                return a.first < b.first;
              });

    for (std::size_t i = 0; i < y_site.size() - 1; ++i) {
      const auto* slm = entanglement_zone[y_site[i].second].front().get();
      const auto low_y = y_site[i].first + slm->site_seperation_y * (slm->n_r - 1);
      const auto high_y = y_site[i + 1].first;
      entanglement_site_row_space.emplace_back((high_y + low_y) / 2, slm);
    }
    entanglement_site_row_space.emplace_back(std::numeric_limits<std::size_t>::max(), entanglement_zone[y_site[-1].second].front().get());
    // split the column area for SLM sites
    for (const auto& list_idx: entanglement_zone) {
      const auto* idx = list_idx[0].get();
      entanglement_site_col_space.emplace(idx, std::vector<std::size_t>{});
      const auto x = idx->location_x + idx->site_seperation_x / 2;
      for (std::size_t c = 0; c < idx->n_c - 1; ++c) {
        entanglement_site_col_space[idx].emplace_back(x + c * idx->site_seperation_x);
      }
      entanglement_site_col_space[idx].emplace_back(std::numeric_limits<std::size_t>::max());
    }
    // compute the nearest Rydberg site for each storage site
    storage_site_nearest_Rydberg_site.clear();
    storage_site_nearest_Rydberg_site_dis.clear();

    Rydberg_site_nearest_storage_site.clear();
    for (const auto& list_idx : entanglement_zone) {
      const auto* idx = list_idx.front().get();
      Rydberg_site_nearest_storage_site.emplace(idx, std::vector(2, std::vector(idx->n_c, std::tuple<const SLM*, std::size_t, std::size_t>{nullptr, 0, 0})));
    }


  for (const auto& idx : storage_zone) {
    storage_site_nearest_Rydberg_site.emplace(idx.get(), std::vector(idx->n_r, std::vector(idx->n_c, std::tuple<const SLM*, std::size_t, std::size_t>{nullptr, 0, 0})));
    storage_site_nearest_Rydberg_site_dis.emplace(idx.get(), std::vector(idx->n_r, std::vector(idx->n_c, 0.0)));
    const auto x = idx->location_x;
    auto y = idx->location_y;
    const auto* nearest_slm = entanglement_site_row_space.back().second;
    const auto nearest_slm_half_r = nearest_slm->n_r / 2;
    const SLM* next_nearest_slm = nullptr;
    auto y_lim = entanglement_site_row_space.back().first;
    const auto row_y_l = nearest_slm->location_y;
    const auto row_y = row_y_l + (nearest_slm->n_r - 1) * nearest_slm->site_seperation_y;
    std::size_t row = (y > row_y_l ? y - row_y_l : row_y_l - y) < (y > row_y ? y - row_y : row_y - y) ? 0 : nearest_slm->n_r - 1;
    bool has_increase_y = false;
    // find the entanglement slm for the row
    for (std::size_t i = 0; entanglement_site_row_space.size() - 1; ++i) {
      if (y < entanglement_site_row_space[i].first) {
        nearest_slm = entanglement_site_row_space[i].second;
        next_nearest_slm = entanglement_site_row_space[i+1].second;
        y_lim = entanglement_site_row_space[i].first;
        row = nearest_slm->n_r - 1;
        has_increase_y = true;
        break;
      }
    }
    const auto init_x = x;
    auto init_x_lim = entanglement_site_col_space[nearest_slm].back();
    auto init_col = nearest_slm->n_c - 1;
    for (std::size_t i = 0; i < entanglement_site_col_space[nearest_slm].size(); ++i) {
      if (x < entanglement_site_col_space[nearest_slm][i]) {
        init_x_lim = entanglement_site_col_space[nearest_slm][i];
        init_col = i;
        break;
      }
    }
    for (std::size_t r = 0; r < idx->n_r; ++r) {
      auto x_lim = init_x_lim;
      auto col = init_col;
      auto x = init_x;
      for (std::size_t c = 0; c < idx->n_c; ++c) {
        storage_site_nearest_Rydberg_site[idx.get()][r][c] = {nearest_slm, row, col};
        storage_site_nearest_Rydberg_site_dis[idx.get()][r][c] = distance(idx.get(), r, c, nearest_slm, row, col);
        const auto r_idx = r < nearest_slm_half_r ? 0 : 1;
        if (std::get<0>(Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col]) == nullptr) {
          Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col] = {idx.get(), r, c};
        } else {
          const auto& [prev_idx, prev_r, prev_c] = Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col];
          const auto prev_dis = storage_site_nearest_Rydberg_site_dis[prev_idx][prev_r][prev_c];
          if (prev_dis > storage_site_nearest_Rydberg_site_dis[idx.get()][r][c]) {
            Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col] = {idx.get(), r, c};
          }
        }

        x += idx->site_seperation_x;
        if (x > x_lim && col + 1 < nearest_slm->n_c) {
          col += 1;
          x_lim = entanglement_site_col_space[nearest_slm][col];
        }
      }
      y += idx->site_seperation_y;
      if (has_increase_y && y > y_lim && next_nearest_slm != nullptr) {
        has_increase_y = false;
        nearest_slm = next_nearest_slm;
        row = 0;
      }
    }
  }

  for (const auto& [key, _] : Rydberg_site_nearest_storage_site) {
      // check if row_0 is all -1:
      std::vector<std::optional<std::size_t>> idx_first_none_empty{std::nullopt, std::nullopt};
      std::vector<std::optional<std::size_t>> idx_last_none_empty{std::nullopt, std::nullopt};
      for (std::size_t i = 0; i < Rydberg_site_nearest_storage_site[key].size(); ++i) {
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
      for (std::size_t i = 0; i< Rydberg_site_nearest_storage_site[key].size(); ++i) {
        if (idx_first_none_empty[i] && idx_first_none_empty[i] > 0) {
          for (std::size_t j = *idx_first_none_empty[i] - 1; j >= 0; --j) {
            Rydberg_site_nearest_storage_site[key][i][j] = Rydberg_site_nearest_storage_site[key][i][*idx_first_none_empty[i]];
          }
        }
        if (idx_last_none_empty[i]) {
          for (std::size_t j = *idx_last_none_empty[i] + 1; j < Rydberg_site_nearest_storage_site[key][i].size(); ++j) {
            Rydberg_site_nearest_storage_site[key][i][j] = Rydberg_site_nearest_storage_site[key][i][*idx_last_none_empty[i]];
          }
        }
      }
      if (!idx_first_none_empty[0] && !idx_last_none_empty[0]) {
        Rydberg_site_nearest_storage_site[key][0] = Rydberg_site_nearest_storage_site[key][1];
      } else if (!idx_first_none_empty[1] && !idx_last_none_empty[1]) {
        Rydberg_site_nearest_storage_site[key][1] = Rydberg_site_nearest_storage_site[key][0];
      }
    }
  }

  auto distance(const SLM* idx1, std::size_t r1, std::size_t c1,
                const SLM* idx2, std::size_t r2, std::size_t c2) -> double {
    const auto& p1 = exact_SLM_location(idx1, r1, c1);
    const auto& p2 = exact_SLM_location(idx2, r2, c2);
    const auto dx =
        static_cast<double>(p1.first) - static_cast<double>(p2.first);
    const auto dy =
        static_cast<double>(p1.second) - static_cast<double>(p2.second);
    return std::sqrt((dx * dx) + (dy * dy));
  }

  auto nearest_storage_site(const SLM* idx, std::size_t r, std::size_t c)
      -> std::tuple<const SLM*, std::size_t, std::size_t> {
    const auto* slm = idx->entanglement_id->front().get();
    const auto nearest_slm_half_r = slm->n_r / 2;
    if (r < nearest_slm_half_r) {
      return Rydberg_site_nearest_storage_site[slm][0][c];
    }
    return Rydberg_site_nearest_storage_site[slm][1][c];
  }

  /// return the nearest Rydberg site for a qubit in the storage zone
  auto nearest_entanglement_site(const SLM* idx, std::size_t r, std::size_t c)
      -> std::tuple<const SLM*, std::size_t, std::size_t> {
    return storage_site_nearest_Rydberg_site[idx][r][c];
  }

  /// return the distance nearest Rydberg site for a qubit in the storage zone
  auto nearest_entanglement_site_distance(const SLM* idx, std::size_t r,
                                          std::size_t c) -> double {
    return storage_site_nearest_Rydberg_site_dis[idx][r][c];
  }

  /// return the nearest Rydberg site for two qubit in the storage zone
  /// based on the position of two qubits
  auto nearest_entanglement_site(const SLM* idx1, std::size_t r1,
                                 std::size_t c1, const SLM* idx2,
                                 std::size_t r2, std::size_t c2)
      -> std::tuple<const SLM*, std::size_t, std::size_t> {
    const auto& storage_site1 = exact_SLM_location(idx1, r1, c1);
    const auto& storage_site2 = exact_SLM_location(idx2, r2, c2);
    const auto& site1 = storage_site_nearest_Rydberg_site[idx1][r1][c1];
    const auto& site2 = storage_site_nearest_Rydberg_site[idx2][r2][c2];
    // the nearest zone for both qubits are in the same entanglement zone
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

  /// return the sum of the distance to move two qubits to one rydberg site
  auto nearest_entanglement_site_dis(const SLM* idx1, std::size_t r1,
                                     std::size_t c1, const SLM* idx2,
                                     std::size_t r2, std::size_t c2) -> double {
    const auto& storage_site1 = exact_SLM_location(idx1, r1, c1);
    const auto& storage_site2 = exact_SLM_location(idx2, r2, c2);
    const auto& list_site =
        nearest_entanglement_site(idx1, r1, c1, idx2, r2, c2);
    const auto& dis = std::numeric_limits<double>::max();
    for (const auto site : list_site) {
      const std::pair<size_t, size_t>& exact_site =
          std::apply(exact_SLM_location, site);
      if (r1 == r2 and idx1 == idx2) {
        const auto dx1 = static_cast<double>(storage_site1.first) -
                         static_cast<double>(exact_site.first);
        const auto dy1 = static_cast<double>(storage_site1.second) -
                         static_cast<double>(exact_site.second);
        const auto dx2 = static_cast<double>(storage_site2.first) -
                         static_cast<double>(exact_site.first);
        const auto dy2 = static_cast<double>(storage_site2.second) -
                         static_cast<double>(exact_site.second);
        dis = std::min(std::max(std::sqrt((dx1 * dx1) + (dy1 * dy1)),
                                std::sqrt((dx2 * dx2) + (dy2 * dy2))),
                       dis);
      } else {
        const auto dx1 = static_cast<double>(storage_site1.first) -
                         static_cast<double>(exact_site.first);
        const auto dy1 = static_cast<double>(storage_site1.second) -
                         static_cast<double>(exact_site.second);
        const auto dx2 = static_cast<double>(storage_site2.first) -
                         static_cast<double>(exact_site.first);
        const auto dy2 = static_cast<double>(storage_site2.second) -
                         static_cast<double>(exact_site.second);
        dis = std::min(std::sqrt((dx1 * dx1) + (dy1 * dy1)) +
                           std::sqrt((dx2 * dx2) + (dy2 * dy2)),
                       dis);
      }
    }
    return dis;
  }

  auto movement_duration(std::size_t x1, std::size_t y1, std::size_t x2,
                         std::size_t y2) -> double {
    // todo: add reference for constant
    constexpr double a = 0.00275;
    const auto dx = static_cast<double>(x1) - static_cast<double>(x2);
    const auto dy = static_cast<double>(y1) - static_cast<double>(y2);
    const auto d = std::sqrt((dx * dx) + (dy * dy));
    const auto t = std::sqrt(d / a);
    return t;
  }
};
} // namespace na