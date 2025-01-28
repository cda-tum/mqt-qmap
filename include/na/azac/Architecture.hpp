#pragma once
#include <filesystem>
#include <fstream>
#include <list>
#include <nlohmann/json.hpp>
#include <unordered_set>
#include <utility>

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
  std::optional<std::size_t> entanglement_id = std::nullopt;

  explicit SLM(nlohmann::json slm_spec, bool storage = true) {
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
    if (!storage) {
      if (slm_spec.contains("entanglement_id")) {
        entanglement_id = slm_spec["entanglement_id"];
      } else {
        throw std::invalid_argument(
            "SLM entanglement id is missed in architecture spec");
      }
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
          const auto slm = std::make_unique<SLM>(slm_spec, false);
          if (y_slm.find(slm->location_y) == y_slm.end()) {
            y_slm[slm->location_y] = entanglement_zone.emplace_back();
          }
          y_slm[slm->location_y].emplace_back(slm);
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
  auto preprocessing() -> void {}

  /*
  #split the row area for SLM sites
          self.entanglement_site_row_space = [] // 2d array. [[y, idx]] => if y'
  < y, the nearest row to y is the row in zone idx
          self.entanglement_site_col_space = dict()
  #print("self.entanglement_zone")
  #print(self.entanglement_zone)
          y_site = []
          for i, idx in enumerate(self.entanglement_zone):
              slm = self.dict_SLM[idx[0]]
              y_site.append((slm.location[1], i))
          y_site = sorted(y_site, key=lambda site: site[0])

          for i in range(len(y_site) - 1):
              slm = self.dict_SLM[self.entanglement_zone[y_site[i][1]][0]]
              low_y = y_site[i][0] + slm.site_seperation[1] * (slm.n_r - 1)
              high_y = y_site[i+1][0]
              self.entanglement_site_row_space.append(((high_y + low_y) / 2,
  slm.idx)) self.entanglement_site_row_space.append((sys.maxsize,
  self.dict_SLM[self.entanglement_zone[y_site[-1][1]][0]].idx)) #split the
  column area for SLM sites for list_idx in self.entanglement_zone: idx =
  list_idx[0] self.entanglement_site_col_space[idx] = [] slm =
  self.dict_SLM[idx] x = slm.location[0] + slm.site_seperation[0] / 2
              self.entanglement_site_col_space[idx] = [x + c *
  slm.site_seperation[0] for c in range(slm.n_c - 1)]
              self.entanglement_site_col_space[idx].append(sys.maxsize)
  #print("self.entanglement_site_row_space")
  #print(self.entanglement_site_row_space)
  #print("self.entanglement_site_col_space")
  #print(self.entanglement_site_col_space)
  #compute the nearest Rydberg site for each storage site
          self.storage_site_nearest_Rydberg_site = dict()
          self.storage_site_nearest_Rydberg_site_dis = dict()

          self.Rydberg_site_nearest_storage_site = dict()
          for list_idx in self.entanglement_zone:
              idx = list_idx[0]
              slm = self.dict_SLM[idx]
              self.Rydberg_site_nearest_storage_site[idx] = [[-1 for j in
  range(slm.n_c)] for i in range(2)]



          for idx in self.storage_zone:
              slm = self.dict_SLM[idx]
              self.storage_site_nearest_Rydberg_site[idx] = [[0 for j in
  range(slm.n_c)] for i in range(slm.n_r)]
              self.storage_site_nearest_Rydberg_site_dis[idx] = [[0 for j in
  range(slm.n_c)] for i in range(slm.n_r)] x, y = slm.location nearest_slm =
  self.entanglement_site_row_space[-1][1] nearest_slm_half_r =
  self.dict_SLM[nearest_slm].n_r // 2 next_nearest_slm = -1 y_lim =
  self.entanglement_site_row_space[-1][0] row_y_l =
  self.dict_SLM[nearest_slm].location[1] row_y = row_y_l +
  (self.dict_SLM[nearest_slm].n_r - 1) * slm.site_seperation[1] if abs(y -
  row_y_l) < abs(y - row_y): row = 0 else: row = self.dict_SLM[nearest_slm].n_r
  - 1 has_increase_y = False #find the entanglement slm for the row for i in
  range(len(self.entanglement_site_row_space) - 1): if y <
  self.entanglement_site_row_space[i][0]: nearest_slm =
  self.entanglement_site_row_space[i][1] next_nearest_slm =
  self.entanglement_site_row_space[i+1][1] y_lim =
  self.entanglement_site_row_space[i][0] row = self.dict_SLM[nearest_slm].n_r -
  1 has_increase_y = True break init_x = x init_x_lim =
  self.entanglement_site_col_space[nearest_slm][-1] init_col =
  self.dict_SLM[nearest_slm].n_c - 1 for i in
  range(len(self.entanglement_site_col_space[nearest_slm])): if x <
  self.entanglement_site_col_space[nearest_slm][i]: init_x_lim =
  self.entanglement_site_col_space[nearest_slm][i] init_col = i break

              for r in range(slm.n_r):
                  x_lim = init_x_lim
                  col = init_col
                  x = init_x
                  for c in range(slm.n_c):
                      self.storage_site_nearest_Rydberg_site[idx][r][c] =
  (nearest_slm, row, col) self.storage_site_nearest_Rydberg_site_dis[idx][r][c]
  = self.distance(idx, r, c, nearest_slm, row, col)

                      if row < nearest_slm_half_r:
                          r_idx = 0
                      else:
                          r_idx = 1
                      if
  self.Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col] == -1:
                          self.Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col]
  = (idx, r, c) else: (prev_idx, prev_r, prev_c) =
  self.Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col] prev_dis =
  self.storage_site_nearest_Rydberg_site_dis[prev_idx][prev_r][prev_c] if
  prev_dis > self.storage_site_nearest_Rydberg_site_dis[idx][r][c]:
                              self.Rydberg_site_nearest_storage_site[nearest_slm][r_idx][col]
  = (idx, r, c)

                      x += slm.site_seperation[0]
                      if x > x_lim and col + 1 < self.dict_SLM[nearest_slm].n_c:
                          col += 1
                          x_lim =
  self.entanglement_site_col_space[nearest_slm][col] y += slm.site_seperation[1]
                  if has_increase_y and y > y_lim and next_nearest_slm > -1:
                      has_increase_y = False
                      nearest_slm = next_nearest_slm
                      row = 0

          for key in self.Rydberg_site_nearest_storage_site:
  #check if row_0 is all - 1:
              idx_first_none_empty = [-1, -1]
              idx_last_none_empty = [-1, -1]
              for i, row in
  enumerate(self.Rydberg_site_nearest_storage_site[key]): for j, site in
  enumerate(row): if site != -1: if idx_first_none_empty[i] == -1:
                              idx_first_none_empty[i] = j
                          idx_last_none_empty[i] = j
  #assume at least one row is non empty
              if idx_first_none_empty[0] == -1 and idx_last_none_empty[0] == -1\
                  and idx_first_none_empty[1] == -1 and idx_last_none_empty[1]
  == -1: assert(0) for i in
  range(len(self.Rydberg_site_nearest_storage_site[key])): if
  idx_first_none_empty[i] != -1: for j in range(idx_first_none_empty[i] - 1, -1,
  -1): self.Rydberg_site_nearest_storage_site[key][i][j] =
  self.Rydberg_site_nearest_storage_site[key][i][idx_first_none_empty[i]] if
  idx_last_none_empty[i] != -1: for j in range(idx_last_none_empty[i] + 1,
  len(self.Rydberg_site_nearest_storage_site[key][i])):
                          self.Rydberg_site_nearest_storage_site[key][i][j] =
  self.Rydberg_site_nearest_storage_site[key][i][idx_last_none_empty[i]] if
  idx_first_none_empty[0] == -1 and idx_last_none_empty[0] == -1:
                  self.Rydberg_site_nearest_storage_site[key][0] =
  self.Rydberg_site_nearest_storage_site[key][1] elif idx_first_none_empty[1] ==
  -1 and idx_last_none_empty[1] == -1:
                  self.Rydberg_site_nearest_storage_site[key][1] =
  self.Rydberg_site_nearest_storage_site[key][0]

  #print("self.storage_site_nearest_Rydberg_site")
  #for key in self.storage_site_nearest_Rydberg_site:
  #print(key)
  #for site in self.storage_site_nearest_Rydberg_site[key]:
  #print(site)
  #print("self.storage_site_nearest_Rydberg_site_dis")
  #print(self.storage_site_nearest_Rydberg_site_dis)
  #print("self.Rydberg_site_nearest_storage_site")
  #print(self.Rydberg_site_nearest_storage_site)

      def distance(self, idx1, r1, c1, idx2, r2, c2):
          p1 = self.exact_SLM_location(idx1, r1, c1)
          p2 = self.exact_SLM_location(idx2, r2, c2)
          return math.dist(p1, p2) # Euclidean distance

      def nearest_storage_site(self, idx, r, c):
          slm = self.dict_SLM[idx]
          slm_idx = self.entanglement_zone[slm.entanglement_id][0]
          slm = self.dict_SLM[slm_idx]
          nearest_slm_half_r = slm.n_r // 2
          if r < nearest_slm_half_r:
              return self.Rydberg_site_nearest_storage_site[slm_idx][0][c]
          else:
              return self.Rydberg_site_nearest_storage_site[slm_idx][1][c]

      def nearest_entanglement_site(self, idx, r, c):
  #return the nearest Rydberg site for a qubit in the storage zone
          return self.storage_site_nearest_Rydberg_site[idx][r][c]

      def nearest_entanglement_site_distance(self, idx, r, c):
  #return the distance nearest Rydberg site for a qubit in the           \
              storage zone
          return self.storage_site_nearest_Rydberg_site_dis[idx][r][c]

      def nearest_entanglement_site(self, idx1, r1, c1, idx2, r2, c2):
  #return the nearest Rydberg site for two qubit in the storage zone
  #based on the position of two qubits
          storage_site1 = self.exact_SLM_location(idx1, r1, c1)
          storage_site2 = self.exact_SLM_location(idx2, r2, c2)
          site1 = self.storage_site_nearest_Rydberg_site[idx1][r1][c1]
          site2 = self.storage_site_nearest_Rydberg_site[idx2][r2][c2]
  #the nearest zone for both qubits are in the same entanglement zone
          if site1 == site2:
              return [site1]
          elif site1[0] == site2[0]:
              near_x = (storage_site1[0] + storage_site1[1]) // 2
              middle_site_c = (site1[2] + site2[2])// 2
              middle_site_x = self.exact_SLM_location(site1[0], site1[1],
  middle_site_c)[0] near_site_idx = middle_site_c #near_site_dis =
  abs(middle_site_x - near_x) #slm = self.dict_SLM[site1[0]] #next_dis =
  abs(middle_site_x + slm.seperation[0] - near_x) #if near_site_idx < slm.n_c -
  1 and next_dis < near_site_dis: #near_site_idx += 1 #elif near_site_idx > 0:
  #next_dis = abs(middle_site_x - slm.seperation[0] - near_x)
  #if next_dis < near_site_dis:
  #near_site_idx -= 1

              return [(site1[0], site1[1], near_site_idx)]
          else:
              return [site1, site2]
              slm1 = self.dict_SLM[site1[0]]
              slm2 = self.dict_SLM[site2[0]]
              row_y_1 = slm1.location[1] + site1[1] * slm1.site_seperation[1]
              row_y_2 = slm2.location[1] + site2[1] * slm2.site_seperation[1]
              diff_y_1 = abs(row_y_1 - storage_site1[1]) + abs(row_y_1 -
  storage_site2[1]) diff_y_2 = abs(row_y_2 - storage_site1[1]) + abs(row_y_2 -
  storage_site2[1]) if diff_y_1 < diff_y_2: return (site1[0], site1[1],
  (site1[2] + site2[2])// 2) else: return (site2[0], site2[1], (site1[2] +
  site2[2])// 2)

      def nearest_entanglement_site_dis(self, idx1, r1, c1, idx2, r2, c2):
  #return the sum of the distance to move two qubits to one rydberg site
          storage_site1 = self.exact_SLM_location(idx1, r1, c1)
          storage_site2 = self.exact_SLM_location(idx2, r2, c2)
          list_site = self.nearest_entanglement_site(idx1, r1, c1, idx2, r2, c2)
          dis = sys.maxsize
          for site in list_site:
              exact_site = self.exact_SLM_location(site[0], site[1], site[2])
  #return math.dist(storage_site1, exact_site) +                         \
              math.dist(storage_site2, exact_site)
              if r1 == r2 and idx1 == idx2:
                  dis = min(max(math.dist(storage_site1, exact_site),
  math.dist(storage_site2, exact_site)), dis) else: dis = min(
  math.dist(storage_site1, exact_site) + math.dist(storage_site2, exact_site),
  dis ) return dis


      def movement_duration(self, x1, y1, x2, y2):
  #d / t ^ 2 = a = 2750m / s
          a = 0.00275
          d = math.dist((x1, y1), (x2, y2))
  #d = 15
          t = math.sqrt(d/a)
          return t
          */
};
} // namespace na