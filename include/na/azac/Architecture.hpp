#pragma once

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <istream>
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace na {
/// An 2D-Array of AOD traps
struct AOD {
  std::size_t id = 0;
  std::pair<std::size_t, std::size_t> site_separation{0, 0};
  std::size_t n_r = 0;
  std::size_t n_c = 0;

  explicit AOD(nlohmann::json aodSpec);
};

/// An 2D-array of SLM traps
struct SLM {
  std::size_t id = 0; ///< SLM id
  /// separation of individual sites in x and y direction
  std::pair<std::size_t, std::size_t> site_separation{0, 0};
  std::size_t n_r = 0; ///< number of rows
  std::size_t n_c = 0; ///< number of columns
  /// x,y-coordinate of the left uppermost SLM
  std::pair<std::size_t, std::size_t> location{0, 0};
  /// if the SLM is used in entanglement zone, a pointer to all entanglement
  /// SLMs in the same group
  std::vector<std::unique_ptr<SLM>>* entanglement_id = nullptr;

  explicit SLM(nlohmann::json slmSpec,
               decltype(entanglement_id) entanglementId = nullptr);
};
} // namespace na

template<>
struct std::hash<const na::SLM* const>
{
  std::size_t operator()(const na::SLM* const slm) const noexcept
  {
    const std::size_t h1 = std::hash<std::size_t>{}(slm->location.first);
    const std::size_t h2 = std::hash<std::size_t>{}(slm->location.second);
    return h1 ^ (h2 << 1);
  }
};

namespace na {

/// Class to define zone architecture
struct Architecture {
  std::string name;
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
  /// A map from an entanglement site to the nearest storage site.
  /// To get the storage site that is expressed as a triple of
  /// (SLM*, row, column), use the following code:
  /// @code
  /// Rydberg_site_nearest_storage_site[slm][0/1][column];
  /// @endcode
  /// The second index denotes the slm in a pair of two slms forming an
  /// entanglement zone.
  /// @see storage_site_nearest_Rydberg_site_dis
  std::unordered_map<
      const SLM* const,
      std::vector<
          std::vector<std::tuple<const SLM* const, std::size_t, std::size_t>>>>
      entanglementToNearestStorageSite;
  /// A map from a storage site to the nearest Rydberg site.
  /// @see Rydberg_site_nearest_storage_site
  std::unordered_map<
      const SLM* const,
      std::vector<
          std::vector<std::tuple<const SLM* const, std::size_t, std::size_t>>>>
      storageToNearestEntanglementSite;
  /// A map from a storage site to the distance to the nearest Rydberg site.
  /// @see Rydberg_site_nearest_storage_site
  std::unordered_map<const SLM* const, std::vector<std::vector<double>>>
      storageToNearestEntanglementSiteDistance;

  Architecture() = default;
  explicit Architecture(const std::string& filename)
      : Architecture(std::filesystem::path(filename)) {}
  explicit Architecture(const std::filesystem::path& filepath)
      : Architecture(std::ifstream(filepath)) {}
  explicit Architecture(std::istream& is) : Architecture(std::move(is)) {}
  explicit Architecture(std::istream&& is) { load(std::move(is)); }
  auto load(const std::string& filename) -> void {
    load(std::filesystem::path(filename));
  }
  auto load(const std::filesystem::path& filepath) -> void {
    load(std::ifstream(filepath));
  }
  auto load(std::istream& is) -> void { load(std::move(is)); }
  auto load(std::istream&& is) -> void;
  //===--------------------------------------------------------------------===//
  /// @see is_valid_SLM_position
  /// Check if the given position is a valid SLM position, i.e., whether the
  /// given row and column are within the range of the SLM.
  auto is_valid_SLM_position(const SLM& slm, std::size_t r, std::size_t c) const
      -> bool;
  auto
  is_valid_SLM_position(const std::tuple<const SLM* const, const std::size_t,
                                         const std::size_t>& t) const -> bool {
    return is_valid_SLM_position(std::get<0>(t), std::get<1>(t),
                                 std::get<2>(t));
  }
  /// @see is_valid_SLM_position
  auto is_valid_SLM_position(const SLM* const slm, const std::size_t r,
                             const std::size_t c) const -> bool {
    return is_valid_SLM_position(*slm, r, c);
  }
  //===--------------------------------------------------------------------===//
  /// Compute the exact location of the SLM site given the row and column
  /// indices expressed in coordinates in the global coordinate system.
  auto exact_SLM_location(const SLM& slm, std::size_t r, std::size_t c) const
      -> std::pair<std::size_t, std::size_t>;
  /// @see exact_SLM_location
  auto exact_SLM_location(const std::tuple<const SLM* const, const std::size_t,
                                           const std::size_t>& t) const
      -> std::pair<std::size_t, std::size_t> {
    return exact_SLM_location(std::get<0>(t), std::get<1>(t), std::get<2>(t));
  }
  /// @see exact_SLM_location
  auto exact_SLM_location(const SLM* const slm, const std::size_t r,
                          const std::size_t c) const
      -> std::pair<std::size_t, std::size_t> {
    return exact_SLM_location(*slm, r, c);
  }
  //===--------------------------------------------------------------------===//
  /// Compute the site region for entanglement zone and the nearest Rydberg site
  /// for each storage site.
  /// @note We assume we only have one storage zone or one entanglement zone per
  /// row.
  auto preprocessing() -> void;
  //===--------------------------------------------------------------------===//
  /// Compute the distance between two specific SLM sites
  auto distance(const SLM& idx1, std::size_t r1, std::size_t c1,
                const SLM& idx2, std::size_t r2, std::size_t c2) const
      -> double;
  /// @see distance
  auto distance(const std::tuple<const SLM* const, const std::size_t,
                                 const std::size_t>& t1,
                const std::tuple<const SLM* const, const std::size_t,
                                 const std::size_t>& t2) const -> double {
    return distance(std::get<0>(t1), std::get<1>(t1), std::get<2>(t1),
                    std::get<0>(t2), std::get<1>(t2), std::get<2>(t2));
  }
  /// @see distance
  auto distance(const SLM* const idx1, const std::size_t r1,
                const std::size_t c1, const SLM* const idx2,
                const std::size_t r2, const std::size_t c2) const -> double {
    return distance(*idx1, r1, c1, *idx2, r2, c2);
  }
  //===--------------------------------------------------------------------===//
  /// return the nearest storage site for an entanglement site
  auto nearest_storage_site(const SLM& slm, std::size_t r, std::size_t c) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t>;
  /// @see nearest_storage_site
  auto nearest_storage_site(
      const std::tuple<const SLM* const, const std::size_t, const std::size_t>&
          t) const -> std::tuple<const SLM* const, std::size_t, std::size_t> {
    return nearest_storage_site(std::get<0>(t), std::get<1>(t), std::get<2>(t));
  }
  /// @see nearest_storage_site
  auto nearest_storage_site(const SLM* const slm, const std::size_t r,
                            const std::size_t c) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t> {
    return nearest_storage_site(*slm, r, c);
  }
  //===--------------------------------------------------------------------===//
  /// return the nearest Rydberg site for a qubit in the storage zone
  auto nearest_entanglement_site(const SLM* idx, std::size_t r,
                                 std::size_t c) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t>;
  /// @see nearest_entanglement_site
  auto nearest_entanglement_site(
      const std::tuple<const SLM* const, std::size_t, std::size_t>& t) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t> {
    return nearest_entanglement_site(std::get<0>(t), std::get<1>(t),
                                     std::get<2>(t));
  }
  /// @see nearest_entanglement_site
  auto nearest_entanglement_site(const SLM& slm, const std::size_t r,
                                 const std::size_t c) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t> {
    return nearest_entanglement_site(&slm, r, c);
  }
  //===--------------------------------------------------------------------===//
  /// return the distance nearest Rydberg site for a qubit in the storage zone
  auto nearest_entanglement_site_distance(const SLM* idx, std::size_t r,
                                          std::size_t c) const -> double;
  /// @see nearest_entanglement_site_distance
  auto nearest_entanglement_site_distance(
      const std::tuple<const SLM* const, std::size_t, std::size_t>& t1) const
      -> double {
    return nearest_entanglement_site_distance(std::get<0>(t1), std::get<1>(t1),
                                              std::get<2>(t1));
  }
  /// @see nearest_entanglement_site_distance
  auto nearest_entanglement_site_distance(const SLM& slm, const std::size_t r,
                                          const std::size_t c) const -> double {
    return nearest_entanglement_site_distance(&slm, r, c);
  }
  //===--------------------------------------------------------------------===//
  /// return the nearest Rydberg site for two qubit in the storage zone
  /// based on the position of two qubits
  auto nearest_entanglement_site(const SLM* idx1, std::size_t r1,
                                 std::size_t c1, const SLM* idx2,
                                 std::size_t r2, std::size_t c2) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t>;
  /// @see nearest_entanglement_site
  auto nearest_entanglement_site(
      const std::tuple<const SLM* const, std::size_t, std::size_t>& t1,
      const std::tuple<const SLM* const, std::size_t, std::size_t>& t2) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t> {
    return nearest_entanglement_site(std::get<0>(t1), std::get<1>(t1),
                                     std::get<2>(t1), std::get<0>(t2),
                                     std::get<1>(t2), std::get<2>(t2));
  }
  /// @see nearest_entanglement_site
  auto nearest_entanglement_site(const SLM& idx1, const std::size_t r1,
                                 const std::size_t c1, const SLM& idx2,
                                 const std::size_t r2,
                                 const std::size_t c2) const
      -> std::tuple<const SLM* const, std::size_t, std::size_t> {
    return nearest_entanglement_site(&idx1, r1, c1, &idx2, r2, c2);
  }
  //===--------------------------------------------------------------------===//
  /// return the maximum/sum of the distance to move two qubits to one rydberg
  /// site. If the two qubits are in the same row, i.e., can be picked up
  /// simultaneously, the maximum distance is returned. Otherwise, the
  /// sum of the distances is returned.
  auto nearest_entanglement_site_distance(const SLM* slm1, std::size_t r1,
                                          std::size_t c1, const SLM* slm2,
                                          std::size_t r2, std::size_t c2) const
      -> double;
  /// @see nearest_entanglement_site_distance
  auto nearest_entanglement_site_distance(
      const std::tuple<const SLM* const, std::size_t, std::size_t>& t1,
      const std::tuple<const SLM* const, std::size_t, std::size_t>& t2) const
      -> double {
    return nearest_entanglement_site_distance(std::get<0>(t1), std::get<1>(t1),
                                              std::get<2>(t1), std::get<0>(t2),
                                              std::get<1>(t2), std::get<2>(t2));
  }
  /// @see nearest_entanglement_site_distance
  auto nearest_entanglement_site_distance(const SLM& slm1, const std::size_t r1,
                                          const std::size_t c1, const SLM& slm2,
                                          const std::size_t r2,
                                          const std::size_t c2) const
      -> double {
    return nearest_entanglement_site_distance(&slm1, r1, c1, &slm2, r2, c2);
  }
  //===--------------------------------------------------------------------===//
  /// Returns the time to move from one location to another location
  static auto movement_duration(std::size_t x1, std::size_t y1, std::size_t x2,
                                std::size_t y2) -> double;
  /// @see movement_duration
  auto movement_duration(const std::pair<std::size_t, std::size_t>& p1,
                         const std::pair<std::size_t, std::size_t>& p2) const
      -> double {
    return movement_duration(p1.first, p1.second, p2.first, p2.second);
  }
};
} // namespace na

