#pragma once

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <istream>
#include <memory>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace na {
/// An 2D-Array of AOD traps
struct AOD {
  std::size_t id = 0;
  std::size_t siteSeparation = 0;
  std::size_t nRows = 0;
  std::size_t nCols = 0;

  explicit AOD(nlohmann::json aodSpec);
};

/// An 2D-array of SLM traps
struct SLM {
  std::size_t id = 0; ///< SLM id, used only in output
  /// separation of individual sites in x and y direction
  std::pair<std::size_t, std::size_t> siteSeparation{0, 0};
  std::size_t nRows = 0; ///< number of rows
  std::size_t nCols = 0; ///< number of columns
  /// x,y-coordinate of the left uppermost SLM
  std::pair<std::size_t, std::size_t> location{0, 0};
  /// if the SLM is used in entanglement zone, a pointer to all entanglement
  /// SLMs in the same group
  std::vector<std::unique_ptr<SLM>>* entanglementZone = nullptr;
  std::optional<std::size_t> entanglementId = std::nullopt;

  explicit SLM(nlohmann::json slmSpec);
  explicit SLM(nlohmann::json slmSpec,
               decltype(entanglementZone) entanglementZone,
               std::size_t entanglementId);
  [[nodiscard]] auto isStorage() const -> bool {
    return entanglementZone == nullptr;
  }
  [[nodiscard]] auto isEntanglement() const -> bool { return !isStorage(); }
};
} // namespace na

namespace std {

template <> struct hash<tuple<const na::SLM*, size_t, size_t>> {
  size_t
  operator()(const tuple<const na::SLM*, size_t, size_t>& t) const noexcept {
    const auto& [slm, a, b] = t;
    const size_t h1 = hash<const na::SLM*>{}(slm);
    const size_t h2 = hash<size_t>{}(a);
    const size_t h3 = hash<size_t>{}(b);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

template <> struct hash<tuple<const na::SLM*, size_t, const na::SLM*, size_t>> {
  size_t operator()(const tuple<const na::SLM*, size_t, const na::SLM*, size_t>&
                        t) const noexcept {
    const auto& [slm1, x1, slm2, x2] = t;
    const size_t h1 = hash<const na::SLM*>{}(slm1);
    const size_t h2 = hash<size_t>{}(x1);
    const size_t h3 = hash<const na::SLM*>{}(slm2);
    const size_t h4 = hash<size_t>{}(x2);
    return h1 ^ (h3 << 1) ^ (h2 << 2) ^ (h4 << 3);
  }
};
} // namespace std

namespace na {

/// Class to define zone architecture
struct Architecture {
  std::string name;
  std::vector<std::unique_ptr<SLM>> storageZones;
  std::vector<std::vector<std::unique_ptr<SLM>>> entanglementZones;
  std::vector<std::unique_ptr<AOD>> aods;
  double timeAtomTransfer = 15; ///< µs
  double timeRydberg = 0.36;    ///< µs
  double time1QGate = 0.625;    ///< µs
  std::size_t archRangeMinX = 0;
  std::size_t archRangeMaxX = 0;
  std::size_t archRangeMinY = 0;
  std::size_t archRangeMaxY = 0;
  std::vector<std::size_t> rydbergRangeMinX;
  std::vector<std::size_t> rydbergRangeMaxX;
  std::vector<std::size_t> rydbergRangeMinY;
  std::vector<std::size_t> rydbergRangeMaxY;
  /// A map from an entanglement site to the nearest storage sites in ascending
  /// order by their distance.
  /// To get the nearest storage site expressed as a triple of
  /// (SLM*, row, column), use the following code:
  /// @code
  /// entanglementToNearestStorageSite[slm][0/1][column][0];
  /// @endcode
  /// The second index denotes the slm in a pair of two slms forming an
  /// entanglement zone.
  /// The third index denotes the i-th nearest storage site.
  /// @see storageToNearestEntanglementSite
  std::unordered_map<const SLM*, std::vector<std::vector<std::tuple<
                                     const SLM*, std::size_t, std::size_t>>>>
      entanglementToNearestStorageSite;
  /// A map from a storage site to the nearest Rydberg sites.
  /// @see entanglementToNearestStorageSite
  std::unordered_map<
      const SLM*,
      std::vector<std::vector<std::unordered_map<
          const SLM*, std::vector<std::vector<
                          std::tuple<const SLM*, std::size_t, std::size_t>>>>>>>
      storageToNearestEntanglementSite;

  Architecture() = default;
  explicit Architecture(const std::string& filename)
      : Architecture(std::filesystem::path(filename)) {}
  explicit Architecture(const std::filesystem::path& filepath)
      : Architecture(std::ifstream(filepath)) {}
  explicit Architecture(std::istream& is) : Architecture(std::move(is)) {}
  explicit Architecture(std::istream&& is) {
    load(std::move(is));
    preprocessing();
  }
  explicit Architecture(nlohmann::json& json) : Architecture(std::move(json)) {}
  explicit Architecture(nlohmann::json&& json) {
    load(std::move(json));
    preprocessing();
  }
  auto load(const std::string& filename) -> void {
    load(std::filesystem::path(filename));
  }
  auto load(const std::filesystem::path& filepath) -> void {
    load(std::ifstream(filepath));
  }
  auto load(std::istream& is) -> void { load(std::move(is)); }
  auto load(std::istream&& is) -> void {
    nlohmann::json architectureSpec{};
    std::move(is) >> architectureSpec;
    load(std::move(architectureSpec));
  }
  auto load(const nlohmann::json& architectureSpec) -> void {
    load(std::move(architectureSpec));
  }
  auto load(const nlohmann::json&& architectureSpec) -> void;
  auto exportNAVizMachine() const -> std::string;
  auto exportNAVizMachine(std::ostream&& os) const -> void {
    os << exportNAVizMachine();
  }
  auto exportNAVizMachine(std::ostream& os) const -> void {
    exportNAVizMachine(std::move(os));
  }
  auto exportNAVizMachine(const std::filesystem::path& path) const -> void {
    exportNAVizMachine(std::ofstream(path));
  }
  auto exportNAVizMachine(const std::string& filename) const -> void {
    exportNAVizMachine(std::filesystem::path(filename));
  }
  //===--------------------------------------------------------------------===//
  /// Check if the given position is a valid SLM position, i.e., whether the
  /// given row and column are within the range of the SLM.
  auto isValidSlmPosition(const SLM& slm, std::size_t r, std::size_t c) const
      -> bool;
  auto
  /// @see isValidSlmPosition
  isValidSlmPosition(const std::tuple<const SLM*, const std::size_t,
                                      const std::size_t>& t) const -> bool {
    return isValidSlmPosition(*std::get<0>(t), std::get<1>(t), std::get<2>(t));
  }
  //===--------------------------------------------------------------------===//
  /// Compute the exact location of the SLM site given the row and column
  /// indices expressed in coordinates in the global coordinate system.
  auto exactSlmLocation(const SLM& slm, std::size_t r, std::size_t c) const
      -> std::pair<std::size_t, std::size_t>;
  /// @see exactSlmLocation
  auto exactSlmLocation(const std::tuple<const SLM*, const std::size_t,
                                         const std::size_t>& t) const
      -> std::pair<std::size_t, std::size_t> {
    return exactSlmLocation(*std::get<0>(t), std::get<1>(t), std::get<2>(t));
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
  auto distance(
      const std::tuple<const SLM*, const std::size_t, const std::size_t>& t1,
      const std::tuple<const SLM*, const std::size_t, const std::size_t>& t2)
      const -> double {
    return distance(*std::get<0>(t1), std::get<1>(t1), std::get<2>(t1),
                    *std::get<0>(t2), std::get<1>(t2), std::get<2>(t2));
  }
  //===--------------------------------------------------------------------===//
  /// return the nearest storage site for an entanglement site
  auto nearestStorageSite(const SLM& slm, std::size_t r, std::size_t c) const
      -> const std::tuple<const SLM*, std::size_t, std::size_t>&;
  /// @see nearestStorageSite
  auto nearestStorageSite(const std::tuple<const SLM*, const std::size_t,
                                           const std::size_t>& t) const
      -> const std::tuple<const SLM*, std::size_t, std::size_t>& {
    return nearestStorageSite(*std::get<0>(t), std::get<1>(t), std::get<2>(t));
  }
  //===--------------------------------------------------------------------===//
  /// return the nearest Rydberg site for two qubit in the storage zone
  /// based on the position of two qubits
  auto nearestEntanglementSite(const SLM& idx1, std::size_t r1, std::size_t c1,
                               const SLM& idx2, std::size_t r2,
                               std::size_t c2) const
      -> const std::tuple<const SLM*, std::size_t, std::size_t>&;
  /// @see nearestEntanglementSite
  auto nearestEntanglementSite(
      const std::tuple<const SLM*, std::size_t, std::size_t>& t1,
      const std::tuple<const SLM*, std::size_t, std::size_t>& t2) const
      -> const std::tuple<const SLM*, std::size_t, std::size_t>& {
    return nearestEntanglementSite(*std::get<0>(t1), std::get<1>(t1),
                                   std::get<2>(t1), *std::get<0>(t2),
                                   std::get<1>(t2), std::get<2>(t2));
  }
  //===--------------------------------------------------------------------===//
  /// return the maximum/sum of the distance to move two qubits to one rydberg
  /// site. If the two qubits are in the same row, i.e., can be picked up
  /// simultaneously, the maximum distance is returned. Otherwise, the
  /// sum of the distances is returned.
  auto nearestEntanglementSiteDistance(const SLM& slm1, std::size_t r1,
                                       std::size_t c1, const SLM& slm2,
                                       std::size_t r2, std::size_t c2) const
      -> double;
  /// @see nearestEntanglementSiteDistance
  auto nearestEntanglementSiteDistance(
      const std::tuple<const SLM*, std::size_t, std::size_t>& t1,
      const std::tuple<const SLM*, std::size_t, std::size_t>& t2) const
      -> double {
    return nearestEntanglementSiteDistance(*std::get<0>(t1), std::get<1>(t1),
                                           std::get<2>(t1), *std::get<0>(t2),
                                           std::get<1>(t2), std::get<2>(t2));
  }
  //===--------------------------------------------------------------------===//
  /// Returns the time to move from one location to another location
  static auto movementDuration(std::size_t x1, std::size_t y1, std::size_t x2,
                               std::size_t y2) -> double;
  /// @see movementDuration
  auto movementDuration(const std::pair<std::size_t, std::size_t>& p1,
                        const std::pair<std::size_t, std::size_t>& p2) const
      -> double {
    return movementDuration(p1.first, p1.second, p2.first, p2.second);
  }
  //===--------------------------------------------------------------------===//
  /// Returns the other site of a pair of entanglement sites
  auto otherEntanglementSite(const SLM& slm, std::size_t r, std::size_t c) const
      -> std::tuple<const SLM*, std::size_t, std::size_t>;
  /// @see otherEntanglementSite
  auto otherEntanglementSite(
      const std::tuple<const SLM*, std::size_t, std::size_t>& t) const
      -> std::tuple<const SLM*, std::size_t, std::size_t> {
    return otherEntanglementSite(*std::get<0>(t), std::get<1>(t),
                                 std::get<2>(t));
  }
};
} // namespace na
