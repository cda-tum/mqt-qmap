#pragma once

#include "Definitions.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <istream>
#include <memory>
#include <nlohmann/json.hpp>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
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
  const std::array<SLM, 2>* entanglementZone_ = nullptr;
  /// only used for printing
  std::optional<std::size_t> entanglementId_ = std::nullopt;

  explicit SLM(nlohmann::json slmSpec);
  [[nodiscard]] auto isEntanglement() const -> bool {
    return entanglementZone_ != nullptr;
  }
  [[nodiscard]] auto isStorage() const -> bool { return !isEntanglement(); }
  [[nodiscard]] auto operator==(const SLM& other) const -> bool;
};
} // namespace na
template <class F, class S> struct std::hash<std::pair<F, S>> {
  size_t operator()(const std::pair<F, S>& p) const noexcept {
    const auto h1 =
        std::hash<std::remove_cv_t<std::remove_reference_t<F>>>{}(p.first);
    const auto h2 =
        std::hash<std::remove_cv_t<std::remove_reference_t<S>>>{}(p.second);
    return qc::combineHash(h1, h2);
  }
};
template <class S, class T, class U> struct std::hash<std::tuple<S, T, U>> {
  size_t operator()(const std::tuple<S, T, U>& t) const noexcept {
    const auto h1 = std::hash<std::remove_cv_t<std::remove_reference_t<S>>>{}(
        std::get<0>(t));
    const auto h2 = std::hash<std::remove_cv_t<std::remove_reference_t<T>>>{}(
        std::get<1>(t));
    const auto h3 = std::hash<std::remove_cv_t<std::remove_reference_t<U>>>{}(
        std::get<2>(t));
    return qc::combineHash(qc::combineHash(h1, h2), h3);
  }
};
template <> struct std::hash<na::SLM> {
  size_t operator()(const na::SLM& slm) const noexcept {
    return std::hash<std::pair<size_t, size_t>>{}(slm.location);
  }
};
template <class T> struct std::hash<std::array<T, 2>> {
  size_t operator()(const std::array<T, 2>& p) const noexcept {
    const auto h1 =
        std::hash<std::remove_cv_t<std::remove_reference_t<T>>>{}(p.front());
    const auto h2 =
        std::hash<std::remove_cv_t<std::remove_reference_t<T>>>{}(p.back());
    return qc::combineHash(h1, h2);
  }
};

namespace na {

/// Class to define zone architecture
struct Architecture {
  std::string name;
  std::vector<std::unique_ptr<SLM>> storageZones;
  std::vector<std::unique_ptr<std::array<SLM, 2>>> entanglementZones;
  std::vector<std::unique_ptr<AOD>> aods;
  struct OperationDurations {
    double timeAtomTransfer = 15; ///< µs
    double timeRydberg = 0.36;    ///< µs
    double time1QGate = 0.625;    ///< µs
  };
  std::optional<OperationDurations> operationDurations;
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
  std::unordered_map<
      std::reference_wrapper<const SLM>,
      std::vector<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                         std::size_t, std::size_t>>>,
      std::hash<SLM>, std::equal_to<SLM>>
      entanglementToNearestStorageSite;
  /// A map from a storage site to the nearest Rydberg sites.
  /// @see entanglementToNearestStorageSite
  std::unordered_map<
      std::reference_wrapper<const SLM>,
      std::vector<std::vector<std::unordered_map<
          std::reference_wrapper<const SLM>,
          std::vector<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                             std::size_t, std::size_t>>>,
          std::hash<SLM>, std::equal_to<SLM>>>>,
      std::hash<SLM>, std::equal_to<SLM>>
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
  // Delete copy constructor and copy assignment operator
  Architecture(const Architecture&) = delete;
  Architecture& operator=(const Architecture&) = delete;
  // Default move constructor and move assignment operator
  Architecture(Architecture&&) noexcept = default;
  Architecture& operator=(Architecture&&) noexcept = default;
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
  [[nodiscard]] auto isValidSlmPosition(const SLM& slm, std::size_t r,
                                        std::size_t c) const -> bool;
  /// Compute the exact location of the SLM site given the row and column
  /// indices expressed in coordinates in the global coordinate system.
  [[nodiscard]] auto exactSlmLocation(const SLM& slm, std::size_t r,
                                      std::size_t c) const
      -> std::pair<std::size_t, std::size_t>;
  /**
   * In the loop, we will calculate a lower bound of the distance
   * between the entanglement site and a storage SLM. Any site in the
   * storage SLM will have at least this distance to the entanglement
   * site. This distance will be the variable @c minimalDistance.
   * Among all storage SLMs, we will find the one that has the minimum
   * distance to the entanglement site, the @c minimumDistance.
   */
  [[nodiscard]] auto findNearestStorageSLM(size_t x, size_t y) const
      -> const SLM&;
  /**
   * In the loop, we will calculate a lower bound of the distance
   * between the entanglement site and a storage SLM. Any site in
   * the storage SLM will have at least this distance to the
   * entanglement site. This distance will be the variable @c
   * minimalDistance. Among all storage SLMs, we will find the one
   * that has the minimum distance to the entanglement site, the @c
   * minimumDistance.
   */
  [[nodiscard]] auto findNearestEntanglementSLM(size_t x, size_t y,
                                                size_t otherX,
                                                size_t otherY) const
      -> const SLM&;
  /// Compute the site region for entanglement zone and the nearest Rydberg site
  /// for each storage site.
  /// @note We assume we only have one storage zone or one entanglement zone per
  /// row.
  auto preprocessing() -> void;
  /// Compute the distance between two specific SLM sites
  auto distance(const SLM& idx1, std::size_t r1, std::size_t c1,
                const SLM& idx2, std::size_t r2, std::size_t c2) const
      -> double;
  /// return the nearest storage site for an entanglement site
  auto nearestStorageSite(const SLM& slm, std::size_t r, std::size_t c) const
      -> const
      std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t>&;
  /// return the nearest Rydberg site for two qubit in the storage zone
  /// based on the position of two qubits
  auto nearestEntanglementSite(const SLM& idx1, std::size_t r1, std::size_t c1,
                               const SLM& idx2, std::size_t r2,
                               std::size_t c2) const -> const
      std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t>&;
  /// return the maximum/sum of the distance to move two qubits to one rydberg
  /// site. If the two qubits are in the same row, i.e., can be picked up
  /// simultaneously, the maximum distance is returned. Otherwise, the
  /// sum of the distances is returned.
  auto nearestEntanglementSiteDistance(const SLM& slm1, std::size_t r1,
                                       std::size_t c1, const SLM& slm2,
                                       std::size_t r2, std::size_t c2) const
      -> double;
  /// Returns the time to move from one location to another location
  static auto movementDuration(std::size_t x1, std::size_t y1, std::size_t x2,
                               std::size_t y2) -> double;
  /// Returns the other site of a pair of entanglement sites
  auto otherEntanglementSite(const SLM& slm, std::size_t r, std::size_t c) const
      -> std::tuple<std::reference_wrapper<const SLM>, std::size_t,
                    std::size_t>;
};
} // namespace na
