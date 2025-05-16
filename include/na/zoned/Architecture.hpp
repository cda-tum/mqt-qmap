/*
 * Copyright (c) 2023 - 2025 Chair for Design Automation, TUM
 * Copyright (c) 2025 Munich Quantum Software Company GmbH
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "ir/Definitions.hpp"

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

namespace na::zoned {
/// An 2D-Array of AOD traps
struct AOD {
  std::size_t id = 0;
  std::size_t siteSeparation = 0;
  std::size_t nRows = 0;
  std::size_t nCols = 0;

  /// Creates an AOD from a JSON specification.
  explicit AOD(nlohmann::json aodSpec);
};

/// An 2D-array of SLM traps.
struct SLM {
  std::size_t id = 0; ///< SLM id, used only in output
  /// separation of individual sites in x and y direction.
  std::pair<std::size_t, std::size_t> siteSeparation{0, 0};
  std::size_t nRows = 0; ///< number of rows
  std::size_t nCols = 0; ///< number of columns
  /// x,y-coordinate of the left uppermost SLM.
  std::pair<std::size_t, std::size_t> location{0, 0};
  /// if the SLM is used in entanglement zone, a pointer to all entanglement
  /// SLMs in the same group.
  const std::array<SLM, 2>* entanglementZone_ = nullptr;
  /// only used for printing.
  std::optional<std::size_t> entanglementId_ = std::nullopt;
  /// Creates an SLM array from a JSON specification.
  explicit SLM(nlohmann::json slmSpec);
  /// @returns true if the SLM is part of an entanglement zone.
  [[nodiscard]] auto isEntanglement() const -> bool {
    return entanglementZone_ != nullptr;
  }
  /// @returns true, if the SLM is part of a storage zone.
  [[nodiscard]] auto isStorage() const -> bool { return !isEntanglement(); }
  /// @returns true, if both SLMs are equal, i.e., they have the same
  /// location and dimensions.
  [[nodiscard]] auto operator==(const SLM& other) const -> bool;
};
/// An element of type Site identifies a concrete site in an SLM array.
using Site =
    std::tuple<std::reference_wrapper<const SLM>, std::size_t, std::size_t>;
/// An unordered map from an SLM to a value of type V.
/// @tparam V the type of the value
template <class V>
using SLMMap = std::unordered_map<std::reference_wrapper<const SLM>, V,
                                  std::hash<SLM>, std::equal_to<SLM>>;
} // namespace na::zoned
/// Hash function for pairs where both unqualified types have a hash function.
template <class F, class S> struct std::hash<std::pair<F, S>> {
  size_t operator()(const std::pair<F, S>& p) const noexcept {
    const auto h1 =
        std::hash<std::remove_cv_t<std::remove_reference_t<F>>>{}(p.first);
    const auto h2 =
        std::hash<std::remove_cv_t<std::remove_reference_t<S>>>{}(p.second);
    return qc::combineHash(h1, h2);
  }
};
/// Hash function for 3-tuples where all unqualified types have a hash function.
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
/// Hash function for an SLM based on its location.
template <> struct std::hash<na::zoned::SLM> {
  size_t operator()(const na::zoned::SLM& slm) const noexcept {
    return std::hash<std::pair<size_t, size_t>>{}(slm.location);
  }
};
/// An array of length 2 where the unqualified type has a hash function.
template <class T> struct std::hash<std::array<T, 2>> {
  size_t operator()(const std::array<T, 2>& p) const noexcept {
    const auto h1 =
        std::hash<std::remove_cv_t<std::remove_reference_t<T>>>{}(p.front());
    const auto h2 =
        std::hash<std::remove_cv_t<std::remove_reference_t<T>>>{}(p.back());
    return qc::combineHash(h1, h2);
  }
};

namespace na::zoned {

/// Class to define zone architecture.
struct Architecture {
  std::string name; ///< name of the architecture
  /// All storage zones of the architecture. The objects are owned by the
  /// Architecture class. Outside the class, the SLMs are only referenced.
  std::vector<std::unique_ptr<SLM>> storageZones;
  /// All entanglement zones of the architecture. Each entanglement zone
  /// consists of two SLMs. The objects are owned by the Architecture class.
  /// Outside the class, the SLMs are only referenced.
  std::vector<std::unique_ptr<std::array<SLM, 2>>> entanglementZones;
  /// All AODs of the architecture. The objects are owned by the Architecture
  /// class. Outside the class, the AODs are only referenced.
  std::vector<std::unique_ptr<AOD>> aods;
  /// A struct to define the operation durations.
  struct OperationDurations {
    double timeAtomTransfer = 15;       ///< µs
    double timeRydbergGate = 0.36;      ///< µs
    double timeSingleQubitGate = 0.625; ///< µs
  };
  /// Operation durations.
  std::optional<OperationDurations> operationDurations;
  /// A struct to define the operation fidelities.
  struct OperationFidelities {
    double fidelityRydbergGate = 0.995;
    double fidelitySingleQubitGate = 0.9997;
    double fidelityAtomTransfer = 0.999;
  };
  /// Operation fidelities.
  std::optional<OperationFidelities> operationFidelities;
  std::optional<double> qubitT1; ///< T1 time of the qubit in µs
  /// Minimum X coordinates of the different Rydberg zones, i.e., where the
  /// Rydberg laser can affect the atoms.
  std::vector<std::size_t> rydbergRangeMinX;
  /// Maximum X coordinates of the different Rydberg zones, i.e., where the
  /// Rydberg laser can affect the atoms.
  std::vector<std::size_t> rydbergRangeMaxX;
  /// Minimum Y coordinates of the different Rydberg zones, i.e., where the
  /// Rydberg laser can affect the atoms.
  std::vector<std::size_t> rydbergRangeMinY;
  /// Maximum Y coordinates of the different Rydberg zones, i.e., where the
  /// Rydberg laser can affect the atoms.
  std::vector<std::size_t> rydbergRangeMaxY;
  /// A map from an entanglement site to the nearest storage site.
  /// @see storageToNearestEntanglementSite
  SLMMap<std::vector<std::vector<Site>>> entanglementToNearestStorageSite;
  /// A map from a pair of storage sites to the nearest Rydberg site.
  /// @see entanglementToNearestStorageSite
  SLMMap<std::vector<std::vector<SLMMap<std::vector<std::vector<Site>>>>>>
      storageToNearestEntanglementSite;

  /// Creates an Architecture from a file containing a JSON specification.
  /// @param filename the name of the file given as a string
  static auto fromJSONFile(const std::string& filename) -> Architecture {
    return fromJSONFile(std::filesystem::path(filename));
  }
  /// Creates an Architecture from a file containing a JSON specification.
  /// @param filepath the path to the file
  static auto fromJSONFile(const std::filesystem::path& filepath)
      -> Architecture {
    return fromJSON(std::ifstream(filepath));
  }
  /// Creates an Architecture from a JSON specification.
  /// @param ifs the input stream containing the JSON specification
  static auto fromJSON(std::istream&& ifs) -> Architecture {
    return fromJSON(nlohmann::json::parse(std::move(ifs)));
  }
  /// Creates an Architecture from a JSON specification.
  /// @param json the JSON string
  static auto fromJSON(std::string&& json) -> Architecture {
    return fromJSON(nlohmann::json::parse(json));
  }
  /// Creates an Architecture from a JSON specification.
  /// @param json the JSON specification
  static auto fromJSON(nlohmann::json json) -> Architecture;
  Architecture() = default;
  // Explicitly delete copy constructor and copy assignment operator because the
  // SLMs and AODs owned by the Architecture cannot be copied.
  Architecture(const Architecture&) = delete;
  Architecture& operator=(const Architecture&) = delete;
  // Default move constructor and move assignment operator
  Architecture(Architecture&&) noexcept = default;
  Architecture& operator=(Architecture&&) noexcept = default;
  /// Export the architecture for the NAViz tool.
  /// @returns a string containing the NAViz-compatible architecture
  /// specification
  auto exportNAVizMachine() const -> std::string;
  /// Export the architecture for the NAViz tool.
  /// @param os the output stream to write the NAViz-compatible architecture
  auto exportNAVizMachine(std::ostream&& os) const -> void {
    os << exportNAVizMachine();
  }
  /// Export the architecture for the NAViz tool.
  /// @param os the output stream to write the NAViz-compatible architecture
  auto exportNAVizMachine(std::ostream& os) const -> void {
    exportNAVizMachine(std::move(os));
  }
  /// Export the architecture for the NAViz tool.
  /// @param path the path to the `.namachine`-file to write the
  /// NAViz-compatible architecture
  auto exportNAVizMachine(const std::filesystem::path& path) const -> void {
    exportNAVizMachine(std::ofstream(path));
  }
  /// Export the architecture for the NAViz tool.
  /// @param filename the name of the `.namachine`-file to write the
  /// NAViz-compatible architecture
  auto exportNAVizMachine(const std::string& filename) const -> void {
    exportNAVizMachine(std::filesystem::path(filename));
  }
  /// Check if the given position is a valid SLM position, i.e., whether the
  /// given row and column are within the range of the SLM.
  [[nodiscard]] auto isValidSlmPosition(const SLM& slm, std::size_t r,
                                        std::size_t c) const -> bool;
  /// Compute the exact location of the SLM site given the row and column
  /// indices expressed in coordinates in the global coordinate system.
  [[nodiscard]] auto exactSlmLocation(const SLM& slm, std::size_t r,
                                      std::size_t c) const
      -> std::pair<std::size_t, std::size_t>;
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
  /// Returns the other site of a pair of entanglement sites
  auto otherEntanglementSite(const SLM& slm, std::size_t r, std::size_t c) const
      -> std::tuple<std::reference_wrapper<const SLM>, std::size_t,
                    std::size_t>;

private:
  /// Compute the site region for entanglement zone and the nearest Rydberg site
  /// for each storage site.
  /// @note We assume we only have one storage zone or one entanglement zone per
  /// row.
  auto preprocessing() -> void;
  /**
   * In the loop, we will calculate a lower bound of the distance
   * between the entanglement site and a storage SLM. Any site in the
   * storage SLM will have at least this distance to the entanglement
   * site. This distance will be the variable @c minimalDistance.
   * Among all storage SLMs, we will find the one that has the minimum
   * distance to the entanglement site, the @c minimumDistance.
   * @note those functions are meant to be used in @ref preprocessing.
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
   * @note those functions are meant to be used in @ref preprocessing.
   */
  [[nodiscard]] auto findNearestEntanglementSLM(size_t x, size_t y,
                                                size_t otherX,
                                                size_t otherY) const
      -> const SLM&;
};
} // namespace na::zoned
