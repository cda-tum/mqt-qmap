#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace na {

template <typename T1, typename T2>
auto distance(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b)
    -> double {
  return std::sqrt(
      std::pow(static_cast<double>(a.first) - static_cast<double>(b.first), 2) +
      std::pow(static_cast<double>(a.second) - static_cast<double>(b.second),
               2));
}

/// Computes a maximum matching in a bipartite graph
/// @note implemented pseudocode from
/// https://epubs.siam.org/doi/pdf/10.1137/0202019?download=true
auto maximumBipartiteMatching(
    const std::vector<std::vector<std::size_t>>& sparseMatrix,
    bool inverted = false) -> std::vector<std::optional<std::size_t>>;

/// @note implemented following pseudocode in
/// https://www2.eecs.berkeley.edu/Pubs/TechRpts/1978/ERL-m-78-67.pdf
auto minimumWeightFullBipartiteMatching(
    const std::vector<std::vector<std::optional<double>>>& costMatrix)
    -> std::vector<std::size_t>;

} // namespace na