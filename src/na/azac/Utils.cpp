#include "na/azac/Utils.hpp"

#include <algorithm>
#include <cstddef>
#include <deque>
#include <functional>
#include <numeric>
#include <optional>
#include <queue>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace na {
auto minimumWeightFullBipartiteMatching(
    const std::vector<std::vector<std::optional<double>>>& costMatrix)
    -> std::vector<std::size_t> {
  const std::size_t sizeX = costMatrix.size();
  if (sizeX == 0) {
    return {};
  }
  auto it = costMatrix.cbegin();
  const std::size_t sizeY = it->size();
  if (sizeX > sizeY) {
    throw std::invalid_argument(
        "Input matrix must have more columns than rows");
  }
  if (std::all_of(it->cbegin(), it->cend(), [](const std::optional<double>& o) {
        return o == std::nullopt;
      })) {
    throw std::invalid_argument("Input matrix must not contain empty rows");
  }
  // check the rectangular shape of input matrix, i.e., check whether all
  // consecutive rows have the same size
  for (++it; it != costMatrix.cend(); ++it) {
    if (it->size() != sizeY) {
      throw std::invalid_argument("Input matrix must be rectangular");
    }
    if (std::all_of(
            it->cbegin(), it->cend(),
            [](const std::optional<double>& o) { return o == std::nullopt; })) {
      throw std::invalid_argument("Input matrix must not contain empty rows");
    }
  }
  // for all x lists all neighbors y in increasing order of c(x, y)
  std::vector list(sizeX, std::vector<std::size_t>{});
  for (std::size_t x = 0; x < sizeX; ++x) {
    for (std::size_t y = 0; y < sizeY; ++y) {
      if (costMatrix[x][y]) {
        list[x].emplace_back(y);
      }
    }
    std::sort(list[x].begin(), list[x].end(),
              [x, &costMatrix](const std::size_t a, const std::size_t b) {
                return costMatrix[x][a] < costMatrix[x][b];
              });
  }
  // initialize the set of free sources
  std::vector freeSources(sizeX, true);
  // initialize the set of free targets
  std::vector freeDestinations(sizeY, true);
  // initialize the matching
  std::vector<std::optional<std::size_t>> invMatching(sizeY, std::nullopt);
  std::size_t sizeMatching = 0;
  std::vector quantitiesX(sizeX, 0.0);
  std::vector quantitiesY(sizeY, 0.0);
  std::vector potentialsX(sizeX, 0.0);
  std::vector potentialsY(sizeY, 0.0);
  double maxPotential = 0.0;
  while (sizeMatching < sizeX) {
    std::vector<std::size_t> pathSetX(sizeX, 0);
    std::vector<std::size_t> pathSetY(sizeY, 0);
    // items have the form (<special>, <x>, <y>, <cost>)
    // if special, the iterator to the edge in list is stored in the
    // optional
    std::priority_queue<
        std::tuple<double, std::size_t, std::size_t,
                   std::optional<std::vector<std::size_t>::const_iterator>>,
        std::vector<std::tuple<
            double, std::size_t, std::size_t,
            std::optional<std::vector<std::size_t>::const_iterator>>>,
        std::greater<std::tuple<
            double, std::size_t, std::size_t,
            std::optional<std::vector<std::size_t>::const_iterator>>>>
        queue{};
    std::vector residueSetX = freeSources;
    std::vector residueSetY(sizeY, false);
    for (std::size_t x = 0; x < sizeX; ++x) {
      if (residueSetX[x]) {
        quantitiesX[x] = 0.0;
        const auto listIt = list[x].cbegin();
        const auto y = *listIt;
        queue.emplace(quantitiesX[x] + *costMatrix[x][y] + potentialsX[x] -
                          maxPotential,
                      x, y, listIt);
      }
    }
    // intersection of `remainder_set` and `freeDestinations` is empty
    std::vector intersection(sizeY, false);
    std::size_t x = 0;
    std::size_t y = 0;
    while (std::all_of(intersection.cbegin(), intersection.cend(),
                       [](const bool b) { return !b; })) {
      // select regular item from queue
      bool special = false;
      do {
        const auto& itm = queue.top();
        auto optIt = std::get<3>(itm);
        special = optIt.has_value();
        x = std::get<1>(itm);
        y = std::get<2>(itm);
        queue.pop();
        if (special) {
          if (list[x].back() != y) {
            const auto w = *(++(*optIt));
            queue.emplace(quantitiesX[x] + *costMatrix[x][w] + potentialsX[x] -
                              maxPotential,
                          x, w, *optIt);
          }
          queue.emplace(quantitiesX[x] + *costMatrix[x][y] + potentialsX[x] -
                            potentialsY[y],
                        x, y, std::nullopt);
        }
      } while (special || invMatching[y] == x);
      // select regular item from queue - done
      if (!residueSetY[y]) {
        pathSetY[y] = x;
        residueSetY[y] = true;
        intersection[y] = freeDestinations[y];
        quantitiesY[y] = quantitiesX[x] + *costMatrix[x][y] + potentialsX[x] -
                         potentialsY[y];
        if (!freeDestinations[y]) {
          const auto v = *invMatching[y];
          pathSetX[v] = y;
          residueSetX[v] = true;
          quantitiesX[v] = quantitiesY[y];
          const auto itW = list[v].cbegin();
          const auto w = *itW;
          queue.emplace(quantitiesX[v] + *costMatrix[v][w] + potentialsX[v] -
                            maxPotential,
                        v, w, itW);
        }
      }
    }
    // reset maxPotential and update quantities and potentials
    // afterward maxPotential will have the correct value again
    maxPotential = std::numeric_limits<double>::min();
    for (std::size_t v = 0; v < sizeY; ++v) {
      if (!residueSetY[v]) {
        quantitiesY[v] = quantitiesY[y];
      }
      potentialsY[v] += quantitiesY[v];
      maxPotential = std::max(potentialsY[v], maxPotential);
    }
    for (std::size_t v = 0; v < sizeX; ++v) {
      if (!residueSetX[v]) {
        quantitiesX[v] = quantitiesY[y];
      }
      potentialsX[v] += quantitiesX[v];
      maxPotential = std::max(potentialsX[v], maxPotential);
    }
    freeDestinations[y] = false;
    ++sizeMatching;
    while (true) {
      x = pathSetY[y];
      invMatching[y] = x;
      if (freeSources[x]) {
        freeSources[x] = false;
        break;
      }
      y = pathSetX[x];
    }
  }
  std::vector<std::size_t> matching(sizeX, 0);
  for (std::size_t y = 0; y < sizeY; ++y) {
    if (const auto optX = invMatching[y]) {
      matching[*optX] = y;
    }
  }
  return matching;
}
} // namespace na
