#include "na/azac/Utils.hpp"

#include <algorithm>
#include <cstddef>
#include <deque>
#include <numeric>
#include <optional>
#include <queue>
#include <stack>
#include <tuple>
#include <stdexcept>
#include <utility>
#include <vector>

namespace na {
auto maximumBipartiteMatching(
    const std::vector<std::vector<std::size_t>>& sparseMatrix,
    const bool inverted) -> std::vector<std::optional<std::size_t>> {
  const auto maxSink = std::accumulate(
      sparseMatrix.cbegin(), sparseMatrix.cend(), static_cast<std::size_t>(0),
      [](const std::size_t max, const std::vector<std::size_t>& row) {
        return std::max(max, *std::max_element(row.cbegin(), row.cend()));
      });
  std::vector freeSources(sparseMatrix.size(), true);
  std::vector<std::optional<std::size_t>> invMatching(maxSink + 1, std::nullopt);
  while (true) {
    // find the reachable free sinks on shortest augmenting paths via bfs
    // for all distances, std::nullopt means "not visited yet", i.e., infinite
    // distance
    std::vector<std::optional<std::size_t>> distance(sparseMatrix.size(),
                                                     std::nullopt);
    for (std::size_t s = 0; s < freeSources.size(); ++s) {
      if (freeSources[s]) {
        distance[s] = 0;
      }
    }
    std::queue<std::size_t> queue{};
    for (std::size_t source = 0; source < freeSources.size(); ++source) {
      if (freeSources[source]) {
        queue.push(source);
      }
    }
    std::optional<std::size_t> maxDistance = std::nullopt;
    while (!queue.empty()) {
      const auto source = queue.front();
      queue.pop();
      if (!maxDistance || *distance[source] <= *maxDistance) {
        for (const auto sink : sparseMatrix[source]) {
          if (invMatching[sink]) { // a matched sink is found
            const auto nextSource = *invMatching[sink];
            if (!distance[nextSource]) { // nextSource is not visited yet
              distance[nextSource] = *distance[source] + 1;
              queue.push(nextSource);
            }
          } else { // a free sink is found
            maxDistance = distance[source];
          }
        }
      }
    }
    if (!maxDistance) { // no augmenting path exists
      break;
    }
    // find the augmenting paths via dfs and update the matching
    for (std::size_t freeSource = 0; freeSource < freeSources.size(); ++freeSource) {
      if (freeSources[freeSource]) {
        std::stack stack(std::deque{freeSource});
        // this vector tracks the predecessors of each source, i.e., the sink
        // AND the source coming before the source in the augmenting path
        std::vector<std::optional<std::pair<std::size_t, std::size_t>>> parents(
            sparseMatrix.size(), std::nullopt);
        std::optional<std::pair<std::size_t, std::size_t>> freeSinkFound =
            std::nullopt;
        while (!freeSinkFound && !stack.empty()) {
          const auto source = stack.top();
          stack.pop();
          for (const auto sink : sparseMatrix[source]) {
            if (invMatching[sink]) { // a matched sink is found
              const auto nextSource = *invMatching[sink];
              if (distance[nextSource] &&
                  *distance[nextSource] == *distance[source] + 1) {
                // the edge from source to sink is a valid edge that was
                // encountered during the bfs
                parents[nextSource] = {source, sink};
                stack.push(nextSource);
                  }
            } else { // a free sink is found
              freeSinkFound = {source, sink};
              break;
            }
          }
          distance[source] = std::nullopt; // mark source as visited
        }
        if (freeSinkFound) {
          // augment the matching
          auto source = freeSinkFound->first;
          auto sink = freeSinkFound->second;
          invMatching[sink] =
              source; // that is the additional edge in the matching
          while (source != freeSource) {
            sink = parents[source]->second;
            source = parents[source]->first;
            invMatching[sink] =
                source; // update the matching, i.e., flip the edge from the
            // successor to the predecessor
          }
          freeSources[freeSource] = false;
        }
      }
    }
  }
  // ===-------------------------------===
  if (inverted) {
    return invMatching;
  }
  // invert the matching
  std::vector<std::optional<std::size_t>> matching(sparseMatrix.size(),
                                                   std::nullopt);
  for (std::size_t i = 0; i < invMatching.size(); ++i) {
    if (invMatching[i]) {
      matching[*invMatching[i]] = i;
    }
  }
  return matching;
}

auto minimumWeightFullBipartiteMatching(
    const std::vector<std::vector<std::optional<double>>>& cost_matrix)
    -> std::vector<std::size_t> {
  const std::size_t sizeX = cost_matrix.size();
  auto it = cost_matrix.cbegin();
  const std::size_t sizeY = it->size();
  // check the rectangular shape of input matrix, i.e., check whether all
  // consecutive rows have the same size
  for (++it; it != cost_matrix.cend(); ++it) {
    if (it->size() != sizeY) {
      throw std::invalid_argument("Input matrix must be rectangular");
    }
  }
  // for all x lists all neighbors y in increasing order of c(x, y)
  std::vector list(sizeX, std::vector<std::size_t>{});
  for (std::size_t x = 0; x < sizeX; ++x) {
    for (std::size_t y = 0; y < sizeY; ++y) {
      list[x].emplace_back(y);
    }
    std::sort(list[x].begin(), list[x].end(),
              [x, &cost_matrix](const std::size_t a, const std::size_t b) {
                return cost_matrix[x][a] < cost_matrix[x][b];
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
        std::tuple<std::optional<std::vector<std::size_t>::const_iterator>,
                   std::size_t, std::size_t, double>>
        queue{};
    std::vector residueSetX = freeSources;
    std::vector residueSetY(sizeY, false);
    for (std::size_t x = 0; x < sizeX; ++x) {
      if (residueSetX[x]) {
        quantitiesX[x] = 0.0;
        const auto listIt = list[x].cbegin();
        const auto y = *listIt;
        queue.emplace(listIt, x, y,
                      quantitiesX[x] + *cost_matrix[x][y] + potentialsX[x] -
                          maxPotential);
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
        auto optIt = std::get<0>(itm);
        special = optIt.has_value();
        x = std::get<1>(itm);
        y = std::get<2>(itm);
        queue.pop();
        if (special) {
          if (list[x].back() != y) {
            const auto itW = ++(*optIt);
            const auto w = *itW;
            queue.emplace(itW, x, w,
                          quantitiesX[x] + *cost_matrix[x][w] +
                              potentialsX[x] - maxPotential);
          }
          queue.emplace(std::nullopt, x, y,
                        quantitiesX[x] + *cost_matrix[x][y] +
                            potentialsX[x] - potentialsY[y]);
        }
      } while (!special && invMatching[y] != x);
      // select regular item from queue - done
      if (!residueSetY[y]) {
        pathSetY[y] = x;
        residueSetY[y] = true;
        intersection[y] = freeDestinations[y];
        quantitiesY[y] = quantitiesX[x] + *cost_matrix[x][y] +
                         potentialsX[x] - potentialsY[y];
        if (!freeDestinations[y]) {
          const auto v = *invMatching[y];
          pathSetX[v] = y;
          residueSetX[v] = true;
          quantitiesX[v] = quantitiesY[y];
          const auto itW = list[v].cbegin();
          const auto w = *itW;
          queue.emplace(itW, v, w,
                        quantitiesX[v] + *cost_matrix[v][w] +
                            potentialsX[v] - maxPotential);
        }
      }
    }
    for (std::size_t v = 0; v < sizeX + sizeY; ++v) {
      if (residueSetX[v]) {
        potentialsX[v] = quantitiesY[y];
      } else {
        potentialsX[v] += potentialsY[y];
      }
      maxPotential = std::max(potentialsX[v], maxPotential);
    }
    while (true) {
      x = pathSetY[y];
      const bool freeSourceFound = freeSources[x];
      invMatching[y] = x;
      ++sizeMatching;
      freeSources[x] = false;
      freeDestinations[y] = false;
      intersection[y] = false;
      if (freeSourceFound) {
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
