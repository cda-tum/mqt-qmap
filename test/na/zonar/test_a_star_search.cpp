#include "na/zonar/AStarSearch.hpp"

#include <cmath>
#include <cstddef>
#include <gmock/gmock-matchers.h>
#include <gtest/gtest.h>
#include <memory>
#include <utility>
#include <vector>

TEST(AStarSearch, Grid) {
  struct Node {
    std::size_t x;
    std::size_t y;
    std::vector<const Node*> neighbors;
    Node(const std::size_t x, const std::size_t y,
         std::vector<const Node*> neighbors)
        : x(x), y(y), neighbors(std::move(neighbors)) {}
  };
  std::vector<std::vector<std::unique_ptr<Node>>> nodes(4);
  nodes[0].emplace_back(
      std::make_unique<Node>(3, 3, std::vector<const Node*>{}));
  nodes[0].reserve(4);
  for (std::size_t i = 1; i < 4; ++i) {
    nodes[0].emplace_back(std::make_unique<Node>(
        3, 3 - i, std::vector<const Node*>{nodes[0].back().get()}));
    nodes[i].reserve(4);
    nodes[i].emplace_back(std::make_unique<Node>(
        3 - i, 3, std::vector<const Node*>{nodes[i - 1][0].get()}));
  }
  for (std::size_t i = 1; i < 4; ++i) {
    for (std::size_t j = 1; j < 4; ++j) {
      nodes[i].emplace_back(std::make_unique<Node>(
          3 - i, 3 - j,
          std::vector<const Node*>{nodes[i].back().get(),
                                   nodes[i - 1][j].get()}));
    }
  }
  const auto path =
      na::aStarTreeSearch<Node, std::vector<const Node*>::const_iterator>(
          /* start: */
          nodes[3][3].get(),
          /* getNeighbors: */
          [](const Node* node)
              -> std::pair<std::vector<const Node*>::const_iterator,
                           std::vector<const Node*>::const_iterator> {
            return {node->neighbors.cbegin(), node->neighbors.cend()};
          },
          /* isGoal: */
          [](const Node* node) -> bool { return node->x == 3 && node->y == 1; },
          /* getCost: */
          [](const Node* /* unused */, const Node* /* unused */) -> double {
            return 1.0;
          },
          /* getHeuristic: */
          [](const Node* node) -> double {
            return std::sqrt(std::pow(3.0 - static_cast<double>(node->x), 2) +
                             std::pow(1.0 - static_cast<double>(node->y), 2));
          });
  EXPECT_EQ(path.size(), 5);
  EXPECT_THAT(path,
              testing::AnyOf(
                  std::vector<const Node*>{nodes[3][3].get(), nodes[2][3].get(),
                                           nodes[1][3].get(), nodes[0][3].get(),
                                           nodes[0][2].get()},
                  std::vector<const Node*>{nodes[3][3].get(), nodes[2][3].get(),
                                           nodes[1][3].get(), nodes[1][2].get(),
                                           nodes[0][2].get()},
                  std::vector<const Node*>{nodes[3][3].get(), nodes[2][3].get(),
                                           nodes[2][2].get(), nodes[1][2].get(),
                                           nodes[0][2].get()},
                  std::vector<const Node*>{nodes[3][3].get(), nodes[3][2].get(),
                                           nodes[2][2].get(), nodes[1][2].get(),
                                           nodes[0][2].get()}));
}