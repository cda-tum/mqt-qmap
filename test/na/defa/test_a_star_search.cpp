#include "na/defa/AStarSearch.hpp"

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
  std::vector<std::unique_ptr<Node>> nodes{};
  nodes.reserve(31);
  for (std::size_t j = 0; j < 16; ++j) {
    nodes.emplace_back(
        std::make_unique<Node>(2 * j, 4, std::vector<const Node*>{}));
  }
  for (std::size_t i = 4; i > 0; --i) {
    for (std::size_t j = 0; j < 1 << (i - 1); ++j) {
      nodes.emplace_back(std::make_unique<Node>(
          (1 << (5 - i)) - 1 + (j * (1 << (6 - i))), i - 1,
          std::vector<const Node*>{nodes[32 - (1 << (i + 1)) + (2 * j)].get(),
                                   nodes[32 - (1 << (i + 1)) + (2 * j) + 1].get()}));
    }
  }
  ASSERT_NO_THROW({
    const auto path =
        (na::aStarTreeSearch<Node, std::vector<const Node*>::const_iterator>(
            /* start: */
            nodes[30].get(),
            /* getNeighbors: */
            [](const Node* node)
                -> std::pair<std::vector<const Node*>::const_iterator,
                             std::vector<const Node*>::const_iterator> {
              return {node->neighbors.cbegin(), node->neighbors.cend()};
            },
            /* isGoal: */
            [](const Node* node) -> bool {
              return node->x == 8 && node->y == 4;
            },
            /* getCost: */
            [](const Node* /* unused */, const Node* /* unused */) -> double {
              return 1.0;
            },
            /* getHeuristic: */
            [](const Node* node) -> double {
              return std::sqrt(std::pow(8.0 - static_cast<double>(node->x), 2) +
                               std::pow(4.0 - static_cast<double>(node->y), 2));
            }));
    EXPECT_EQ(path.size(), 5);
    EXPECT_THAT(path, testing::AnyOf(std::vector<const Node*>{
                          nodes[30].get(), nodes[28].get(), nodes[25].get(),
                          nodes[18].get(), nodes[4].get()}));
  });
}