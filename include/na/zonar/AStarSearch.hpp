#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace na {
template <class Node, class ForwardIt>

/**
 * @brief A* search algorithm
 * @details A* is a graph traversal and path search algorithm that finds the
 * shortest path between a start node and a goal node. It evaluates nodes by
 * combining the cost to reach the node and the cost to get from the node to the
 * goal estimated by a heuristic function.
 * @param start the start node
 * @param getNeighbors a function that returns the neighbors of a node
 * @param isGoal a function that returns true if a node is the goal
 * @param getCost a function that returns the cost of moving from one node to
 * another
 * @param getHeuristic a function that returns the heuristic cost from the node
 * to any goal.
 * @return a vector of nodes representing the path from the start to a goal
 * @note @p getHeuristic must be admissible, meaning that it never
 * overestimates the cost to reach the goal from the current node calculated by
 * @p getCost for every edge on the path.
 * @note Note that there must be a hash function for the Node class defined.
 * It is also useful to define the == operator for the Node class.
 */
auto aStarTreeSearch(
    const Node* start,
    std::function<std::pair<ForwardIt, ForwardIt>(const Node*)> getNeighbors,
    std::function<bool(const Node*)> isGoal,
    std::function<double(const Node*, const Node*)> getCost,
    std::function<double(const Node*)> getHeuristic)
    -> std::vector<const Node*> {
  //===--------------------------------------------------------------------===//
  // Setup open set structure
  //===--------------------------------------------------------------------===//
  // struct for items in the open set
  struct Item {
    double priority; //< sum of cost and heuristic
    double cost;
    const Node* node;
    std::optional<std::pair<ForwardIt, ForwardIt>> siblings;
    Item* parent;
    Item(const double priority, const double cost, const Node* node,
         ForwardIt siblings, ForwardIt end, Item* parent)
        : priority(priority), cost(cost), node(node),
          siblings(std::pair{siblings, end}), parent(parent) {}
    Item(const double priority, const double cost, const Node* node,
         Item* parent)
        : priority(priority), cost(cost), node(node), siblings(std::nullopt),
          parent(parent) {}
  };
  // compare function for the open set
  struct ItemCompare {
    bool operator()(const Item* a, const Item* b) const {
      return a->priority > b->priority;
    }
  };
  std::vector<std::unique_ptr<Item>> items;
  // open list of nodes to be evaluated as a minimum heap based on the priority
  std::priority_queue<Item*, std::vector<Item*>, ItemCompare> openSet;
  openSet.emplace(
      items.emplace_back(std::make_unique<Item>(getHeuristic(start), 0.0, start, nullptr))
          .get());
  //===--------------------------------------------------------------------===//
  // Perform A* search
  //===--------------------------------------------------------------------===//
  while (!openSet.empty()) {
    Item* itm = openSet.top();
    openSet.pop();
    // if a goal is reached, that is the shortest path to a goal under the
    // assumption that the heuristic is admissible
    if (isGoal(itm->node)) {
      // reconstruct the path from the goal to the start and then reverse it
      std::vector<const Node*> path;
      for (; itm != nullptr; itm = itm->parent) {
        path.emplace_back(itm->node);
      }
      std::reverse(path.begin(), path.end());
      return path;
    }
    // replace the entry in the open set representing the popped item including
    // all its siblings with higher cost
    if (itm->siblings) {
      const Node* nextSibling = *itm->siblings->first++;
      const auto cost = itm->parent->cost + getCost(itm->node, nextSibling);
      const auto heuristic = getHeuristic(nextSibling);
      if (itm->siblings->first == itm->siblings->second) {
        openSet.emplace(
            items
                .emplace_back(std::make_unique<Item>(cost + heuristic, cost,
                                                     nextSibling, itm->parent))
                .get());
      } else {
        openSet.emplace(
            items
                .emplace_back(std::make_unique<Item>(
                    cost + heuristic, cost, nextSibling, itm->siblings->first,
                    itm->siblings->second, itm->parent))
                .get());
      }
    }
    // expand the current node by adding all neighbors to the open set in the
    // form of one representative for all neighbors with the cost of the
    // neighbor with the lowest cost
    auto [it, end] = getNeighbors(itm->node);
    if (it != end) {
      const Node* firstNeighbor = *it++;
      const auto cost = itm->cost + getCost(itm->node, firstNeighbor);
      const auto heuristic = getHeuristic(firstNeighbor);
      if (it == end) {
      openSet.emplace(
          items
              .emplace_back(std::make_unique<Item>(cost + heuristic, cost,
                                                   firstNeighbor, itm))
              .get());
      } else {
      openSet.emplace(
          items
              .emplace_back(std::make_unique<Item>(cost + heuristic, cost,
                                                   firstNeighbor, it, end, itm))
              .get());
      }
    }
  }
  throw std::runtime_error("No path from start to any goal found.");
}
} // namespace na