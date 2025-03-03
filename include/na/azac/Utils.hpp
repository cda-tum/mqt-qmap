#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <queue>
#include <stdexcept>
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

/**
 * @brief A* search algorithm for trees
 * @details A* is a graph traversal and path search algorithm that finds the
 * shortest path between a start node and a goal node. It evaluates nodes by
 * combining the cost to reach the node and the cost to get from the node to the
 * goal estimated by a heuristic function.
 * @par
 * This implementation of the A* search algorithm has some particularities:
 * - To increase performance for the special case of a tree, where there cannot
 * be any cycles and a node can only be reached by one path, it does not keep
 * visited nodes. This would require a hash set or similar data structure to
 * store visited nodes and check if a node has already been visited. This
 * check would take at least O(log(n)) time for a hash set and is superfluous
 * for trees.
 * - As a consequence of the first point, this implementation also does not
 * check whether a node is already in the open set. This would also require an
 * O(log(n)) check operation which is not necessary for trees as one path can
 * only reach a node.
 * @note This implementation of A* search can only handle trees and not general
 * graphs. This is because it does not keep track of visited nodes and therefore
 * cannot detect cycles. Also for DAGs it may expand nodes multiple times when
 * they can be reached by different paths from the start node.
 * @note @p getHeuristic must be admissible, meaning that it never
 * overestimates the cost to reach the goal from the current node calculated by
 * @p getCost for every edge on the path.
 * @note The calling program has to make sure that the pointers passed to this
 * function are valid and that the iterators are not invalidated during the
 * search, e.g., by calling one of the passed functions like @p getNeighbors.
 * @param start a pointer to the start node
 * @param getNeighbors a function that returns the neighbors of a node
 * @param isGoal a function that returns true if a node is one of potentially
 * multiple goals
 * @param getCost a function that returns the total cost to reach that
 * particular node from the start node
 * @param getHeuristic a function that returns the heuristic cost from the node
 * to any goal.
 * @return a vector of node pointers representing the path from the start to a
 * goal
 */
template <class Node>
auto aStarTreeSearch(
    const Node* start,
    std::function<std::vector<const Node*>(const Node*)> getNeighbors,
    std::function<bool(const Node*)> isGoal,
    std::function<double(const Node*)> getCost,
    std::function<double(const Node*)> getHeuristic)
    -> std::vector<const Node*> {
  //===--------------------------------------------------------------------===//
  // Setup open set structure
  //===--------------------------------------------------------------------===//
  // struct for items in the open set
  struct Item {
    double priority;  //< sum of cost and heuristic
    const Node* node; //< pointer to the node
    // pointer to the parent item to reconstruct the path in the end
    Item* parent;
    Item(const double priority, const Node* node, Item* parent)
        : priority(priority), node(node), parent(parent) {}
  };
  // compare function for the open set
  struct ItemCompare {
    bool operator()(const Item* a, const Item* b) const {
      // this way, the item with the lowest priority is on top of the heap
      return a->priority > b->priority;
    }
  };
  // vector of items to store all items and keep them alive also after they are
  // popped from the open set.
  // they are required alive to reconstruct the path in the end.
  std::vector<std::unique_ptr<Item>> items;
  // open list of nodes to be evaluated as a minimum heap based on the priority.
  // whenever an item is placed in the queue it is created in the vector `items`
  // before and only a reference is placed in the queue
  std::priority_queue<Item*, std::vector<Item*>, ItemCompare> openSet;
  openSet.emplace(items
                      .emplace_back(std::make_unique<Item>(getHeuristic(start),
                                                           start, nullptr))
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
    // expand the current node by adding all neighbors to the open set
    const auto& neighbors = getNeighbors(itm->node);
    if (!neighbors.empty()) {
      for (const auto* neighbor : neighbors) {
        // getCost returns the total cost to reach the current node
        const auto cost = getCost(neighbor);
        const auto heuristic = getHeuristic(neighbor);
        openSet.emplace(items
                            .emplace_back(std::make_unique<Item>(
                                cost + heuristic, neighbor, itm))
                            .get());
      }
    }
  }
  throw std::runtime_error("No path from start to any goal found.");
}

} // namespace na
