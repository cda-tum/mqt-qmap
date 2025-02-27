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
 * @brief A* search algorithm for trees where neighbors are sorted by cost
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
 * - To keep the maintenance of the open set as simple as possible, the open set
 * only stores one neighbor of its parent node at a time. To achieve this, the
 * @p getNeighbors function must return the neighbors in increasing order of
 * cost. The first neighbor is then the one with the lowest cost and is
 * placed in the priority queue. If this neighbor is not the only one, it is
 * placed in the queue as a special item meaning that it represents also all
 * other neighbors with higher cost. When this item is popped from the queue,
 * it is replaced by the next neighbor with higher cost.
 * @note This implementation of A* search can only handle trees and not general
 * graphs. This is because it does not keep track of visited nodes and therefore
 * cannot detect cycles. Also for DAGs it may expand nodes multiple times when
 * they can be reached by different paths from the start node.
 * @note The function @p getNeighbors must return the neighbors of a node in a
 * sorted order by cost. The neighbor with the lowest cost must be the first
 * element in the range of neighbors.
 * @note @p getHeuristic must be admissible, meaning that it never
 * overestimates the cost to reach the goal from the current node calculated by
 * @p getCost for every edge on the path.
 * @note The calling program has to make sure that the pointers and iterators
 * passed to this function are valid and that the iterators are not invalidated
 * during the search, e.g., by calling one of the passed functions like @p
 * getNeighbors.
 * @param start a pointer to the start node
 * @param getNeighbors a function that returns the neighbors of a node as an
 * iterator pair, the first iterator must be the 'begin' iterator and the second
 * iterator must be the 'end' iterator. The neighbors must be sorted by cost,
 * see the note above.
 * @param isGoal a function that returns true if a node is one of potentially
 * multiple goals
 * @param getCost a function that returns the cost of moving from one node to
 * another
 * @param getHeuristic a function that returns the heuristic cost from the node
 * to any goal.
 * @return a vector of node pointers representing the path from the start to a
 * goal
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
    double priority;  //< sum of cost and heuristic
    double cost;      //< actual cost to reach the node
    const Node* node; //< pointer to the node
    // iterator pair to the siblings of a node with a 'begin' and 'end' iterator
    // if the iterator is 'empty' the optional is set to 'nullopt' and then
    // this item represents a regular item. If it is not empty, the item
    // is a special item and represents the contained node together with all
    // its more costly siblings.
    std::optional<std::pair<ForwardIt, ForwardIt>> siblings;
    // pointer to the parent item to reconstruct the path in the end
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
                                                           0.0, start, nullptr))
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
      // if the sibling is the last one, create a regular item otherwise create
      // a special item together with the advanced iterator
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
      // getCost returns the cost for the edge from the current node (itm->node)
      // to the neighbor (firstNeighbor).
      // hence, the total cost is the cost to reach the current node (itm->cost)
      // plus the cost to reach the neighbor (getCost(itm->node, firstNeighbor))
      const auto cost = itm->cost + getCost(itm->node, firstNeighbor);
      const auto heuristic = getHeuristic(firstNeighbor);
      // if the neighbor is the last one, create a regular item otherwise create
      // a special item together with the advanced iterator
      if (it == end) {
        openSet.emplace(items
                            .emplace_back(std::make_unique<Item>(
                                cost + heuristic, cost, firstNeighbor, itm))
                            .get());
      } else {
        openSet.emplace(
            items
                .emplace_back(std::make_unique<Item>(
                    cost + heuristic, cost, firstNeighbor, it, end, itm))
                .get());
      }
    }
  }
  throw std::runtime_error("No path from start to any goal found.");
}
} // namespace na