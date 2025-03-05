#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <queue>
#include <stdexcept>
#include <unordered_map>
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

template <class Priority, class T, class Compare = std::less<Priority>>
/// A class that implements a heap data structure.
/// @details The heap is a container that provides constant time lookup of the
/// largest (by default) element, at the expense of logarithmic insertion and
/// extraction. A user-provided `Compare` can be supplied to change the
/// ordering, e.g., using `std::greater<T>` would cause the smallest element to
/// appear as the top(). Opposed to `std::priority_queue`, this heap allows for
/// updating the priority of an element in O(log(n)) time. Additionally, it
/// allows for erasing an element in O(log(n)) time and elements are unique.
/// This means, if an element is pushed that is already in the heap, the
/// priority of the existing element is updated.
/// @tparam Priority The type of the priority of the elements.
/// @tparam T The type of the elements.
/// @tparam Compare The type of the comparison function.
class Heap {
public:
  using PriorityType = Priority;
  using ElementType = T;
  using ValueType = std::pair<PriorityType, ElementType>;
  using SizeType = typename std::vector<ValueType>::size_type;
  using PriorityCompare = Compare;
  using Reference = ValueType&;

private:
  std::vector<ValueType> heap_;
  std::unordered_map<ElementType, SizeType> keyToIndex_;

  /// Moves the element at index i up the heap until the heap property is
  /// satisfied.
  /// @param i The index of the element to move up.
  void heapifyUp(size_t i) {
    while (i > 0) {
      size_t parent = (i - 1) / 2;
      if (PriorityCompare{}(heap_[i].first, heap_[parent].first)) {
        std::swap(heap_[i], heap_[parent]);
        keyToIndex_[heap_[i]] = i;
        keyToIndex_[heap_[parent]] = parent;
        i = parent;
      } else {
        break;
      }
    }
  }

  /// Moves the element at index i down the heap until the heap property is
  /// satisfied.
  /// @param i The index of the element to move down.
  void heapifyDown(size_t i) {
    while (true) {
      size_t leftChild = (2 * i) + 1;
      size_t rightChild = (2 * i) + 2;
      size_t smallest = i;

      if (leftChild < heap_.size() &&
          PriorityCompare{}(heap_[leftChild], heap_[smallest])) {
        smallest = leftChild;
      }
      if (rightChild < heap_.size() &&
          PriorityCompare{}(heap_[rightChild], heap_[smallest])) {
        smallest = rightChild;
      }
      if (smallest != i) {
        std::swap(heap_[i], heap_[smallest]);
        keyToIndex_[heap_[i]] = i;
        keyToIndex_[heap_[smallest]] = smallest;
        i = smallest;
      } else {
        break;
      }
    }
  }

public:
  /// Returns a reference to the top element of the heap.
  /// @details Takes O(1) time.
  /// @return A reference to the top element of the heap.
  [[nodiscard]] auto top() const -> const Reference { return heap_.front(); }
  /// Removes the top element from the heap.
  /// @details Takes O(log(n)) time.
  auto pop() -> void {
    keyToIndex_.erase(heap_.front().second);
    std::swap(heap_.front(), heap_.back());
    heap_.pop_back();
    keyToIndex_.at(heap_.front().second) = 0;
    heapifyDown(0);
  }
  /// Checks if the heap is empty.
  /// @return True if the heap is empty, false otherwise.
  [[nodiscard]] auto empty() const -> bool { return heap_.empty(); }
  /// Returns the number of elements in the heap.
  /// @return The number of elements in the heap.
  [[nodiscard]] auto size() const -> SizeType { return heap_.size(); }
  /// Adds an element to the heap.
  /// @details Takes O(log(n)) time.
  /// @param value The priority-element pair to add to the heap.
  auto push(const ValueType& value) -> void {
    if (keyToIndex_.find(value.second) != keyToIndex_.end()) {
      update(value.second, value.first);
    }
    heap_.push_back(value);
    keyToIndex_.at(value.second) = heap_.size() - 1;
    heapifyUp(heap_.size() - 1);
  }
  /// Constructs a new element in place in the heap.
  /// @details Takes O(log(n)) time.
  /// @param args The arguments to construct the priority-element pair.
  /// @return A reference to the newly constructed priority-element pair.
  template <class... Args> auto emplace(Args&&... args) -> const Reference {
    auto value = value_type(std::forward<Args>(args)...);
    if (keyToIndex_.find(value.second) != keyToIndex_.end()) {
      return update(std::move(value));
    }
    push(std::move(value));
    return heap_[keyToIndex_.at(value.second)];
  }
  /// Updates the priority of an element in the heap.
  /// @details Takes O(log(n)) time.
  /// @param value The priority-element pair with the new priority for the
  /// element to be updated.
  auto update(const ValueType& value) -> const Reference {
    const auto i = keyToIndex_.find(value.second)->second;
    heap_[i] = value;
    keyToIndex_.at(value.second) = i;
    // for the case that the priority is increased
    heapifyUp(i);
    // for the case that the priority is decreased
    heapifyDown(i);
    return heap_[keyToIndex_.at(value.second)];
  }
  /// Eraes an element from the heap.
  /// @details Takes O(log(n)) time.
  /// @param element The element to erase.
  auto erase(const ElementType& element) -> void {
    const auto i = keyToIndex_.find(element)->second;
    keyToIndex_.erase(element);
    std::swap(heap_[i], heap_.back());
    heap_.pop_back();
    keyToIndex_.at(heap_[i].second) = i;
    heapifyDown(i);
  }
};

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
