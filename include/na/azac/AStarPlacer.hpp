#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na {
class AStarPlacer {
  using DiscreteSite = std::pair<uint16_t, uint16_t>;

  std::reference_wrapper<const Architecture> architecture_;
  /// If true, during the initial placement the atoms are placed starting in the
  /// last row instead of the first row in the first SLM
  bool reverseInitialPlacement_ = false;
  /// this flag indicates whether the  placement should use a window when
  /// selecting potential free sites
  bool useWindow_ = true;
  /// If the window is used, this denotes the height in
  /// terms of rows of the window centered at the nearest site
  size_t windowHeight_ = 4;
  /// If the window is used, this denotes the width in terms of columns of the
  /// window centered at the nearest site
  size_t windowWidth_ = 4;

  /// A node representing one stage in the process of placing all atoms
  /// that must be moved for the next stage starting from the last mapping
  /// until a new mapping is found satisfying all constraints of the next
  /// stage
  struct Node {
    /// the maximum distance an already placed atom must travel to its
    /// target location
    float maxDistanceOfPlacedAtom = 0.0;
    /// a set of all sites that are already occupied by an atom due to the
    /// current placement
    std::unordered_set<DiscreteSite> consumedFreeSites;
    /// a binary search tree representing the horizontal groups
    /// @see getNeighbors for more details
    std::vector<std::map<uint16_t, uint16_t>> hGroups;
    /// the maximum distance of placed atoms in every group to their
    /// target location
    std::vector<float> maxDistancesOfPlacedAtomsPerHGroup;
    /// @see hGroups
    std::vector<std::map<uint16_t, uint16_t>> vGroups;
    /// @see maxDistancesOfPlacedAtomsPerHGroup
    std::vector<float> maxDistancesOfPlacedAtomsPerVGroup;
  };

  /// A list of all nodes that have been created so far.
  /// This list is dynamically extended when new nodes are created.
  /// This happens when a node is expanded by calling getNeighbors.
  std::vector<std::unique_ptr<Node>> nodes_;

  /// When placing atoms after a rydberg layer bayk in the storage zone, this
  /// struct stores for every such atom all required information, i.e., the
  /// current site and potential target sites ordered by distance (ascending).
  struct AtomJob {
    /// the current site of the atom
    DiscreteSite currentSite;
    /// a struct describing one potential target site
    struct Option {
      /// the target site
      DiscreteSite site;
      /// the distance the atom must travel to reach the target site
      float distance;
    };
    /// a list of all potential target sites ordered by distance (ascending)
    std::vector<Option> options;
  };
  /// a list of all atoms that must be placed in the storage zone after a
  /// rydberg layer
  std::vector<AtomJob> atomJobs_;

  /// When placing gates in the entanglement zone before a rydberg layer, this
  /// struct stores for every such gate all required information, i.e., the
  /// current sites of the corresponding atoms and potential target sites
  /// ordered by distance (ascending).
  struct GateJob {
    /// the current sites of the two atoms
    std::pair<DiscreteSite, DiscreteSite> currentDiscreteSites;
    /// a struct describing one potential target site for each atom
    struct Option {
      /// the target sites for the two atoms
      std::pair<DiscreteSite, DiscreteSite> sites;
      /// the max distance the atoms must travel to reach the target sites
      std::pair<float, float> distance;
    };
    /// a list of all potential target sites ordered by distance (ascending)
    std::vector<Option> options;
  };
  /// a list of all gates that must be placed in the entanglement zone before a
  /// rydberg layer
  std::vector<GateJob> gateJobs_;

  /// The number of either atom or gate jobs that must be performed
  size_t nJobs_ = 0;

protected:
  AStarPlacer(const Architecture& architecture, const nlohmann::json& config);

  /**
   * This function defines the interface of the placer and delegates the
   * placement of the qubits to the other functions.
   * @param nQubits the number of qubits to be placed
   * @param twoQubitGateLayers the qubits that must be placed for each layer
   * @param reuseQubits the qubits that are reused in the next stage
   */
  [[nodiscard]] auto
  place(size_t nQubits,
        const std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>&
            twoQubitGateLayers,
        const std::vector<std::unordered_set<qc::Qubit>>& reuseQubits)
      -> std::vector<std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>;

private:
  /**
   * @brief A* search algorithm for trees
   * @details A* is a graph traversal and path search algorithm that finds the
   * shortest path between a start node and a goal node. It evaluates nodes by
   * combining the cost to reach the node and the cost to get from the node to
   * the goal estimated by a heuristic function.
   * @par
   * This implementation of the A* search algorithm has some particularities:
   * - To increase performance for the special case of a tree, where there
   * cannot be any cycles and a node can only be reached by one path, it does
   * not keep visited nodes. This would require a hash set or similar data
   * structure to store visited nodes and check if a node has already been
   * visited. This check would take at least O(log(n)) time for a hash set and
   * is superfluous for trees.
   * - As a consequence of the first point, this implementation also does not
   * check whether a node is already in the open set. This would also require an
   * O(log(n)) check operation which is not necessary for trees as one path can
   * only reach a node.
   * @note This implementation of A* search can only handle trees and not
   * general graphs. This is because it does not keep track of visited nodes and
   * therefore cannot detect cycles. Also for DAGs it may expand nodes multiple
   * times when they can be reached by different paths from the start node.
   * @note @p getHeuristic must be admissible, meaning that it never
   * overestimates the cost to reach the goal from the current node calculated
   * by
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
   * @param getHeuristic a function that returns the heuristic cost from the
   * node to any goal.
   * @return a vector of node pointers representing the path from the start to a
   * goal
   */
  [[nodiscard]] static auto aStarTreeSearch(
      const Node& start,
      const std::function<std::vector<std::reference_wrapper<const Node>>(
          const Node&)>& getNeighbors,
      const std::function<bool(const Node&)>& isGoal,
      const std::function<double(const Node&)>& getCost,
      const std::function<double(const Node&)>& getHeuristic)
      -> std::vector<std::reference_wrapper<const Node>>;

  /**
   * This function takes a list of atoms together with their current placement
   * and returns two maps from concrete columns and rows to their discrete
   * indices.
   * @param placement is a list of atoms together with their current placement
   * @param atoms is a list of all atoms that must be placed
   * @return a pair of two maps, the first one maps rows to their discrete
   * indices and the second one maps columns to their discrete indices
   */
  [[nodiscard]] auto discretizePlacementOfAtoms(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& placement,
      const std::vector<qc::Qubit>& atoms) const
      -> std::pair<std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint16_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>,
                   std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint16_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>;

  /**
   * This function discretizes the storage zone of the architecture and returns
   * two maps from concrete columns and rows to their discrete indices.
   * @param occupiedSites is a set of occupied sites in the storage zone
   * @return a pair of two maps, the first one maps rows to their discrete
   * indices and the second one maps columns to their discrete indices
   */
  [[nodiscard]] auto discretizeNonOccupiedStorageSites(
      const std::unordered_set<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>,
          std::hash<std::tuple<const SLM&, size_t, size_t>>,
          std::equal_to<std::tuple<const SLM&, size_t, size_t>>>& occupiedSites)
      const
      -> std::pair<std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint16_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>,
                   std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint16_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>;

  /**
   * This function discretizes the entanglement zone of the architecture and
   * returns two maps from concrete columns and rows to their discrete indices.
   * @param occupiedSites is a set of occupied sites in the storage zone
   * @return a pair of two maps, the first one maps rows to their discrete
   * indices and the second one maps columns to their discrete indices
   */
  [[nodiscard]] auto discretizeNonOccupiedEntanglementSites(
      const std::unordered_set<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>,
          std::hash<std::tuple<const SLM&, size_t, size_t>>,
          std::equal_to<std::tuple<const SLM&, size_t, size_t>>>& occupiedSites)
      const
      -> std::pair<std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint16_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>,
                   std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint16_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>;
  /**
   * Given the discretion of the previous placement of the atoms and the
   * discretization of the target slms, this function updates the placement of
   * based on the previous placement for the given list of atoms.
   * @param atoms is a list of atoms that must be placed
   * @param previousPlacement is the previous placement of the atoms
   * @param discreteRows is a map from the previous rows to their discrete
   *      indices
   * @param discreteColumns is a map from the previous columns to their discrete
   *    indices
   * @param rowMapping is a map from the discrete rows to their indices
   * @param columnMapping is a map from the discrete columns to their indices
   * @param discreteTargetRows is a map from the target rows to their discrete
   *  indices
   *  @param discreteTargetColumns is a map from the target columns to their
   *   discrete indices
   *  @param placement is the updated placement of the atoms
   */
  auto updatePlacement(
      const std::vector<qc::Qubit>& atoms,
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousPlacement,
      const std::unordered_map<
          std::pair<std::reference_wrapper<const SLM>, size_t>, uint16_t,
          std::hash<std::pair<const SLM&, size_t>>,
          std::equal_to<std::pair<const SLM&, size_t>>>& discreteRows,
      const std::unordered_map<
          std::pair<std::reference_wrapper<const SLM>, size_t>, uint16_t,
          std::hash<std::pair<const SLM&, size_t>>,
          std::equal_to<std::pair<const SLM&, size_t>>>& discreteColumns,
      const std::unordered_map<uint16_t, uint16_t>& rowMapping,
      const std::unordered_map<uint16_t, uint16_t>& columnMapping,
      const std::unordered_map<
          std::pair<std::reference_wrapper<const SLM>, size_t>, uint16_t,
          std::hash<std::pair<const SLM&, size_t>>,
          std::equal_to<std::pair<const SLM&, size_t>>>& discreteTargetRows,
      const std::unordered_map<
          std::pair<std::reference_wrapper<const SLM>, size_t>, uint16_t,
          std::hash<std::pair<const SLM&, size_t>>,
          std::equal_to<std::pair<const SLM&, size_t>>>& discreteTargetColumns,
      std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                             size_t>>& placement) const -> void;

  /// This function generates a trivial initial placement for the qubits and
  /// just places in order.
  [[nodiscard]] auto makeInitialPlacement(size_t nQubits) const -> std::vector<
      std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

  [[nodiscard]] auto makeIntermediatePlacement(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousPlacement,
      const std::unordered_set<qc::Qubit>& previousReuseQubits,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& twoQubitGates)
      -> std::pair<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                          size_t, size_t>>,
                   std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                          size_t, size_t>>>;

  /**
   * This function places the qubits corresponding to gates in the entanglement
   * zone. After this placement has been performed, the activation of the
   * Rydberg beam will execute the gates in the given layer. Afterward, the
   * next placement for moving (non-reuse) qubits back to the storage zone is
   * determined by @ref placeQubitsInStorageZone.
   */
  [[nodiscard]] auto placeGatesInEntanglementZone(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousPlacement,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& twoQubitGates)
      -> std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

  /**
   * This function places qubits from the entanglement zone in the storage
   * zone after a rydberg gate has been performed.
   * It initializes the graph structure for the A* algorithm.
   * Afterward, the A* algorithm is called to find the optimal mapping.
   */
  auto placeQubitsInStorageZone(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousPlacement,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::pair<qc::Qubit, qc::Qubit>>& twoQubitGates)
      -> std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

  auto isGoal(const Node& node) const -> bool {
    return node.consumedFreeSites.size() == nJobs_;
  }

  /// @brief Returns the cost of a node, i.e., the total cost to reach that node
  /// from the start node.
  /// @details Different groups cannot be rearranged concurrently in one step.
  /// Hence, we add up the time it takes to perform the rearrangements of one
  /// group in one step and sum it up over all groups.
  /// This will not resemble the exact time to rearrange all costs because at
  /// this point it is not clear how the horizontal and vertical groups can be
  /// combined.
  [[nodiscard]] auto getCost(const Node& node) const -> float;

  /// @brief Return the estimated cost still required to reach a goal node.
  /// @details To yield an optimal results, the heuristic must be admissible,
  /// i.e., never overestimating the cost.
  /// The heuristic returns the estimated costs that are still added to the
  /// current actual cost to reach a goal node.
  /// Hence, the heuristic must always be less or equal to the additional cost
  /// needed to reach a goal.
  /// In the best case, all atoms that are not placed yet are compatible with
  /// an existing group and can just be added to that group.
  /// Hence, the sum in the cost function does not get an additional summand
  /// just the existing summands may increase.
  /// In the case of minimal increase in the overall cost, only one summand
  /// increases its value.
  /// This increase is bounded from below by the maximal distance of an atom to
  /// its nearest potential target site minus the maximum distance already
  /// placed atoms must travel to their determined target site.
  [[nodiscard]] auto getAtomPlacementHeuristic(const Node& node) const -> float;

  /// @brief Return the estimated cost still required to reach a goal node.
  /// @details To yield an optimal results, the heuristic must be admissible,
  /// i.e., never overestimating the cost.
  /// The heuristic returns the estimated costs that are still added to the
  /// current actual cost to reach a goal node.
  /// Hence, the heuristic must always be less or equal to the additional cost
  /// needed to reach a goal.
  /// In the best case, all atoms that are not placed yet are compatible with
  /// an existing group and can just be added to that group.
  /// Hence, the sum in the cost function does not get an additional summand
  /// just the existing summands may increase.
  /// In the case of minimal increase in the overall cost, only one summand
  /// increases its value.
  /// This increase is bounded from below by the maximal distance of an atom to
  /// its nearest potential target site minus the maximum distance already
  /// placed atoms must travel to their determined target site.
  [[nodiscard]] auto getGatePlacementHeuristic(const Node& node) const -> float;

  /// @brief Return pointers to all neighbors of the given node.
  /// @details When calling this function, the neighbors are allocated
  /// permanently such that (1) the returned pointers remain valid when the
  /// execution returned from this function and (2) not all nodes in the tree
  /// have to be created before they are needed.
  /// Hence, nodes are only created on demand in this function.
  /// Consequently, this function must only be called once per node.
  /// Otherwise, neighbors for the same node are created twice.
  /// @par
  /// When creating a new node, the horizontal and vertical groups are checked
  /// whether the new corresponding placement is compatible with any of the
  /// existing groups.
  /// If yes, the new placement is added to the respective group and otherwise,
  /// a new group is formed with the new placement.
  [[nodiscard]] auto getAtomPlacementNeighbors(const Node& node)
      -> std::vector<std::reference_wrapper<const Node>>;

  /// @brief Return pointers to all neighbors of the given node.
  /// @details When calling this function, the neighbors are allocated
  /// permanently such that (1) the returned pointers remain valid when the
  /// execution returned from this function and (2) not all nodes in the tree
  /// have to be created before they are needed.
  /// Hence, nodes are only created on demand in this function.
  /// Consequently, this function must only be called once per node.
  /// Otherwise, neighbors for the same node are created twice.
  /// @par
  /// When creating a new node, the horizontal and vertical groups are checked
  /// whether the new corresponding placement is compatible with any of the
  /// existing groups.
  /// If yes, the new placement is added to the respective group and otherwise,
  /// a new group is formed with the new placement.
  [[nodiscard]] auto getGatePlacementNeighbors(const Node& node)
      -> std::vector<std::reference_wrapper<const Node>>;

  /// Checks for the new placement of the atom whether it is compatible with
  /// one of the existing groups. If yes, the new placement is added to the
  /// respective group. Otherwise, a new group is formed with the new placement.
  /// @param key the start index of the new placement
  /// @param value the target index of the new placement
  /// @param groups the groups to which the new placement can be added
  /// @param maxDistances the maximum distances of placed atoms in each group
  /// @return true if the new placement could be added to an existing group
  static auto checkCompatibilityAndAddPlacement(
      uint16_t key, uint16_t value, float distance,
      std::vector<std::map<uint16_t, uint16_t>>& groups,
      std::vector<float>& maxDistances) -> bool;
};
} // namespace na
