#pragma once

#include "Definitions.hpp"
#include "na/azac/Architecture.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::azac {
class AStarPlacer {
  friend class AStarPlacerTest_AStarSearch_Test;
  using DiscreteSite = std::array<uint8_t, 2>;

  std::reference_wrapper<const Architecture> architecture_;
  /// If true, during the initial placement the atoms are placed starting in the
  /// last row instead of the first row in the first SLM
  bool reverseInitialPlacement_ = false;
  /// this flag indicates whether the  placement should use a window when
  /// selecting potential free sites
  bool useWindow_ = true;
  /// If the window is used, this denotes the minimum width in terms of rows of
  /// the window centered at the nearest site
  size_t windowMinWidth_ = 4;
  size_t windowMinHeight_ = 6;
  /// If the window is used, this denotes the ration between the height and the
  /// width of the window. A value greater than 1 means that the window is
  /// higher than wide. A value of exactly 1 means that the window is square. A
  /// value smaller than 1 means that the window is wider than high.
  double windowRatio_ = 1.5;
  /// If the window is used, this denotes the share of free sites in the window
  /// in relation to the number of atoms to be moved in this step. The window is
  /// extended according to the ratio as long as the share of free sites is
  /// smaller than this value. A value greater or equal to 1 ensures that there
  /// exists a solution. However, a smaller value might be a reasonable good
  /// guess since it is almost certain that not all atoms to be moved will end
  /// in the same window.
  double windowShare_ = 0.4;
  /// The heuristic used in the A* search contains a term that resembles the
  /// standard deviation of the differences between the current and target sites
  /// of the atoms to be moved in every orientation. This factor is multiplied
  /// with the standard deviation to adjust the influence of this term. Setting
  /// it to 0.0 disables this term resulting in an admissible heuristic.
  /// However, this leads to a vast exploration of the search tree and usually
  /// results in a huge number of nodes visited.
  float deepeningFactor_ = 1.0F;
  /// The cost function can consider the distance of atoms to their interaction
  /// partner in the next layer. This factor is multiplied with the distance to
  /// adjust the influence of this term. Setting it to 0.0 disables the
  /// lookahead entirely. A factor of 1.0 implies that the lookahead is as
  /// important as the distance to the target site, which is usually not
  /// desired.
  float lookaheadFactor_ = 0.2F;
  /// When placing atoms after a rydberg layer back in the storage zone, this
  /// struct stores for every such atom all required information, i.e., the
  /// current site and potential target sites ordered by distance (ascending).
  struct AtomJob {
    qc::Qubit qubit;
    /// the current site of the atom
    DiscreteSite currentSite;
    /// a struct describing one potential target site
    struct Option {
      /// the target site
      DiscreteSite site;
      /// the distance the atom must travel to reach the target site
      float distance;
      /// additional lookahead distance to next interaction partner
      float lookaheadCost = 0.0F;
    };
    /// a list of all potential target sites ordered by distance (ascending)
    std::vector<Option> options;
    /// minimum lookahead distance
    float minLookaheadCost = 0.0F;
  };

  /// When placing gates in the entanglement zone before a rydberg layer, this
  /// struct stores for every such gate all required information, i.e., the
  /// current sites of the corresponding atoms and potential target sites
  /// ordered by distance (ascending).
  struct GateJob {
    std::array<qc::Qubit, 2> qubits;
    /// the current sites of the two atoms
    std::array<DiscreteSite, 2> currentSites;
    /// a struct describing one potential target site for each atom
    struct Option {
      /// the target sites for the two atoms
      std::array<DiscreteSite, 2> sites;
      /// the max distance the atoms must travel to reach the target sites
      std::array<float, 2> distance;
      /// additional lookahead distance to next interaction partner
      float lookaheadCost = 0.0F;
    };
    /// a list of all potential target sites ordered by distance (ascending)
    std::vector<Option> options;
    /// minimum lookahead distance
    float minLookaheadCost = 0.0F;
  };

  /// A node representing one stage in the process of placing all atoms
  /// that must be moved for the next stage starting from the last mapping
  /// until a new mapping is found satisfying all constraints of the next
  /// stage
  struct AtomNode {
    const AtomJob::Option* option = nullptr;
    /// a set of all sites that are already occupied by an atom due to the
    /// current placement
    std::unordered_set<DiscreteSite> consumedFreeSites;
    /// a binary search tree representing the horizontal and vertical group,
    /// respectively
    /// @see getNeighbors for more details
    std::vector<std::array<std::map<uint8_t, uint8_t>, 2>> groups;
    /// the maximum distance of placed atoms in every group to their
    /// target location
    std::vector<float> maxDistancesOfPlacedAtomsPerGroup;
    /// accumulated lookahead cost
    float lookaheadCost = 0.0F;
  };

  /// A node representing one stage in the process of placing all atoms
  /// that must be moved for the next stage starting from the last mapping
  /// until a new mapping is found satisfying all constraints of the next
  /// stage
  struct GateNode {
    const GateJob::Option* option = nullptr;
    /// a set of all sites that are already occupied by an atom due to the
    /// current placement
    std::unordered_set<DiscreteSite> consumedFreeSites;
    /// a binary search tree representing the horizontal and vertical group,
    /// respectively
    /// @see getNeighbors for more details
    std::vector<std::array<std::map<uint8_t, uint8_t>, 2>> groups;
    /// the maximum distance of placed atoms in every group to their
    /// target location
    std::vector<float> maxDistancesOfPlacedAtomsPerGroup;
    /// accumulated lookahead cost
    float lookaheadCost = 0.0F;
  };

public:
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
        const std::vector<std::vector<std::array<qc::Qubit, 2>>>&
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
  template <class Node>
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
                       uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>,
                   std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint8_t, std::hash<std::pair<const SLM&, size_t>>,
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
                       uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>,
                   std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint8_t, std::hash<std::pair<const SLM&, size_t>>,
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
                       uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>,
                   std::unordered_map<
                       std::pair<std::reference_wrapper<const SLM>, size_t>,
                       uint8_t, std::hash<std::pair<const SLM&, size_t>>,
                       std::equal_to<std::pair<const SLM&, size_t>>>>;

  /// This function generates a trivial initial placement for the qubits and
  /// just places in order.
  [[nodiscard]] auto makeInitialPlacement(size_t nQubits) const -> std::vector<
      std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

  [[nodiscard]] auto makeIntermediatePlacement(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousPlacement,
      const std::unordered_set<qc::Qubit>& previousReuseQubits,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::array<qc::Qubit, 2>>& twoQubitGates,
      const std::vector<std::array<qc::Qubit, 2>>& nextTwoQubitGates)
      -> std::pair<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                          size_t, size_t>>,
                   std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                          size_t, size_t>>>;
  auto addGateOption(
      const std::unordered_map<
          std::pair<std::reference_wrapper<const SLM>, unsigned long>,
          unsigned char, std::hash<std::pair<const SLM&, unsigned long>>,
          std::equal_to<std::pair<const SLM&, unsigned long>>>&
          discreteTargetRows,
      const std::unordered_map<
          std::pair<std::reference_wrapper<const SLM>, unsigned long>,
          unsigned char, std::hash<std::pair<const SLM&, unsigned long>>,
          std::equal_to<std::pair<const SLM&, unsigned long>>>&
          discreteTargetColumns,
      const SLM& leftSlm, size_t leftRow, size_t leftCol, const SLM& rightSlm,
      size_t rightRow, size_t rightCol, const SLM& nearestSlm, size_t r,
      size_t c, GateJob& job) const -> void;

  /**
   * This function places the qubits corresponding to gates in the entanglement
   * zone. After this placement has been performed, the activation of the
   * Rydberg beam will execute the gates in the given layer. Afterward, the
   * next placement for moving (non-reuse) qubits back to the storage zone is
   * determined by @ref placeAtomsInStorageZone.
   */
  [[nodiscard]] auto placeGatesInEntanglementZone(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousPlacement,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::array<qc::Qubit, 2>>& twoQubitGates,
      const std::unordered_set<qc::Qubit>& nextReuseQubits,
      const std::vector<std::array<qc::Qubit, 2>>& nextTwoQubitGates)
      -> std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

  /**
   * This function places qubits from the entanglement zone in the storage
   * zone after a rydberg gate has been performed.
   * It initializes the graph structure for the A* algorithm.
   * Afterward, the A* algorithm is called to find the optimal mapping.
   */
  auto placeAtomsInStorageZone(
      const std::vector<std::tuple<std::reference_wrapper<const SLM>, size_t,
                                   size_t>>& previousPlacement,
      const std::unordered_set<qc::Qubit>& reuseQubits,
      const std::vector<std::array<qc::Qubit, 2>>& twoQubitGates,
      const std::vector<std::array<qc::Qubit, 2>>& nextTwoQubitGates)
      -> std::vector<
          std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>;

  template <class Node>
  [[nodiscard]] static auto isGoal(const size_t nAtoms, const Node& node)
      -> bool {
    return node.consumedFreeSites.size() == nAtoms;
  }

  /// @brief Returns the cost of a node, i.e., the total cost to reach that node
  /// from the start node.
  /// @details Different groups cannot be rearranged concurrently in one step.
  /// Hence, we add up the time it takes to perform the rearrangements of one
  /// group in one step and sum it up over all groups.
  /// This will not resemble the exact time to rearrange all costs because at
  /// this point it is not clear how the horizontal and vertical groups can be
  /// combined.
  [[nodiscard]] static auto getCost(const GateNode& node) -> float;
  [[nodiscard]] static auto getCost(const AtomNode& node) -> float;

  /// @brief Calculates the standard deviation of the differences value - key
  /// and sums them up over all horizontal and vertical groups.
  [[nodiscard]] static auto sumStdDeviationForGroups(
      const std::array<float, 2>& scaleFactors,
      const std::vector<std::array<std::map<uint8_t, uint8_t>, 2>>& groups)
      -> float;

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
  [[nodiscard]] static auto
  getHeuristic(const std::vector<AtomJob>& atomJobs, float deepeningFactor,
               const std::array<float, 2>& scaleFactors, const AtomNode& node)
      -> float;
  [[nodiscard]] static auto
  getHeuristic(const std::vector<GateJob>& gateJobs, float deepeningFactor,
               const std::array<float, 2>& scaleFactors, const GateNode& node)
      -> float;

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
  [[nodiscard]] static auto
  getNeighbors(std::deque<std::unique_ptr<AtomNode>>& nodes,
               const std::vector<AtomJob>& atomJobs, const AtomNode& node)
      -> std::vector<std::reference_wrapper<const AtomNode>>;
  [[nodiscard]] static auto
  getNeighbors(std::deque<std::unique_ptr<GateNode>>& nodes,
               const std::vector<GateJob>& gateJobs, const GateNode& node)
      -> std::vector<std::reference_wrapper<const GateNode>>;

  /// Checks the compatibility with a new assignment, i.e., a key-value pair,
  /// whether t is compatible with an existing group. The group can either be
  /// a horizontal or vertical group. In case, the new assignment is compatible
  /// with the group, an iterator is returned pointing to the assignment in the
  /// group, if it already exists or to the element directly following the
  /// new key. Additionally, a boolean is returned indicating whether the new
  /// exists in the group.
  /// @param key the key of the new assignment
  /// @param value the value of the new assignment
  /// @param group the group to which the new assignment is compared
  /// @return an optional pair of an iterator pointing to the key equal or
  /// directly following the new key in the group and a boolean indicating
  /// whether the new assignment exists in the group. If the new assignment is
  /// not compatible with the group, an empty optional is returned.
  [[nodiscard]] static auto
  checkCompatibilityWithGroup(uint8_t key, uint8_t value,
                              const std::map<uint8_t, uint8_t>& group)
      -> std::optional<
          std::pair<std::map<uint8_t, uint8_t>::const_iterator, bool>>;

  /// Checks for the new placement of the atom whether it is compatible with
  /// one of the existing groups. If yes, the new placement is added to the
  /// respective group. Otherwise, a new group is formed with the new placement.
  /// @param hKey the horizontal start index of the new placement
  /// @param hValue the horizontal target index of the new placement
  /// @param vKey the vertical start index of the new placement
  /// @param vValue the vertical target index of the new placement
  /// @param distance the distance the atom must travel to reach the target site
  /// @param groups the groups to which the new placement can be added
  /// @param maxDistances the maximum distances of placed atoms in each group
  /// @return true if the new placement could be added to an existing group
  static auto checkCompatibilityAndAddPlacement(
      uint8_t hKey, uint8_t hValue, uint8_t vKey, uint8_t vValue,
      float distance,
      std::vector<std::array<std::map<uint8_t, uint8_t>, 2>>& groups,
      std::vector<float>& maxDistances) -> bool;
};
} // namespace na::azac
