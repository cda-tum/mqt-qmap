//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "cstdint"
#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "operations/AodOperation.hpp"
#include "operations/OpType.hpp"
#include "utility"
#include "vector"

#include <memory>

namespace qc {
// Possible types two Move combination can be combined to
enum class ActivationMergeType : uint8_t { Impossible, Trivial, Merge, Append };

/**
 * @brief Class to convert abstract move operations to AOD movements on a
 * neutral atom architecture
 * @details The scheduler takes a quantum circuit containing abstract move
 * operations and tries to merge them into parallel AO movements. It also
 * manages the small offset movements required while loading or unloading of
 * AODs.
 */
class MoveToAodConverter {
protected:
  /**
   * @brief Struct to store information about specific AOD activations.
   * @details It contains:
   * - the offset moves in x and y direction
   * - the actual moves
   */
  struct AodActivationHelper {
    /**
     * @brief Describes a single AOD movement in either x or y direction
     * @details It contains:
     * - the initial position of the AOD
     * - the delta of the AOD movement
     * - the size of the offset of the AOD movement
     */
    struct AodMove {
      // start of the move
      uint32_t init;
      // need offset move to avoid crossing
      int32_t offset;
      // delta of the actual move
      fp delta;

      AodMove() = default;

      AodMove(uint32_t init, fp delta, int32_t offset)
          : init(init), offset(offset), delta(delta) {}
    };
    /**
     * @brief Manages the activation of an atom using an AOD.
     * @details The same struct is used also to deactivate the AOD, just
     * reversed.
     */
    struct AodActivation {
      // first: x, second: delta x, third: offset x
      std::vector<std::shared_ptr<AodMove>> activateXs;
      std::vector<std::shared_ptr<AodMove>> activateYs;
      std::vector<AtomMove>                 moves;

      AodActivation(const AodMove& activateX, const AodMove& activateY,
                    const AtomMove& move)
          : moves({move}) {
        activateXs.push_back(std::make_unique<AodMove>(activateX));
        activateYs.push_back(std::make_unique<AodMove>(activateY));
      }
      AodActivation(const Dimension dim, const AodMove& activate,
                    const AtomMove& move)
          : moves({move}) {
        if (dim == Dimension::X) {
          activateXs.push_back(std::make_unique<AodMove>(activate));
        } else {
          activateYs.push_back(std::make_unique<AodMove>(activate));
        }
      }

      [[nodiscard]] std::vector<std::shared_ptr<AodMove>>
      getActivates(Dimension dim) const {
        if (dim == Dimension::X) {
          return activateXs;
        }
        return activateYs;
      }
    };

    // AODScheduler
    // NeutralAtomArchitecture to call necessary hardware information
    const NeutralAtomArchitecture& arch;
    std::vector<AodActivation>     allActivations;
    // Differentiate between loading and unloading
    OpType type;

    // Constructor
    AodActivationHelper()                           = delete;
    AodActivationHelper(const AodActivationHelper&) = delete;
    AodActivationHelper(AodActivationHelper&&)      = delete;
    AodActivationHelper(const NeutralAtomArchitecture& arch, OpType type)
        : arch(arch), type(type) {}

    // Methods

    /**
     * @brief Returns all AOD moves in the given dimension/direction which start
     * at the given initial position
     * @param dim The dimension/direction to check
     * @param init The initial position to check
     * @return A vector of AOD moves
     */
    [[nodiscard]] std::vector<std::shared_ptr<AodMove>>
    getAodMovesFromInit(Dimension dim, uint32_t init) const;

    // Activation management
    /**
     * @brief Checks if the move can be added to the current activations
     * @param origin The origin of the move
     * @param v The move vector of the move
     * @return A pair of ActivationMerge, in x and y direction
     */
    [[nodiscard]] std::pair<ActivationMergeType, ActivationMergeType>
    canAddActivation(const Coordinate& origin, MoveVector v) const;
    /**
     * @brief Checks if the move can be added to the current activations in the
     * given dimension/direction
     * @param dim The dimension/direction to check
     * @param origin The origin of the move
     * @param v The move vector of the move
     * @return The ActivationMerge type
     */
    [[nodiscard]] ActivationMergeType
    canAddActivationDim(Dimension dim, const Coordinate& origin,
                        MoveVector v) const;
    /**
     * @brief Adds the move to the current activations
     * @details The move is merged into the current activations depending on the
     * given merge types
     * @param merge The merge types in x and y direction
     * @param origin The origin of the move
     * @param move The move to add
     * @param v The move vector of the move
     */
    void
    addActivation(std::pair<ActivationMergeType, ActivationMergeType> merge,
                  const Coordinate& origin, const AtomMove& move, MoveVector v);
    /**
     * @brief Merges the given activation into the current activations
     * @param dim The dimension/direction of the activation
     * @param activationDim The activation to merge in the given
     * dimension/direction
     * @param activationOtherDim The activation to merge/add in the other
     * dimension/direction
     */
    void mergeActivationDim(Dimension dim, const AodActivation& activationDim,
                            const AodActivation& activationOtherDim);
    /**
     * @brief Orders the aod offset moves such that they will not cross each
     * other
     * @param aodMoves The aod offset moves to order
     * @param sign The direction of the offset moves (right/left or down/up)
     */
    static void reAssignOffsets(std::vector<std::shared_ptr<AodMove>>& aodMoves,
                                int32_t                                sign);

    /**
     * @brief Returns the maximum offset in the given dimension/direction from
     * the given initial position
     * @param dim The dimension/direction to check
     * @param init The initial position to check
     * @param sign The direction of the offset moves (right/left or down/up)
     * @return The maximum offset
     */
    [[nodiscard]] uint32_t getMaxOffsetAtInit(Dimension dim, uint32_t init,
                                              int32_t sign) const;

    /**
     * @brief Checks if there is still space at the given initial position and
     * the given direction
     * @param dim The dimension/direction to check
     * @param init The initial position to check
     * @param sign The direction of the offset moves (right/left or down/up)
     * @return True if there is still space, false otherwise
     */
    [[nodiscard]] bool checkIntermediateSpaceAtInit(Dimension dim,
                                                    uint32_t  init,
                                                    int32_t   sign) const;

    // Convert activation to AOD operations
    /**
     * @brief Converts activation into AOD operation (activate, move,
     * deactivate)
     * @param activation The activation to convert
     * @param arch The neutral atom architecture to call necessary hardware
     * information
     * @param type The type of the activation (loading or unloading)
     * @return The activation as AOD operation
     */
    [[nodiscard]] std::pair<AodOperation, AodOperation>
    getAodOperation(const AodActivation& activation) const;
    /**
     * @brief Converts all activations into AOD operations
     * @return All activations of the AOD activation helper as AOD operations
     */
    [[nodiscard]] std::vector<AodOperation> getAodOperations() const;
  };

  /**
   * @brief Move operations within a move group can be executed in parallel
   * @details A move group contains:
   * - the moves that can be executed in parallel
   * - the AOD operations to load, shuttle and unload the atoms
   * - the qubits that are used by the gates in the move group
   */
  struct MoveGroup {
    // the moves and the index they appear in the original quantum circuit (to
    // insert them back later)
    std::vector<std::pair<AtomMove, uint32_t>> moves;
    std::vector<AodOperation>                  processedOpsInit;
    std::vector<AodOperation>                  processedOpsFinal;
    AodOperation                               processedOpShuttle;
    std::vector<CoordIndex>                    qubitsUsedByGates;

    // Constructor
    explicit MoveGroup() = default;

    // Methods
    /**
     * @brief Checks if the given move can be added to the move group
     * @param move Move to check
     * @return True if the move can be added, false otherwise
     */
    bool canAdd(const AtomMove& move, const NeutralAtomArchitecture& archArg);
    /**
     * @brief Adds the given move to the move group
     * @param move Move to add
     * @param idx Index of the move in the original quantum circuit
     */
    void add(const AtomMove& move, uint32_t idx);
    /**
     * @brief Returns the circuit index of the first move in the move group
     * @return Circuit index of the first move in the move group
     */
    [[nodiscard]] uint32_t getFirstIdx() const { return moves.front().second; }
    /**
     * @brief Checks if the two moves can be executed in parallel
     * @param v1 The first move
     * @param v2 The second move
     * @return True if the moves can be executed in parallel, false otherwise
     */
    static bool parallelCheck(const MoveVector& v1, const MoveVector& v2);

    /**
     * @brief Helper function to create the actual shuttling operation between
     * the loading at the initial position and the unloading at the final
     * position
     * @param opsInit Loading operations
     * @param opsFinal Unloading operations
     * @return The shuttling operation between the loading and unloading
     * operations
     */
    static AodOperation
    connectAodOperations(const std::vector<AodOperation>& opsInit,
                         const std::vector<AodOperation>& opsFinal);
  };

  const NeutralAtomArchitecture& arch;
  QuantumComputation             qcScheduled;
  std::vector<MoveGroup>         moveGroups;

  /**
   * @brief Assigns move operations into groups that can be executed in parallel
   * @param qc Quantum circuit to schedule
   */
  void initMoveGroups(QuantumComputation& qc);
  /**
   * @brief Converts the move groups into the actual AOD operations
   * @details For this the following steps are performed:
   * - ActivationHelper to manage the loading
   * - ActivationHelper to manage the unloading
   * If not the whole move group can be executed in parallel, a new move group
   * is created for the remaining moves.
   */
  void processMoveGroups();

public:
  MoveToAodConverter()                          = delete;
  MoveToAodConverter(const MoveToAodConverter&) = delete;
  MoveToAodConverter(MoveToAodConverter&&)      = delete;
  explicit MoveToAodConverter(const NeutralAtomArchitecture& archArg)
      : arch(archArg), qcScheduled(arch.getNpositions()) {}

  /**
   * @brief Schedules the given quantum circuit using AODs
   * @param qc Quantum circuit to schedule
   * @return Scheduled quantum circuit, containing AOD operations
   */
  QuantumComputation schedule(QuantumComputation& qc);

  /**
   * @brief Returns the number of move groups
   * @return Number of move groups
   */
  [[nodiscard]] auto getNMoveGroups() const { return moveGroups.size(); }
};

} // namespace qc
