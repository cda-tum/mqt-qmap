//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "NeutralAtomLayer.hpp"
#include "hybridmap/HardwareQubits.hpp"
#include "hybridmap/Mapping.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "hybridmap/NeutralAtomScheduler.hpp"
#include "hybridmap/NeutralAtomUtils.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/Operation.hpp"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace na {

/**
 * @brief Struct to store the runtime parameters of the mapper.
 */
struct MapperParameters {
  qc::fp lookaheadWeightSwaps = 0.1;
  qc::fp lookaheadWeightMoves = 0.1;
  qc::fp decay = 0.1;
  qc::fp shuttlingTimeWeight = 1;
  qc::fp gateWeight = 1;
  qc::fp shuttlingWeight = 1;
  uint32_t seed = 0;
  bool verbose = false;
  InitialCoordinateMapping initialMapping = InitialCoordinateMapping::Trivial;
};

/**
 * @brief Class to map a quantum circuit to a neutral atom architecture.
 * @details The mapping has following important parts:
 * - initial mapping: The initial mapping of the circuit qubits to the hardware
 * qubits.
 * - layer creation: The creation of the front and lookahead layers, done one
 * the fly and taking into account basic commutation rules.
 * - estimation: The estimation of the number of swap gates and moves needed to
 * execute a given gate and the decision which technique is better.
 * - gate based mapping: SABRE based algorithm to choose the bast swap for the
 * given layers.
 * - shuttling based mapping: Computing and evaluation of possible moves and
 * choosing best.
 * - multi-qubit-gates: Additional steps and checks to bring multiple qubits
 * together.
 * -> Final circuit contains abstract SWAP gates and MOVE operations, which need
 * to be decomposed using AODScheduler.
 */
class NeutralAtomMapper {
protected:
  // The considered architecture
  const NeutralAtomArchitecture& arch;
  // The mapped quantum circuit
  qc::QuantumComputation mappedQc;
  // The mapped quantum circuit converted to AOD movements
  qc::QuantumComputation mappedQcAOD;
  // The scheduler to schedule the mapped quantum circuit
  NeutralAtomScheduler scheduler;
  // The gates that have been executed
  std::vector<const qc::Operation*> executedCommutingGates;
  // Gates in the front layer to be executed with swap gates
  GateList frontLayerGate;
  // Gates in the front layer to be executed with move operations
  GateList frontLayerShuttling;
  // Gates in the lookahead layer to be executed with swap gates
  GateList lookaheadLayerGate;
  // Gates in the lookahead layer to be executed with move operations
  GateList lookaheadLayerShuttling;
  // The minimal weight for any multi-qubit gate
  qc::fp twoQubitSwapWeight = 1;
  // The runtime parameters of the mapper
  MapperParameters parameters;
  // The qubits that are blocked by the last swap
  std::deque<std::set<HwQubit>> lastBlockedQubits;
  // The last moves that have been executed
  std::deque<AtomMove> lastMoves;
  // Precomputed decay weights
  std::vector<qc::fp> decayWeights;
  // Counter variables
  uint32_t nSwaps = 0;
  uint32_t nMoves = 0;

  // The current placement of the hardware qubits onto the coordinates
  HardwareQubits hardwareQubits;
  // The current mapping between circuit qubits and hardware qubits
  Mapping mapping;

  // Methods for mapping
  /**
   * @brief Maps the gate to the mapped quantum circuit.
   * @param op The gate to map
   */
  void mapGate(const qc::Operation* op);
  /**
   * @brief Maps all currently possible gates and updated until no more gates
   * can be mapped.
   * @param layer The layer to map all possible gates for
   */
  void mapAllPossibleGates(NeutralAtomLayer& layer);
  /**
   * @brief Returns all gates that can be executed now
   * @param gates The gates to be checked
   * @return All gates that can be executed now
   */
  GateList getExecutableGates(const GateList& gates);
  /**
   * @brief Checks if the given gate can be executed for the given mapping and
   * hardware arrangement.
   * @param opPointer The gate to check
   * @return True if the gate can be executed, false otherwise
   */
  bool isExecutable(const qc::Operation* opPointer);

  /**
   * @brief Update the mapping for the given swap gate.
   * @param swap The swap gate to update the mapping for
   */
  void updateMappingSwap(Swap swap);
  /**
   * @brief Update the mapping for the given move operation.
   * @param move The move operation to update the mapping for
   */
  void updateMappingMove(AtomMove move);

  // Methods for gate vs. shuttling
  /**
   * @brief Assigns the given gates to the gate or shuttling layers.
   * @param frontGates The gates to be assigned to the front layers
   * @param lookaheadGates The gates to be assigned to the lookahead layers
   */
  void reassignGatesToLayers(const GateList& frontGates,
                             const GateList& lookaheadGates);
  /**
   * @brief Estimates the minimal number of swap gates and time needed to
   * execute the given gate.
   * @param opPointer The gate to estimate the number of swap gates and time for
   * @return The minimal number of swap gates and time needed to execute the
   * given gate
   */
  std::pair<uint32_t, qc::fp>
  estimateNumSwapGates(const qc::Operation* opPointer);
  /**
   * @brief Estimates the minimal number of move operations and time needed to
   * execute the given gate.
   * @param opPointer The gate to estimate the number of move operations and
   * time for
   * @return The minimal number of move operations and time needed to execute
   * the given gate
   */
  std::pair<uint32_t, qc::fp> estimateNumMove(const qc::Operation* opPointer);
  /**
   * @brief Uses estimateNumSwapGates and estimateNumMove to decide if a swap
   * gate or move operation is better.
   * @param opPointer The gate to estimate the number of swap gates and move
   * operations for
   * @return True if a swap gate is better, false if a move operation is better
   */
  bool swapGateBetter(const qc::Operation* opPointer);

  // Methods for swap gates mapping
  /**
   * @brief Finds the best swap gate for the front layer.
   * @details The best swap gate is the one that minimizes the cost function.
   * This takes into account close by swaps from two-qubit gates and exact moves
   * from multi-qubit gates.
   * @return The best swap gate for the front layer
   */
  Swap findBestSwap(const Swap& lastSwap);
  /**
   * @brief Returns all possible swap gates for the front layer.
   * @details The possible swap gates are all swaps starting from qubits in the
   * front layer.
   * @return All possible swap gates for the front layer
   */
  [[nodiscard]] std::set<Swap>
  getAllPossibleSwaps(const std::pair<Swaps, WeightedSwaps>& swapsFront) const;

  /**
   * @brief Returns the next best shuttling move operation for the front layer.
   * @return The next best shuttling move operation for the front layer
   */

  // Methods for shuttling operations mapping
  /**
   * @brief Finds the current best move operation based on the cost function.
   * @details Uses getAllMoveCombinations to find all possible move combinations
   * (direct move, move away, multi-qubit moves) and then chooses the best one
   * based on the cost function.
   * @return The current best move operation
   */
  AtomMove findBestAtomMove();
  /**
   * @brief Returns all possible move combinations for the front layer.
   * @details This includes direct moves, move away and multi-qubit moves.
   * Only move combinations with minimal number of moves are kept.
   * @return Vector of possible move combinations for the front layer
   */
  MoveCombs getAllMoveCombinations();
  /**
   * @brief Returns all possible move away combinations for a move from start to
   * target.
   * @details The possible move away combinations are all combinations of move
   * operations that move qubits away and then the performs the actual move
   * operation. The move away is chosen such that it is in the same direction as
   * the second move operation.
   * @param start The start position of the actual move operation
   * @param target The target position of the actual move operation
   * @param excludedCoords Coordinates the qubits should not be moved to
   * @return All possible move away combinations for a move from start to target
   */
  MoveCombs getMoveAwayCombinations(CoordIndex start, CoordIndex target,
                                    const CoordIndices& excludedCoords);

  // Helper methods
  /**
   * @brief Distinguishes between two-qubit swaps and multi-qubit swaps.
   * @details Two-qubit swaps only need to swap next to each other, while
   * multi-qubit swaps need to swap exactly to the multi-qubit gate position.
   * The multi-qubit swaps are weighted depending on their importance to finish
   * the multi-qubit gate.
   * @param layer The layer to distinguish the swaps for (front or lookahead)
   * @return The two-qubit swaps and multi-qubit swaps for the given layer
   */
  std::pair<Swaps, WeightedSwaps> initSwaps(const GateList& layer);
  /**
   * @brief Helper function to set the two-qubit swap weight to the minimal
   * weight of all multi-qubit gates, or 1.
   * @param swapExact The exact moves from multi-qubit gates
   */
  void setTwoQubitSwapWeight(const WeightedSwaps& swapExact);

  /**
   * @brief Returns the best position for the given gate coordinates.
   * @details Recursively calls getMovePositionRec
   * @param gateCoords The coordinates of the gate to find the best position for
   * @return The best position for the given gate coordinates
   */
  CoordIndices getBestMovePos(const CoordIndices& gateCoords);
  MultiQubitMovePos getMovePositionRec(MultiQubitMovePos currentPos,
                                       const CoordIndices& gateCoords,
                                       const size_t& maxNMoves);
  /**
   * @brief Returns possible move combinations to move the gate qubits to the
   * given position.
   * @param gateQubits The gate qubit to be moved
   * @param position The target position of the gate qubits
   * @return Possible move combinations to move the gate qubits to the given
   * position
   */
  MoveCombs getMoveCombinationsToPosition(HwQubits& gateQubits,
                                          CoordIndices& position);

  // Multi-qubit gate based methods
  /**
   * @brief Returns the best position for the given multi-qubit gate.
   * @details Calls getBestMultiQubitPositionRec to find the best position by
   * performing a recursive search in a breadth-first manner.
   * @param opPointer The multi-qubit gate to find the best position for
   * @return The best position for the given multi-qubit gate
   */
  HwQubits getBestMultiQubitPosition(const qc::Operation* opPointer);
  HwQubits getBestMultiQubitPositionRec(HwQubits remainingGateQubits,
                                        std::vector<HwQubit> selectedQubits,
                                        HwQubits remainingNearbyQubits);
  /**
   * @brief Returns the swaps needed to move the given qubits to the given
   * multi-qubit gate position.
   * @param op The multi-qubit gate to find the best position for
   * @param position The target position of the multi-qubit gate
   * @return The swaps needed to move the given qubits to the given multi-qubit
   */
  WeightedSwaps getExactSwapsToPosition(const qc::Operation* op,
                                        HwQubits position);

  // Cost function calculation
  /**
   * @brief Calculates the distance reduction for a swap gate given the
   * necessary close by swaps and exact moves.
   * @details Close by swaps are from two qubit gates, which only require to
   * swap close by. The exact moves are from multi-qubit gates, that require
   * swapping exactly to the multi-qubit gate position.
   * @param swap The swap gate to compute the distance reduction for
   * @param swapCloseBy The close by swaps from two-qubit gates
   * @param moveExact The exact moves from multi-qubit gates
   * @return The distance reduction cost
   */
  qc::fp swapCostPerLayer(const Swap& swap, const Swaps& swapCloseBy,
                          const WeightedSwaps& swapExact);
  /**
   * @brief Calculates the cost of a swap gate.
   * @details The cost of a swap gate is computed with the following terms:
   * - distance reduction for front + lookahead layers using swapCostPerLayer
   * - decay term for blocked qubit from last swaps
   * The cost is negative.
   * @param swap The swap gate to compute the cost for
   * @return The cost of the swap gate
   */
  qc::fp swapCost(const Swap& swap,
                  const std::pair<Swaps, WeightedSwaps>& swapsFront,
                  const std::pair<Swaps, WeightedSwaps>& swapsLookahead);
  /**
   * @brief Calculates the cost of a move operation.
   * @details Assumes the move is executed and computes the distance reduction
   * for the layer.
   * @param move The move operation to compute the cost for
   * @param layer The layer to compute the distance reduction for
   * @return The distance reduction cost
   */
  qc::fp moveCostPerLayer(const AtomMove& move, GateList& layer);

  /**
   * @brief Calculates a parallelization cost if the move operation can be
   * parallelized with the last moves.
   * @param move The move operation to compute the cost for
   * @return The parallelization cost
   */
  qc::fp parallelMoveCost(const AtomMove& move);
  /**
   * @brief Calculates the cost of a move operation.
   * @details The cost of a move operation is computed with the following terms:
   * - distance reduction for front + lookahead layers using moveCostPerLayer
   * - parallelization term based on last moves using parallelMoveCost
   * The three contributions are weighted with the runtime parameters.
   * @param move The move operation to compute the cost for
   * @return The cost of the move operation
   */
  qc::fp moveCost(const AtomMove& move);
  /**
   * @brief Calculates the cost of a series of move operations by summing up the
   * cost of each move.
   * @param moveComb The series of move operations to compute the cost for
   * @return The total cost of the series of move operations
   */
  qc::fp moveCostComb(const MoveComb& moveComb);

  /**
   * @brief Print the current layers for debugging.
   */
  void printLayers();

public:
  // Constructors
  [[maybe_unused]] NeutralAtomMapper(const NeutralAtomMapper&) = delete;
  NeutralAtomMapper& operator=(const NeutralAtomMapper&) = delete;
  NeutralAtomMapper(NeutralAtomMapper&&) = delete;
  explicit NeutralAtomMapper(const NeutralAtomArchitecture& architecture,
                             const MapperParameters& p = MapperParameters())
      : arch(architecture), mappedQc(architecture.getNpositions()),
        mappedQcAOD(architecture.getNpositions()), scheduler(architecture),
        parameters(p), hardwareQubits(architecture, parameters.initialMapping,
                                      parameters.seed) {
    // need at least on free coordinate to shuttle
    if (architecture.getNpositions() - architecture.getNqubits() < 1) {
      this->parameters.gateWeight = 1;
      this->parameters.shuttlingWeight = 0;
    }
  };

  /**
   * @brief Sets the runtime parameters of the mapper.
   * @param p The runtime parameters of the mapper
   */
  void setParameters(const MapperParameters& p) {
    this->parameters = p;
    if (arch.getNpositions() - arch.getNqubits() < 1) {
      this->parameters.gateWeight = 1;
      this->parameters.shuttlingWeight = 0;
    }
    this->reset();
  }

  /**
   * @brief Resets the mapper and the hardware qubits.
   */
  void reset() {
    hardwareQubits =
        HardwareQubits(arch, parameters.initialMapping, parameters.seed);
  }

  // Methods
  /**
   * @brief Maps the given quantum circuit to the given architecture.
   * @details The mapping has following important parts:
   * - initial mapping: The initial mapping of the circuit qubits to the
   * hardware qubits.
   * - layer creation: The creation of the front and lookahead layers, done one
   * the fly and taking into account basic commutation rules.
   * - estimation: The estimation of the number of swap gates and moves needed
   * to execute a given gate and the decision which technique is better.
   * - gate based mapping: SABRE based algorithm to choose the bast swap for the
   * given layers.
   * - shuttling based mapping: Computing and evaluation of possible moves and
   * choosing best.
   * - multi-qubit-gates: Additional steps and checks to bring multiple qubits
   * together.
   * -> Final circuit contains abstract SWAP gates and MOVE operations, which
   * need to be decomposed using convertToAod method.
   *
   * @param qc The quantum circuit to be mapped
   * @param initialMapping The initial mapping of the circuit qubits to the
   * hardware qubits
   * @param verbose If true, prints additional information
   * @return The mapped quantum circuit with abstract SWAP gates and MOVE
   * operations
   */
  qc::QuantumComputation map(qc::QuantumComputation& qc,
                             InitialMapping initialMapping);

  /**
   * @brief Maps the given quantum circuit to the given architecture and
   * converts it to the AOD level.
   * @param qc  The quantum circuit to be mapped
   * @param initialMapping The initial mapping of the circuit qubits to the
   * hardware qubits
   */
  [[maybe_unused]] void mapAndConvert(qc::QuantumComputation& qc,
                                      InitialMapping initialMapping,
                                      bool printInfo) {
    this->parameters.verbose = printInfo;
    map(qc, initialMapping);
    convertToAod(this->mappedQc);
  }

  /**
   * @brief Prints the mapped circuits as an extended OpenQASM string.
   * @return The mapped quantum circuit with abstract SWAP gates and MOVE
   */
  [[maybe_unused]] std::string getMappedQc() {
    std::stringstream ss;
    this->mappedQc.dumpOpenQASM(ss, false);
    return ss.str();
  }

  /**
   * @brief Saves the mapped quantum circuit to a file.
   * @param filename The name of the file to save the mapped quantum circuit to
   */
  [[maybe_unused]] void saveMappedQc(const std::string& filename) {
    std::ofstream ofs(filename);
    this->mappedQc.dumpOpenQASM(ofs, false);
  }

  /**
   * @brief Prints the mapped circuit with AOD operations as an extended
   * OpenQASM
   * @return The mapped quantum circuit with native AOD operations
   */
  [[maybe_unused]] std::string getMappedQcAOD() {
    std::stringstream ss;
    this->mappedQcAOD.dumpOpenQASM(ss, false);
    return ss.str();
  }

  /**
   * @brief Saves the mapped quantum circuit with AOD operations to a file.
   * @param filename The name of the file to save the mapped quantum circuit
   * with AOD operations to
   */
  [[maybe_unused]] void saveMappedQcAOD(const std::string& filename) {
    std::ofstream ofs(filename);
    this->mappedQcAOD.dumpOpenQASM(ofs, false);
  }

  /**
   * @brief Schedules the mapped quantum circuit on the neutral atom
   * architecture.
   * @details For each gate/operation in the input circuit, the scheduler checks
   * the earliest possible time slot for execution. If the gate is a multi qubit
   * gate, also the blocking of other qubits is taken into consideration. The
   * execution times are read from the neutral atom architecture.
   * @param verboseArg If true, prints additional information
   * @param createAnimationCsv If true, creates a csv file for the animation
   * @param shuttlingSpeedFactor The factor to speed up the shuttling time
   * @return The results of the scheduler
   */
  [[maybe_unused]] SchedulerResults
  schedule(bool verboseArg = false, bool createAnimationCsv = false,
           qc::fp shuttlingSpeedFactor = 1.0) {
    return scheduler.schedule(mappedQcAOD, hardwareQubits.getInitHwPos(),
                              verboseArg, createAnimationCsv,
                              shuttlingSpeedFactor);
  }

  /**
   * @brief Saves the animation csv file of the scheduled quantum circuit.
   * @return The animation csv string
   */
  [[maybe_unused]] std::string getAnimationCsv() {
    return scheduler.getAnimationCsv();
  }

  /**
   * @brief Saves the animation csv file of the scheduled quantum circuit.
   * @param filename The name of the file to save the animation csv file to
   */
  [[maybe_unused]] void saveAnimationCsv(const std::string& filename) {
    scheduler.saveAnimationCsv(filename);
  }

  /**
   * @brief Converts a mapped circuit down to the AOD level and CZ level.
   * @details SWAP gates are decomposed into CX gates. Then CnX gates are
   * decomposed into CnZ gates. Move operations are combined if possible and
   * then converted into native AOD operations.
   * @param qc The already mapped quantum circuit with abstract SWAP gates and
   * MOVE operations
   * @return The mapped quantum circuit with native AOD operations
   */
  qc::QuantumComputation convertToAod(qc::QuantumComputation& qc);

  [[maybe_unused]] [[nodiscard]] std::map<HwQubit, HwQubit>
  getInitHwPos() const {
    return hardwareQubits.getInitHwPos();
  }
};

} // namespace na
