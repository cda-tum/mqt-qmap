//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "HybridNeutralAtomMapper.hpp"
#include "NeutralAtomArchitecture.hpp"
#include "NeutralAtomUtils.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"

#include <cstdint>
#include <vector>

namespace qc {

/**
 * @brief Class to manage information exchange between the neutral atom mapper
 * and the ZX extraction.
 * @deteils This class is derived from the HybridNeutralAtomMapper and stores
 * all the information about the neutral atom hardware and the current status of
 * the mapping. The ZX (or another synthesis algorithm) can propose different
 * possible next synthesis steps, which are then evaluated by the
 * HybridNeutralAtomMapper regarding the "effort" to map this synthesis step. It
 * also provides additional functionality to exchange information between the ZX
 * and the HybridNeutralAtomMapper.
 *
 */
class HybridSynthesisMapper : private NeutralAtomMapper {
  using qcs = std::vector<QuantumComputation>;

  AdjacencyMatrix adjacencyMatrix;

public:
  // Constructors
  HybridSynthesisMapper() = delete;
  explicit HybridSynthesisMapper(const NeutralAtomMapper& neutralAtomMapper)
      : NeutralAtomMapper(neutralAtomMapper) {}
  explicit HybridSynthesisMapper(
      const NeutralAtomArchitecture& arch,
      InitialCoordinateMapping       initialCoordinateMapping =
          InitialCoordinateMapping::Trivial,
      const MapperParameters& params = MapperParameters())
      : NeutralAtomMapper(arch, initialCoordinateMapping, params) {}

  // Functions

  /**
   * @brief Remaps the whole circuit again starting from the initial mapping.
   * @param initialMapping The initial mapping to be used.
   * @return The remapped circuit.
   */
  QuantumComputation
  completelyRemap(InitialMapping initialMapping = InitialMapping::Identity);

  /**
   * @brief Evaluates the synthesis steps proposed by the ZX extraction.
   * @param synthesisSteps The synthesis steps proposed by the ZX extraction.
   * @param directlyMap If true, the synthesis steps are directly mapped to the
   * hardware.
   * @return The index of the synthesis step with the lowest effort.
   */
  uint32_t evaluateSynthesisSteps(qcs& synthesisSteps,
                                  bool directlyMap = false);

  /**
   * @brief Directly maps the given QuantumComputation to the hardware by
   * inserting SWAP gates or shuttling move operations.
   * @param qc The gates (QuantumComputation) to be mapped.
   */
  void directlyMap(const QuantumComputation& qc);

  /**
   * @brief Returns the current adjacency matrix of the neutral atom hardware.
   * @return The current adjacency matrix of the neutral atom hardware.
   */
  [[nodiscard]] AdjacencyMatrix getAdjacencyMatrix() const {
    return adjacencyMatrix;
  }
};
} // namespace qc
