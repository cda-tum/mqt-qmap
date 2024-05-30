//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "HybridNeutralAtomMapper.hpp"
#include "NeutralAtomArchitecture.hpp"
#include "NeutralAtomUtils.hpp"

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
};
} // namespace qc
