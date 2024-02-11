//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "NeutralAtomMapper.hpp"

namespace na {
auto NeutralAtomMapper::map(const qc::QuantumComputation& qc) -> void {
  initialQc = qc;
  // get start time
  auto start = std::chrono::high_resolution_clock::now();
  //============================ START MAPPING =================================
  mappedQc = qc;
  //============================= END MAPPING ==================================
  // get end time
  auto end = std::chrono::high_resolution_clock::now();
  // build remaining statistics
  stats.numInitialGates = qc.getNops();
  stats.numMappedGates  = qc.getNops();
  // get the mapping time in milliseconds
  stats.mappingTime =
      std::chrono::duration<qc::fp, std::milli>(end - start).count();
}
} // namespace na