#pragma once

#include "ir/QuantumComputation.hpp"
#include "na/azac/Architecture.hpp"
#include "na/azac/Types.hpp"

#include <functional>
#include <nlohmann/json_fwd.hpp>
#include <utility>
#include <vector>

namespace na::azac {
/**
 * The class ASAPScheduler implements the as-soon-as-possible scheduling
 * strategy for the zoned neutral atom compiler.
 */
class ASAPScheduler {
  std::reference_wrapper<const Architecture> architecture_;
  size_t maxTwoQubitGateNumPerLayer_ = 0;

public:
  /**
   * Create a new ASAPScheduler.
   * @note The second parameter of the constructor is unused.
   * @param architecture is the architecture of the neutral atom system
   */
  ASAPScheduler(const Architecture& architecture, const nlohmann::json& config);
  /**
   * This function schedules the operations of a quantum computation.
   * @details Every operation is scheduled as soon as possible. The function
   * splits the operations into layers. Every layer (except for the last one)
   * contains some one-qubit operations and two-qubit operations. The one-qubit
   * operations are executed before the two-qubit operations. For every layer,
   * all two-qubit operations can be executed in parallel, i.e., every qubit is
   * involved in at most one two-qubit operation. The last layer contains only
   * the remaining one-qubit operations.
   * @param qc is the quantum computation
   * @return a pair of two vectors. The first vector contains the layers of
   * one-qubit operations. The second vector contains the layers of two-qubit
   * operations. A pair of qubits represents every two-qubit operation.
   */
  [[nodiscard]] auto schedule(const qc::QuantumComputation& qc) const
      -> std::pair<std::vector<OneQubitGateLayer>,
                   std::vector<TwoQubitGateLayer>>;
};
} // namespace na::azac
