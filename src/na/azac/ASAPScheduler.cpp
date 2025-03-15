#include "na/azac/ASAPScheduler.hpp"

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/azac/Architecture.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iostream>
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace na {
ASAPScheduler::ASAPScheduler(const Architecture& architecture,
                             const nlohmann::json& config)
    : architecture_(architecture) {
  if (const auto& configIt = config.find("asap_scheduler");
      configIt != config.end() && configIt->is_object()) {
    for (const auto& [key, value] : configIt.value().items()) {
      std::ostringstream oss;
      oss << "[WARN] Configuration for ASAPScheduler contains an unknown key: "
          << key << ". Ignoring.\n";
      std::cout << oss.str();
    }
  }
  // calculate the maximum possible number of two-qubit gates per layer
  for (const auto& zone : architecture_.get().entanglementZones) {
    maxTwoQubitGateNumPerLayer_ += zone->front().nRows * zone->front().nCols;
  }
  if (maxTwoQubitGateNumPerLayer_ == 0) {
    throw std::invalid_argument("Architecture must contain at least one site "
                                "in an entanglement zone");
  }
}
auto ASAPScheduler::schedule(const qc::QuantumComputation& qc) const
    -> std::pair<
        std::vector<std::vector<std::reference_wrapper<const qc::Operation>>>,
        std::vector<std::vector<std::array<qc::Qubit, 2>>>> {
  if (qc.empty()) {
    // early exit if there are no operations to schedule
    return std::pair{
        std::vector<std::vector<std::reference_wrapper<const qc::Operation>>>{},
        std::vector<std::vector<std::array<qc::Qubit, 2>>>{}};
  }
  std::vector<std::vector<std::reference_wrapper<const qc::Operation>>>
      oneQubitGateLayers(1);
  std::vector<std::vector<std::array<qc::Qubit, 2>>> twoQubitGateLayers(0);
  // the following vector contains a mapping from qubits to the layer where
  // the next two-qubit gate can be scheduled for that qubit, i.e., the layer
  // after the last layer with a two-qubit gate acting on that qubit
  std::vector<size_t> nextLayerForQubit(qc.getNqubits(), 0);
  for (const auto& op : qc) {
    if (op->isGlobal(qc.getNqubits()) && qc.getNqubits() > 1) {
      const auto maxNextLayerForQubit = *std::max_element(
          nextLayerForQubit.cbegin(), nextLayerForQubit.cend());
      for (qc::Qubit q = 0; q < qc.getNqubits(); ++q) {
        nextLayerForQubit[q] = maxNextLayerForQubit;
      }
      oneQubitGateLayers[maxNextLayerForQubit].emplace_back(*op);
    } else if (op->isStandardOperation()) {
      const auto& stdOp = dynamic_cast<qc::StandardOperation&>(*op);
      if (stdOp.getNtargets() == 1 && stdOp.getNcontrols() == 0) {
        oneQubitGateLayers[nextLayerForQubit[stdOp.getTargets().front()]]
            .emplace_back(stdOp);
      } else if (stdOp.getType() == qc::Z && stdOp.getNtargets() == 1 &&
                 stdOp.getNcontrols() == 1) {
        const auto qubit1 = stdOp.getTargets().front();
        const auto qubit2 = stdOp.getControls().cbegin()->qubit;
        auto layer =
            std::max(nextLayerForQubit[qubit1], nextLayerForQubit[qubit2]);
        while (layer < twoQubitGateLayers.size() &&
               twoQubitGateLayers[layer].size() >=
                   maxTwoQubitGateNumPerLayer_) {
          ++layer;
        }
        assert(layer <= twoQubitGateLayers.size());
        if (layer == twoQubitGateLayers.size()) {
          // add a new layer
          oneQubitGateLayers.emplace_back();
          twoQubitGateLayers.emplace_back();
        }
        twoQubitGateLayers[layer].emplace_back(
            std::array{std::min(qubit1, qubit2), std::max(qubit1, qubit2)});
        nextLayerForQubit[qubit1] = layer + 1;
        nextLayerForQubit[qubit2] = layer + 1;
      } else {
        throw std::invalid_argument("Operation type not supported");
      }
    } else {
      throw std::invalid_argument("Operation type not supported");
    }
  }
  return std::pair{oneQubitGateLayers, twoQubitGateLayers};
}
} // namespace na
