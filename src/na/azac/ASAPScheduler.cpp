#include "na/azac/ASAPScheduler.hpp"

#include "Definitions.hpp"
#include "ir/QuantumComputation.hpp"
#include "ir/operations/OpType.hpp"
#include "ir/operations/Operation.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "na/azac/Architecture.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace na {
auto ASAPScheduler::schedule(const qc::QuantumComputation& qc) const
    -> std::pair<
        std::vector<std::vector<std::reference_wrapper<const qc::Operation>>>,
        std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>> {
  // calculate the maximum possible number of two-qubit gates per layer
  std::size_t maxTwoQubitGateNumPerLayer = 0;
  for (const auto& zone : architecture_.get().entanglementZones) {
    maxTwoQubitGateNumPerLayer += zone.front()->nRows * zone.front()->nCols;
  }
  std::vector<std::vector<std::reference_wrapper<const qc::Operation>>>
      oneQubitGateLayers(1);
  std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>> twoQubitGateLayers(
      0);
  // the following vector contains a mapping from qubits to the layer where
  // the next two-qubit gate can be scheduled for that qubit, i.e., the layer
  // after the last layer with a two-qubit gate acting on that qubit
  std::vector<size_t> nextLayerForQubit(qc.getNqubits(), 0);
  for (const auto& op : qc) {
    if (op->isStandardOperation()) {
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
               twoQubitGateLayers[layer].size() >= maxTwoQubitGateNumPerLayer) {
          ++layer;
        }
        assert(layer <= twoQubitGateLayers.size());
        if (layer == twoQubitGateLayers.size()) {
          // add a new layer
          oneQubitGateLayers.emplace_back();
          twoQubitGateLayers.emplace_back();
        }
        twoQubitGateLayers[layer].emplace_back(std::min(qubit1, qubit2),
                                               std::max(qubit1, qubit2));
        nextLayerForQubit[qubit1] = layer + 1;
        nextLayerForQubit[qubit2] = layer + 1;
      }
    } else {
      throw std::invalid_argument("Operation type not supported");
    }
  }
  return std::pair{oneQubitGateLayers, twoQubitGateLayers};
}
}; // namespace na
