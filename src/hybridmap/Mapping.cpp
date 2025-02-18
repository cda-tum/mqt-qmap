//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "hybridmap/Mapping.hpp"

#include "hybridmap/NeutralAtomDefinitions.hpp"

#include <stdexcept>

namespace na {
void Mapping::applySwap(Swap swap) {
  auto q1 = swap.first;
  auto q2 = swap.second;
  if (this->isMapped(q1) && this->isMapped(q2)) {
    auto circQ1 = this->getCircQubit(q1);
    auto circQ2 = this->getCircQubit(q2);
    this->setCircuitQubit(circQ2, q1);
    this->setCircuitQubit(circQ1, q2);
  } else if (this->isMapped(q1) && !this->isMapped(q2)) {
    this->setCircuitQubit(this->getCircQubit(q1), q2);
  } else if (this->isMapped(q2) && !this->isMapped(q1)) {
    this->setCircuitQubit(this->getCircQubit(q2), q1);
  } else {
    throw std::runtime_error("Cannot swap unmapped qubits");
  }
}

std::vector<CoordIndex> Mapping::graphMatching() const {
  std::cout << "\n* initialMapping: Graph Matching\n";

  std::vector<CoordIndex> qubitIndices(
      dag.size(), std::numeric_limits<unsigned int>::max());
  std::vector<CoordIndex> hwIndices(hwQubits.getNumQubits(),
                                    std::numeric_limits<unsigned int>::max());

  // make hardware graph
  std::unordered_map<uint32_t, std::vector<uint32_t>> hwGraph;
  for (size_t i = 0; i < hwQubits.getNumQubits(); ++i) {
    auto neighbors = hwQubits.getNearbyQubits(i);
    hwGraph[i] = std::vector<uint32_t>(neighbors.begin(), neighbors.end());
  }
  for (auto& [qubit, neighbors] : hwGraph) {
    std::sort(neighbors.begin(), neighbors.end(),
              [this](uint32_t a, uint32_t b) {
                return hwQubits.getNearbyQubits(a).size() >
                       hwQubits.getNearbyQubits(b).size();
              });
  }

  uint32_t hwCenter = std::numeric_limits<unsigned int>::max();
  size_t maxHwConnections = 0;
  for (const auto& [qubit, neighbors] : hwGraph) {
    if (neighbors.size() > maxHwConnections) {
      maxHwConnections = neighbors.size();
      hwCenter = qubit;
    }
  }

  /*
  //for debug//
  for (size_t i = 0; i < hwQubits.getNumQubits(); ++i){
    auto coordIdx = hwQubits.getCoordIndex(i);
    std::cout << "hwQubit " << i << " -> coordIdx " << coordIdx << "\t";
    std::cout << "-> neighbors: ";
    for(auto j : hwGraph[i]){
      std::cout << j << " ";
    }
    std::cout << "\n";
  }
  std::cout << "=> hwCenter: " << hwCenter << "\n";
  */

  // make circuit graph
  std::vector<std::vector<std::pair<qc::Qubit, double>>> circGraph(dag.size());
  for (qc::Qubit qubit = 0; qubit < dag.size(); ++qubit) {
    std::unordered_map<qc::Qubit, double> weightMap;
    for (const auto& opPtr : dag[qubit]) {
      const auto* op = opPtr->get();
      auto usedQubits = op->getUsedQubits();
      if (usedQubits.size() > 1) {
        for (auto i : usedQubits) {
          if (i != qubit) {
            weightMap[i] += 1.0;
          }
        }
      }
    }
    std::vector<std::pair<qc::Qubit, double>> neighbors(weightMap.begin(),
                                                        weightMap.end());
    std::sort(
        neighbors.begin(), neighbors.end(),
        [](const std::pair<qc::Qubit, double>& a,
           std::pair<qc::Qubit, double>& b) { return a.second > b.second; });
    circGraph[qubit] = std::move(neighbors);
  }

  /*
  //for debug//
  for (qc::Qubit qubit = 0; qubit<circGraph.size(); ++qubit){
    std::cout << "Qubit " << qubit << " -> ";
    for(const auto& [neighbor, weight] : circGraph[qubit]){
      std::cout << "(q" << neighbor << " - weight " << weight << ") ";
    }
    std::cout << "\n";
  }
  */

  // circuit queue for graph matching
  std::vector<std::pair<int, std::pair<int, double>>> nodes;
  for (size_t i = 0; i < circGraph.size(); ++i) {
    int degree = circGraph[i].size();
    double weight_sum = 0;
    for (const auto& neighbor : circGraph[i]) {
      weight_sum += neighbor.second;
    }
    nodes.emplace_back(i, std::make_pair(degree, weight_sum));
  }
  std::sort(nodes.begin(), nodes.end(),
            [](const std::pair<int, std::pair<int, double>>& a,
               const std::pair<int, std::pair<int, double>>& b) {
              if (a.second.first == b.second.first) {
                return a.second.second > b.second.second;
              }
              return a.second.first > b.second.first;
            });
  std::queue<int> circGraphQueue;
  for (const auto& node : nodes) {
    circGraphQueue.push(node.first);
  }

  // graph matching -> return qubit Indices
  uint32_t nMapped = 0;
  bool firstCenter = true;
  int i = 0;
  while (!circGraphQueue.empty() && nMapped != dag.size()) {
    auto qi = circGraphQueue.front();
    HwQubit Qi = std::numeric_limits<unsigned int>::max();
    // std::cout << "*" << ++i << "th mapping: q" << qi << "\n";
    //  center mapping
    if (qubitIndices[qi] == std::numeric_limits<unsigned int>::max()) {
      // first center
      if (firstCenter) {
        Qi = hwCenter;
        firstCenter = false;
      }
      // next..
      else {
        int minDistance = std::numeric_limits<int>::max();
        for (HwQubit Qcandi = 0; Qcandi < hwQubits.getNumQubits(); ++Qcandi) {
          if (hwIndices[Qcandi] != std::numeric_limits<unsigned int>::max()) {
            continue;
          }
          int weightDistance = 0;
          for (auto qn_pair : circGraph[qi]) {
            auto qn = qn_pair.first;
            auto qn_weight = qn_pair.second;
            HwQubit Qn = qubitIndices[qn];
            if (Qn == std::numeric_limits<unsigned int>::max()) {
              continue;
            }
            weightDistance +=
                const_cast<na::HardwareQubits&>(hwQubits).getSwapDistance(
                    Qn, Qcandi, true) *
                qn_weight;
          }
          if (weightDistance < minDistance) {
            minDistance = weightDistance;
            Qi = Qcandi;
          }
        }
      }
      qubitIndices[qi] = Qi;
      hwIndices[Qi] = qi;
      nMapped++;
      // std::cout << "q" << qi << "-> Q" << Qi << "\n";
    } else {
      Qi = qubitIndices[qi];
    }
    // neighbor mapping
    for (auto& qn_pair : circGraph[qi]) {
      uint32_t qn = qn_pair.first;
      if (qubitIndices[qn] != std::numeric_limits<unsigned int>::max()) {
        continue;
      }
      HwQubit Qn = std::numeric_limits<unsigned int>::max();
      for (const auto& Qcandi : hwGraph[Qi]) {
        if (hwIndices[Qcandi] == std::numeric_limits<unsigned int>::max()) {
          Qn = Qcandi;
          break;
        }
      }
      if (Qn != std::numeric_limits<unsigned int>::max()) {
        qubitIndices[qn] = Qn;
        hwIndices[Qn] = qn;
        nMapped++;
        // std::cout << "q" << qn << "-> Q" << Qn << "\n";
      }
    }
    circGraphQueue.pop();
  }

  // for debug
  for (size_t i = 0; i < dag.size(); ++i) {
    if (qubitIndices[i] == std::numeric_limits<unsigned int>::max()) {
      for (HwQubit hw = 0; hw < hwQubits.getNumQubits(); ++hw) {
        if (hwIndices[hw] == std::numeric_limits<unsigned int>::max()) {
          qubitIndices[i] = hw;
          hwIndices[hw] = i;
          break;
        }
      }
    }
  }

  return qubitIndices;
}

} // namespace na
