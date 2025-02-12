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

  // init coord mapping
  std::vector<CoordIndex> qubitIndices(
      dag.size(), std::numeric_limits<unsigned int>::max());
  std::vector<CoordIndex> hwIndices(hwQubits.getNumQubits(),
                                    std::numeric_limits<unsigned int>::max());

  for (size_t i = 0; i < dag.size(); ++i) {
    qubitIndices[i] = i;
  }
  return qubitIndices;

  // use hwQubits instead of positions

  // // interaction graph
  // std::vector<std::vector<double>> circGraph(
  //     dag.size(), std::vector<double>(dag.size(), 0.0));
  // std::vector<std::pair<int, std::pair<int, double>>> circGraph_degree(
  //     dag.size());
  // std::vector<std::vector<std::pair<int, double>>> circGraph_neighbor(
  //     dag.size());
  // for (uint32_t qubit = 0; qubit < dag.size(); ++qubit) {
  //   for (const auto& opPtr : dag[qubit]) {
  //     const auto* op = opPtr->get();
  //     if (op->getUsedQubits().size() > 1) {
  //       for (auto i : op->getUsedQubits()) {
  //         if (i != qubit) {
  //           circGraph[qubit][i] += 1;
  //         }
  //       }
  //     }
  //   }
  // }
  //
  // // generate graph matching queue
  // for (uint32_t qubit = 0; qubit < dag.size(); qubit++) {
  //   int cnt = 0;
  //   double sum = 0;
  //   for (uint32_t i = 0; i < dag.size(); i++) {
  //     double weight = circGraph[qubit][i];
  //     if (weight > 0) {
  //       cnt++;
  //       sum += weight;
  //       circGraph_neighbor[qubit].emplace_back(i, weight);
  //     }
  //   }
  //   circGraph_degree[qubit] = std::make_pair(qubit, std::make_pair(cnt,
  //   sum));
  // }
  // sort(circGraph_degree.begin(), circGraph_degree.end(),
  //      [](std::pair<int, std::pair<int, double>> a,
  //         std::pair<int, std::pair<int, double>> b) {
  //        if (a.second.first == b.second.first)
  //          return a.second.second > b.second.second;
  //        else
  //          return a.second.first > b.second.first;
  //      });
  // std::queue<uint32_t> circGraph_queue;
  // for (const auto i : circGraph_degree) {
  //   circGraph_queue.push(i.first);
  // }
  // for (auto& innerVec : circGraph_neighbor) {
  //   sort(innerVec.begin(), innerVec.end(),
  //        [](std::pair<int, double>& a, std::pair<int, double>& b) {
  //          return a.second > b.second;
  //        });
  // }
  //
  // // graph matching
  // bool firstCenter = true;
  // uint32_t nMapped = 0;
  // uint32_t archCenterX =
  //     (archNcolumns % 2 == 0) ? (archNcolumns / 2 - 1) : (archNcolumns - 1) /
  //     2;
  // uint32_t archCenterY =
  //     (archNrows % 2 == 0) ? (archNrows / 2 - 1) : (archNrows - 1) / 2;
  // uint32_t archCenter = archCenterY * archNcolumns + archCenterX;
  // const auto start = hwQubits.getClosestQubit(archCenter, {});
  //
  // while (!circGraph_queue.empty() && nMapped != dag.size()) {
  //   uint32_t qc = circGraph_queue.front();
  //   uint32_t hc = std::numeric_limits<unsigned int>::max(); // hardwarCenter
  //   // center mapping
  //   if (firstCenter) {
  //     hc = archCenter;
  //     qubitIndices[qc] = hc;
  //     hwIndices[hc] = qc;
  //     firstCenter = false;
  //     nMapped++;
  //   } else if (qubitIndices[qc] == std::numeric_limits<unsigned int>::max())
  //   {
  //
  //     // ref loc
  //     std::vector<int> refLoc;
  //     for (auto i : circGraph_neighbor[qc]) {
  //       if (qubitIndices[i.first] != std::numeric_limits<unsigned
  //       int>::max()) {
  //         refLoc.push_back(qubitIndices[i.first]);
  //       }
  //     }
  //
  //     // candidate loc
  //     std::vector<std::pair<int, int>> distCandiLoc;
  //     std::vector<int> initCandiLoc;
  //     for (int i = 0; i < archNpositions; i++) {
  //       if (hwIndices[i] == std::numeric_limits<unsigned int>::max()) {
  //         initCandiLoc.push_back(i);
  //       }
  //     }
  //     for (auto v : initCandiLoc) {
  //       int distSum = 0;
  //       for (auto r : refLoc) {
  //         int dist = std::abs(v % archNcolumns - r % archNcolumns) +
  //                    std::abs(v / archNcolumns - r / archNcolumns);
  //         distSum += dist;
  //       }
  //       distCandiLoc.emplace_back(v, distSum);
  //     }
  //     sort(distCandiLoc.begin(), distCandiLoc.end(),
  //          [](std::pair<int, int> a, std::pair<int, int> b) {
  //            return a.second < b.second;
  //          });
  //
  //     // find position
  //     hc = distCandiLoc[0].first;
  //     qubitIndices[qc] = hc;
  //     hwIndices[hc] = qc;
  //     nMapped++;
  //   } else {
  //     hc = qubitIndices[qc];
  //   }
  //
  //   // neighbor mapping
  //   if (!circGraph_neighbor[qc].empty()) {
  //     int idx_qc_n = 0;
  //     std::vector<int> qc_n;
  //     for (auto i : circGraph_neighbor[qc]) {
  //       if (qubitIndices[i.first] != std::numeric_limits<unsigned
  //       int>::max())
  //         continue;
  //       if (idx_qc_n >= 4)
  //         continue;
  //       else {
  //         qc_n.push_back(i.first);
  //         idx_qc_n++;
  //       }
  //     }
  //     std::vector<int> hw_n;
  //     if (hc != std::numeric_limits<unsigned int>::max() && qc_n.size() > 0)
  //     {
  //       if ((hc + 1) % archNcolumns != 0 &&
  //           hwIndices[hc + 1] == std::numeric_limits<unsigned int>::max())
  //         hw_n.push_back(hc + 1); // right
  //       if (hc / archNcolumns < (archNrows - 1) &&
  //           hwIndices[hc + archNcolumns] ==
  //               std::numeric_limits<unsigned int>::max())
  //         hw_n.push_back(hc + archNcolumns); // down
  //       if (hc % archNcolumns != 0 &&
  //           hwIndices[hc - 1] == std::numeric_limits<unsigned int>::max())
  //         hw_n.push_back(hc - 1); // left
  //       if (hc > archNcolumns && hwIndices[hc - archNcolumns] ==
  //                                    std::numeric_limits<unsigned
  //                                    int>::max())
  //         hw_n.push_back(hc - archNcolumns); // up
  //     }
  //
  //     int minSize = std::min(qc_n.size(), hw_n.size());
  //     for (int i = 0; i < minSize; i++) {
  //       int qc_i = qc_n[i];
  //       int hw_i = hw_n[i];
  //       qubitIndices[qc_i] = hw_i;
  //       hwIndices[hw_i] = qc_i;
  //       nMapped++;
  //     }
  //   }
  //   circGraph_queue.pop();
  // }
  // return qubitIndices;
}

} // namespace na
