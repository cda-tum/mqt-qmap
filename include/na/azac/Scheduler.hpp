#pragma once

#include "na/azac/Architecture.hpp"
#include "na/azac/CompilerBase.hpp"

#include <type_traits>

namespace na {

template <typename T> class Scheduler {
private:
  /// as soon as possible algorithm for g_q
  auto asap() -> std::vector<
      std::unordered_set<const std::pair<qc::Qubit, qc::Qubit>*>> {
    std::vector<std::unordered_set<const std::pair<qc::Qubit, qc::Qubit>*>>
        gate_scheduling{};
    std::vector<std::size_t> list_qubit_time(static_cast<T*>(this)->n_q, 0);
    for (const std::pair<qc::Qubit, qc::Qubit>& gate :
         static_cast<T*>(this)->g_q) {
      const auto tq0 = list_qubit_time[gate.first];
      const auto tq1 = list_qubit_time[gate.second];
      const auto tg = std::max(tq0, tq1);
      if (tg >= gate_scheduling.size()) {
        gate_scheduling.emplace_back();
      }
      gate_scheduling[tg].insert(&gate);
      list_qubit_time[gate.first] = tg + 1;
      list_qubit_time[gate.second] = tg + 1;
    }
    return gate_scheduling;
  }
  auto graph_coloring() -> std::vector<
      std::unordered_set<const std::pair<qc::Qubit, qc::Qubit>*>> {
    std::vector<std::unordered_set<const std::pair<qc::Qubit, qc::Qubit>*>>
        gate_scheduling{};
    // todo: implement graph coloring algorithm
    throw std::runtime_error("Graph coloring algorithm is not implemented");
    return gate_scheduling;
  }

protected:
  /// solve gate scheduling problem for all-commutable gate cases by graph
  /// coloring algorithm
  auto schedule() -> void {
    const auto t_s = std::chrono::system_clock::now();
    if (static_cast<T*>(this)->hasDependency) {
      switch (static_cast<T*>(this)->scheduling_strategy) {
      case CompilerBase<T>::SchedulingStrategy::ASAP:
        static_cast<T*>(this)->gate_scheduling = asap();
        break;
      case CompilerBase<T>::SchedulingStrategy::TRIVIAL:
      default:
        static_cast<T*>(this)->gate_scheduling = std::accumulate(
            static_cast<T*>(this)->g_q.begin(),
            static_cast<T*>(this)->g_q.end(),
            std::vector<
                std::unordered_set<const std::pair<qc::Qubit, qc::Qubit>*>>{},
            [](std::vector<std::unordered_set<
                   const std::pair<qc::Qubit, qc::Qubit>*>>& acc,
               const std::pair<qc::Qubit, qc::Qubit>& gate) {
              acc.emplace_back(std::unordered_set{&gate});
              return acc;
            });
      }
    } else {
      static_cast<T*>(this)->gate_scheduling = graph_coloring();
    }
    // handle the case that the number of gates in one layer exceed the capacity
    // of rydberg zone
    std::size_t max_gate_num = 0;
    for (const std::vector<SLM>& zone :
         static_cast<T*>(this)->architecture.entanglement_zone) {
      max_gate_num += zone.front().n_r * zone.front().n_c;
    }
    // create a new scheduling where each group of gates has at most
    // max_gate_num
    std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>
        gate_scheduling_split{};
    for (auto gates : static_cast<T*>(this)->gate_scheduling) {
      if (gates.size() < max_gate_num) {
        gate_scheduling_split.emplace_back(gates);
      } else {
        while (!gates.empty()) {
          auto& split_gates = gate_scheduling_split.emplace_back(gates);
        }
      }
    }
    // clear old gate_scheduling and gate_1q_scheduling
    static_cast<T*>(this)->gate_scheduling.clear();
    static_cast<T*>(this)->gate_1q_scheduling.clear();
    // Re-construct valid gate_scheduling and gate_1q_scheduling after splitting
    // gate groups that exceed the capacity of rydberg zone
    for (const auto& gates : gate_scheduling_split) {
      static_cast<T*>(this)->gate_scheduling.emplace_back(gates);
      auto& current_gate_1q = static_cast<T*>(this)->gate_1q_scheduling.emplace_back();
      for (const auto& gate : gates) {
        if (static_cast<T*>(this)->dict_g_1q_parent.contains(gate)) {
          for (const auto& gate_1q :
               static_cast<T*>(this)->dict_g_1q_parent.at(gate)) {
            current_gate_1q.emplace(gate_1q);
          }
        }
      }
    }

    static_cast<T*>(this)->runtime_analysis.scheduling =
        std::chrono::system_clock::now() - t_s;

    std::cout << "[INFO]               Time for scheduling: " << std::fixed
              << std::setprecision(2)
              << static_cast<double>(
                     static_cast<T*>(this)->runtime_analysis.scheduling.count())
              << "µs\n";
  }
};

} // namespace na
