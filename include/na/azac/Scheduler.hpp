#pragma once

#include "na/azac/Architecture.hpp"
#include "na/azac/CompilerBase.hpp"

#include <chrono>
#include <memory>
#include <stdexcept>

namespace na {

template <typename T> class Scheduler {
private:
  /// as soon as possible algorithm for g_q
  auto asap() -> std::vector<std::vector<std::size_t>> {
    std::vector<std::vector<std::size_t>> gate_scheduling{};
    std::vector<std::size_t> list_qubit_time(static_cast<T*>(this)->get_n_q(),
                                             0);
    for (std::size_t i = 0; i < static_cast<T*>(this)->get_n_g(); ++i) {
      const std::pair<qc::Qubit, qc::Qubit>& gate =
          static_cast<T*>(this)->get_g_q()[i];
      const auto tq0 = list_qubit_time[gate.first];
      const auto tq1 = list_qubit_time[gate.second];
      const auto tg = std::max(tq0, tq1);
      if (tg >= gate_scheduling.size()) {
        gate_scheduling.emplace_back();
      }
      gate_scheduling[tg].emplace_back(i);
      list_qubit_time[gate.first] = tg + 1;
      list_qubit_time[gate.second] = tg + 1;
    }
    return gate_scheduling;
  }

  auto graph_coloring() -> std::vector<std::vector<std::size_t>> {
    std::vector<std::vector<std::size_t>> gate_scheduling{};
    // todo: implement graph coloring algorithm
    throw std::runtime_error("Graph coloring algorithm is not implemented");
    return gate_scheduling;
  }

protected:
  /// solve a gate scheduling problem for all-commutable gate cases by graph
  /// coloring algorithm
  auto schedule() -> void {
    const auto t_s = std::chrono::system_clock::now();
    if (static_cast<T*>(this)->is_has_dependency()) {
      switch (static_cast<T*>(this)->get_scheduling_strategy()) {
      case CompilerBase::SchedulingStrategy::ASAP:
        static_cast<T*>(this)->set_gate_scheduling_idx(asap());
        break;
      case CompilerBase::SchedulingStrategy::TRIVIAL:
      default:
        std::vector<std::vector<std::size_t>> scheduling{};
        scheduling.reserve(static_cast<T*>(this)->get_n_g());
        for (std::size_t i = 0; i < static_cast<T*>(this)->get_n_g(); ++i) {
          scheduling.emplace_back(std::vector{i});
        }
        static_cast<T*>(this)->set_gate_scheduling_idx(scheduling);
      }
    } else {
      static_cast<T*>(this)->set_gate_scheduling_idx(graph_coloring());
    }
    // handle the case that the number of gates in one layer exceed the capacity
    // of rydberg zone
    std::size_t maxGateNum = 0;
    for (const std::vector<std::unique_ptr<SLM>>& zone :
         static_cast<T*>(this)->get_architecture().entanglement_zone) {
      maxGateNum += zone.front()->n_r * zone.front()->n_c;
    }
    // create a new scheduling where each group of gates has at most
    // max_gate_num
    std::vector<std::vector<std::size_t>> gateSchedulingSplit{};
    for (const std::vector<std::size_t>& gates :
         static_cast<T*>(this)->get_gate_scheduling_idx()) {
      if (gates.size() < maxGateNum) {
        gateSchedulingSplit.emplace_back(gates);
      } else {
        const std::size_t numGroup = (gates.size() - 1) / maxGateNum + 1;
        for (std::size_t i = 0; i < numGroup; ++i) {
          gateSchedulingSplit.emplace_back();
          for (std::size_t j = i * maxGateNum;
               j < std::min((i + 1) * maxGateNum, gates.size()); ++j) {
            gateSchedulingSplit.back().emplace_back(gates[j]);
          }
        }
      }
    }
    static_cast<T*>(this)->set_gate_scheduling_idx(
        std::move(gateSchedulingSplit));
    // clear old gate_scheduling and gate_1q_scheduling
    static_cast<T*>(this)->get_gate_scheduling().clear();
    static_cast<T*>(this)->get_gate_1q_scheduling().clear();
    // Re-construct valid gate_scheduling and gate_1q_scheduling after splitting
    // gate groups that exceed the capacity of rydberg zone
    for (const std::vector<std::size_t>& gates :
         static_cast<T*>(this)->get_gate_scheduling_idx()) {
      static_cast<T*>(this)->get_gate_scheduling().emplace_back();
      for (const std::size_t gateIdx : gates) {
        static_cast<T*>(this)->get_gate_scheduling().back().emplace_back(
            &static_cast<T*>(this)->get_g_q()[gateIdx]);
      }
      std::vector<const qc::StandardOperation*>& current_gate_1q =
          static_cast<T*>(this)->get_gate_1q_scheduling().emplace_back();
      for (const auto& gateIdx : gates) {
        const std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                                 std::vector<qc::StandardOperation>>& map =
            static_cast<T*>(this)->get_dict_g_1q_parent();
        const auto* gate = &static_cast<T*>(this)->get_g_q()[gateIdx];
        if (map.find(gate) != map.end()) {
          for (const auto& gate_1q : map.at(gate)) {
            current_gate_1q.emplace_back(&gate_1q);
          }
        }
      }
    }

    static_cast<T*>(this)->get_runtime_analysis().scheduling =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - t_s);

    std::cout
        << "[INFO]               Time for scheduling: "
        << static_cast<T*>(this)->get_runtime_analysis().scheduling.count()
        << "µs\n";
  }
};

} // namespace na
