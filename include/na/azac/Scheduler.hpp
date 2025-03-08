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
    std::vector<std::vector<std::size_t>> gateScheduling{};
    std::vector<std::size_t> listQubitTime(static_cast<T*>(this)->getNQubits(),
                                           0);
    for (std::size_t i = 0; i < static_cast<T*>(this)->getNTwoQubitGates();
         ++i) {
      const auto& gate = static_cast<T*>(this)->getTwoQubitGates()[i];
      const auto tq0 = listQubitTime[gate->first];
      const auto tq1 = listQubitTime[gate->second];
      const auto tg = std::max(tq0, tq1);
      if (tg >= gateScheduling.size()) {
        gateScheduling.emplace_back();
      }
      gateScheduling[tg].emplace_back(i);
      listQubitTime[gate->first] = tg + 1;
      listQubitTime[gate->second] = tg + 1;
    }
    return gateScheduling;
  }

  auto graphColoring() -> std::vector<std::vector<std::size_t>> {
    std::vector<std::vector<std::size_t>> gateScheduling{};
    // todo: implement graph coloring algorithm
    throw std::runtime_error("Graph coloring algorithm is not implemented");
    return gateScheduling;
  }

protected:
  /// solve a gate scheduling problem for all-commutable gate cases by graph
  /// coloring algorithm
  auto schedule() -> void {
    const auto t_s = std::chrono::system_clock::now();
    if (static_cast<T*>(this)->isHasDependency()) {
      switch (static_cast<T*>(this)->getSchedulingStrategy()) {
      case CompilerBase::SchedulingStrategy::ASAP:
        static_cast<T*>(this)->setGateSchedulingIdx(asap());
        break;
      case CompilerBase::SchedulingStrategy::TRIVIAL:
      default:
        std::vector<std::vector<std::size_t>> schedule{};
        schedule.reserve(static_cast<T*>(this)->getNTwoQubitGates());
        for (std::size_t i = 0; i < static_cast<T*>(this)->getNTwoQubitGates();
             ++i) {
          schedule.emplace_back(std::vector{i});
        }
        static_cast<T*>(this)->setGateSchedulingIdx(schedule);
      }
    } else {
      static_cast<T*>(this)->setGateSchedulingIdx(graphColoring());
    }
    // handle the case that the number of gates in one layer exceed the capacity
    // of rydberg zone
    std::size_t maxGateNum = 0;
    for (const std::vector<std::unique_ptr<SLM>>& zone :
         static_cast<T*>(this)->getArchitecture().entanglementZones) {
      maxGateNum += zone.front()->nRows * zone.front()->nCols;
    }
    // create a new scheduling where each group of gates has at most
    // maxGatenum
    std::vector<std::vector<std::size_t>> gateSchedulingSplit{};
    for (const std::vector<std::size_t>& gates :
         static_cast<T*>(this)->getGateSchedulingIdx()) {
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
    static_cast<T*>(this)->setGateSchedulingIdx(std::move(gateSchedulingSplit));
    // clear old gateScheduling and gate1qscheduling
    static_cast<T*>(this)->getGateScheduling().clear();
    static_cast<T*>(this)->getGate1QScheduling().clear();
    // Re-construct valid gateScheduling and gate1qscheduling after splitting
    // gate groups that exceed the capacity of rydberg zone
    for (const std::vector<std::size_t>& gates :
         static_cast<T*>(this)->getGateSchedulingIdx()) {
      static_cast<T*>(this)->getGateScheduling().emplace_back();
      for (const std::size_t gateIdx : gates) {
        static_cast<T*>(this)->getGateScheduling().back().emplace_back(
            static_cast<T*>(this)->getTwoQubitGates()[gateIdx].get());
      }
      std::vector<const qc::StandardOperation*>& currentGate1q =
          static_cast<T*>(this)->getGate1QScheduling().emplace_back();
      for (const auto& gateIdx : gates) {
        const std::unordered_map<const std::pair<qc::Qubit, qc::Qubit>*,
                                 std::vector<qc::StandardOperation>>& map =
            static_cast<T*>(this)->getDictG1QParent();
        const auto& gate = static_cast<T*>(this)->getTwoQubitGates()[gateIdx];
        if (map.find(gate.get()) != map.end()) {
          for (const auto& gate1Q : map.at(gate.get())) {
            currentGate1q.emplace_back(&gate1Q);
          }
        }
      }
    }

    static_cast<T*>(this)->getRuntimeAnalysis().scheduling =
        std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::system_clock::now() - t_s);

    std::cout << "[INFO]           Time for scheduling: "
              << static_cast<T*>(this)->getRuntimeAnalysis().scheduling.count()
              << "µs\n";
  }
};

} // namespace na
