//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"
#include "ir/QuantumComputation.hpp"

#include <cstdint>
#include <deque>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace na {
/**
 * @brief Struct to store the results of the scheduler
 */
struct SchedulerResults {
  qc::fp totalExecutionTime;
  qc::fp totalIdleTime;
  qc::fp totalGateFidelities;
  qc::fp totalFidelities;
  uint32_t nCZs = 0;

  SchedulerResults(const qc::fp executionTime, const qc::fp idleTime,
                   const qc::fp gateFidelities, const qc::fp fidelities,
                   const uint32_t cZs)
      : totalExecutionTime(executionTime), totalIdleTime(idleTime),
        totalGateFidelities(gateFidelities), totalFidelities(fidelities),
        nCZs(cZs) {}

  [[nodiscard]] std::string toString() const {
    std::stringstream ss;
    ss << "Total execution time: " << totalExecutionTime;
    ss << "\nTotal idle time: " << totalIdleTime
       << "\nTotal fidelities: " << totalFidelities;
    return ss.str();
  }
  [[nodiscard]] std::string toCsv() const {
    std::stringstream ss;
    ss << totalExecutionTime << ", " << totalIdleTime << "," << totalFidelities;
    return ss.str();
  }

  [[maybe_unused]] [[nodiscard]] std::unordered_map<std::string, qc::fp>
  toMap() const {
    std::unordered_map<std::string, qc::fp> result;
    result["totalExecutionTime"] = totalExecutionTime;
    result["totalIdleTime"] = totalIdleTime;
    result["totalGateFidelities"] = totalGateFidelities;
    result["totalFidelities"] = totalFidelities;
    result["nCZs"] = nCZs;
    return result;
  }
};

/**
 * @brief Class to schedule a quantum circuit on a neutral atom architecture
 * @details For each gate/operation in the input circuit, the scheduler checks
 * the earliest possible time slot for execution. If the gate is a multi qubit
 * gate, also the blocking of other qubits is taken into consideration. The
 * execution times are read from the neutral atom architecture.
 */
class NeutralAtomScheduler {
protected:
  const NeutralAtomArchitecture* arch = nullptr;
  std::string animationCsv;
  std::string animationArchitectureCsv;

public:
  // Constructor
  NeutralAtomScheduler() = default;
  explicit NeutralAtomScheduler(const NeutralAtomArchitecture& architecture)
      : arch(&architecture) {}

  /**
   * @brief Schedules the given quantum circuit on the neutral atom architecture
   * @details For each gate/operation in the input circuit, the scheduler checks
   * the earliest possible time slot for execution. If the gate is a multi qubit
   * gate, also the blocking of other qubits is taken into consideration. The
   * execution times are read from the neutral atom architecture.
   * @param qc Quantum circuit to schedule
   * @param verbose If true, prints additional information
   * @return SchedulerResults
   */
  SchedulerResults schedule(const qc::QuantumComputation& qc,
                            const std::map<HwQubit, HwQubit>& initHwPos,
                            bool verbose, bool createAnimationCsv = false,
                            qc::fp shuttlingSpeedFactor = 1.0);

  std::string getAnimationCsv() { return animationCsv; }
  void saveAnimationCsv(const std::string& filename) const {
    // save animation
    std::ofstream file(filename);
    file << animationCsv;
    file.close();
    // save architecture
    const auto filenameWithoutExtension =
        filename.substr(0, filename.find_last_of('.'));
    file.open(filenameWithoutExtension + "_architecture.csv");
    file << animationArchitectureCsv;
    file.close();
  }

  // Helper Print functions
  static void printSchedulerResults(std::vector<qc::fp>& totalExecutionTimes,
                                    qc::fp totalIdleTime,
                                    qc::fp totalGateFidelities,
                                    qc::fp totalFidelities, uint32_t nCZs);
  static void printTotalExecutionTimes(
      const std::vector<qc::fp>& totalExecutionTimes,
      const std::vector<std::deque<std::pair<qc::fp, qc::fp>>>&
          blockedQubitsTimes);
};

} // namespace na
