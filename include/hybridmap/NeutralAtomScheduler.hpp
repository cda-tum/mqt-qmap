//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "Definitions.hpp"
#include "QuantumComputation.hpp"
#include "hybridmap/NeutralAtomArchitecture.hpp"
#include "hybridmap/NeutralAtomDefinitions.hpp"

#include <cstdint>
#include <deque>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace qc {
/**
 * @brief Struct to store the results of the scheduler
 */
struct SchedulerResults {
  fp       totalExecutionTime;
  fp       totalIdleTime;
  fp       totalGateFidelities;
  fp       totalFidelities;
  uint32_t nCZs = 0;

  SchedulerResults(fp totalExecutionTime, fp totalIdleTime,
                   fp totalGateFidelities, fp totalFidelities, uint32_t nCZs)
      : totalExecutionTime(totalExecutionTime), totalIdleTime(totalIdleTime),
        totalGateFidelities(totalGateFidelities),
        totalFidelities(totalFidelities), nCZs(nCZs) {}

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

  [[maybe_unused]] [[nodiscard]] std::unordered_map<std::string, fp>
  toMap() const {
    std::unordered_map<std::string, fp> result;
    result["totalExecutionTime"]  = totalExecutionTime;
    result["totalIdleTime"]       = totalIdleTime;
    result["totalGateFidelities"] = totalGateFidelities;
    result["totalFidelities"]     = totalFidelities;
    result["nCZs"]                = nCZs;
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
  qc::NeutralAtomArchitecture arch;
  std::string                 animationCsv;
  std::string                 animationArchitectureCsv;

public:
  // Constructor
  explicit NeutralAtomScheduler(qc::NeutralAtomArchitecture arch)
      : arch(std::move(arch)) {}

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
  SchedulerResults schedule(const qc::QuantumComputation&     qc,
                            const std::map<HwQubit, HwQubit>& initHwPos,
                            bool verbose, bool createAnimationCsv = false,
                            fp shuttlingSpeedFactor = 1.0);

  std::string getAnimationCsv() { return animationCsv; }
  void        saveAnimationCsv(const std::string& filename) {
    // save animation
    std::ofstream file(filename);
    file << animationCsv;
    file.close();
    // save architecture
    auto filenameWithoutExtension =
        filename.substr(0, filename.find_last_of('.'));
    file.open(filenameWithoutExtension + "_architecture.csv");
    file << animationArchitectureCsv;
    file.close();
  }

  // Helper Print functions
  static void printSchedulerResults(std::vector<fp>& totalExecutionTimes,
                                    fp totalIdleTime, fp totalGateFidelities,
                                    fp totalFidelities, uint32_t nCZs);
  static void printTotalExecutionTimes(
      std::vector<fp>&                            totalExecutionTimes,
      std::vector<std::deque<std::pair<fp, fp>>>& blockedQubitsTimes);
};

class AnimationAtoms {
  using axesId   = std::uint32_t;
  using marginId = std::uint32_t;

protected:
  uint32_t                  colorSlm    = 0;
  uint32_t                  colorAod    = 1;
  uint32_t                  colorLocal  = 2;
  [[maybe_unused]] uint32_t colorGlobal = 3;
  uint32_t                  colorCz     = 4;

  std::map<CoordIndex, HwQubit>        coordIdxToId;
  std::map<HwQubit, std::pair<fp, fp>> idToCoord;
  std::map<HwQubit, uint32_t>          axesIds;
  std::map<HwQubit, uint32_t>          marginIds;
  uint32_t                             axesIdCounter   = 0;
  uint32_t                             marginIdCounter = 0;

  axesId   addAxis(HwQubit id);
  void     removeAxis(HwQubit id) { axesIds.erase(id); }
  marginId addMargin(HwQubit id);
  void     removeMargin(HwQubit id) { marginIds.erase(id); }

public:
  AnimationAtoms(const std::map<HwQubit, HwQubit>& initHwPos,
                 const NeutralAtomArchitecture&    arch);

  std::string        getInitString();
  std::string        getEndString(fp endTime);
  static std::string createCsvLine(fp startTime, HwQubit id, fp x, fp y,
                                   uint32_t size = 1, uint32_t color = 0,
                                   bool axes = false, axesId axId = 0,
                                   bool margin = false, marginId marginId = 0,
                                   fp marginSize = 0);
  std::string createCsvOp(const std::unique_ptr<Operation>& op, fp startTime,
                          fp endTime, const qc::NeutralAtomArchitecture& arch);
};

} // namespace qc