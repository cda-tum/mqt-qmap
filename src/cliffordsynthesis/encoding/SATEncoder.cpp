//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/SATEncoder.hpp"

#include "LogicUtil/util_logicblock.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "cliffordsynthesis/encoding/MultiGateEncoder.hpp"
#include "cliffordsynthesis/encoding/PhaseCorrectionEncoder.hpp"
#include "cliffordsynthesis/encoding/SingleGateEncoder.hpp"
#include "operations/OpType.hpp"
#include "utils/logging.hpp"

#include <chrono>
#include <cstddef>
#include <string>

namespace cs::encoding {

using namespace logicbase;

void SATEncoder::initializeSolver() {
  DEBUG() << "Initializing solver engine.";
  bool success        = false;
  LogicTerm::termType = TermType::BASE;
  logicutil::Params params;
  for (const auto& [key, value] : config.solverParameters) {
    if (std::holds_alternative<bool>(value)) {
      params.addParam(key, std::get<bool>(value));
    } else if (std::holds_alternative<std::uint32_t>(value)) {
      params.addParam(key, std::get<std::uint32_t>(value));
    } else if (std::holds_alternative<std::string>(value)) {
      params.addParam(key, std::get<std::string>(value));
    } else if (std::holds_alternative<double>(value)) {
      params.addParam(key, std::get<double>(value));
    } else {
      FATAL() << "Unknown parameter type.";
    }
  }

  if (config.useMaxSAT) {
    lb = logicutil::getZ3LogicOptimizer(success, true, params);
  } else {
    lb = logicutil::getZ3LogicBlock(success, true, params);
  }
  if (!success) {
    FATAL() << "Could not initialize solver engine.";
  }
}

void SATEncoder::createFormulation() {
  INFO() << "Creating formulation.";
  const auto start = std::chrono::high_resolution_clock::now();
  initializeSolver();

  S = config.targetTableau->hasDestabilizers() &&
              config.initialTableau->hasDestabilizers()
          ? 2U * N
          : N;

  tableauEncoder =
      std::make_shared<TableauEncoder>(N, S, T, lb, config.ignoreRChanges);
  tableauEncoder->createTableauVariables();
  tableauEncoder->assertTableau(*config.initialTableau, 0U);
  tableauEncoder->assertTableau(*config.targetTableau, T);

  if (config.ignoreRChanges) {
    config.gateSet.removePaulis();
    config.gateSet.emplace_back(qc::OpType::None);
  }

  if (config.useMultiGateEncoding) {
    gateEncoder = std::make_shared<MultiGateEncoder>(
        N, S, T, tableauEncoder->getVariables(), lb, config.gateSet,
        *config.initialTableau, config.ignoreRChanges);
  } else {
    gateEncoder = std::make_shared<SingleGateEncoder>(
        N, S, T, tableauEncoder->getVariables(), lb, config.gateSet,
        *config.initialTableau, config.ignoreRChanges);
  }
  gateEncoder->createSingleQubitGateVariables();
  gateEncoder->createTwoQubitGateVariables();
  gateEncoder->encodeGates();

  if (config.useSymmetryBreaking) {
    gateEncoder->encodeSymmetryBreakingConstraints();
  }

  objectiveEncoder = std::make_shared<ObjectiveEncoder>(
      N, T, gateEncoder->getVariables(), lb, config.gateSet);

  if (config.gateLimit.has_value()) {
    objectiveEncoder->limitGateCount(*config.gateLimit, std::less_equal{});
  }

  if (config.twoQubitGateLimit.has_value()) {
    objectiveEncoder->limitGateCount(*config.twoQubitGateLimit,
                                     std::less_equal{}, false);
  }

  if (config.useMaxSAT) {
    objectiveEncoder->optimizeMetric(config.targetMetric);
  }

  const auto end = std::chrono::high_resolution_clock::now();
  const auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  INFO() << "Formulation created in " << duration << " ms.";
}

Result SATEncoder::solve() const {
  INFO() << "Solving the SAT instance.";

  const auto start  = std::chrono::high_resolution_clock::now();
  const auto result = lb->solve();
  const auto end    = std::chrono::high_resolution_clock::now();
  const auto runtime =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  INFO() << "Instance solved in " << runtime << " ms.";
  return result;
}

void SATEncoder::extractResultsFromModel(Results& res) const {
  auto* const model = lb->getModel();
  tableauEncoder->extractTableauFromModel(res, T, *model);
  auto qc = gateEncoder->extractCircuitFromModel(res, *model);
  if (config.ignoreRChanges) {
    const auto start = std::chrono::high_resolution_clock::now();
    auto       tab   = *config.targetTableau;
    auto       qcTab = Tableau(qc, 0, std::numeric_limits<std::size_t>::max(),
                               config.targetTableau->hasDestabilizers());
    for (std::size_t row = 0U; row < S; ++row) {
      DEBUG() << "Row " << std::to_string(row);
      tab[row][2 * N] = qcTab.at(row)[2 * N];
    }
    PhaseCorrectionEncoder phaseCorrectionEncoder(N, S, tab,
                                                  *config.targetTableau);
    auto                   paulis = phaseCorrectionEncoder.phaseCorrection();

    for (std::size_t row = 0U; row < S; ++row) {
      tab[row][2 * N] = config.targetTableau->at(row)[2 * N];
    }

    for (std::size_t q = 0U; q < N; ++q) {
      if (paulis[q] != qc::OpType::None) {
        DEBUG() << "Phase correction for qubit " << q << ": "
                << qc::toString(paulis[q]);
        qc.emplace_back<qc::StandardOperation>(N, q, paulis[q]);
      }
    }
    const auto end     = std::chrono::high_resolution_clock::now();
    const auto runtime = std::chrono::duration<double>(end - start);
    res.setRuntime(res.getRuntime() + runtime.count());
    res.setResultCircuit(qc);
    res.setResultTableau(tab);
    res.setDepth(qc.getDepth());
  }
}

void SATEncoder::cleanup() const {
  if (lb) {
    lb->reset();
  }
}
Results SATEncoder::run() {
  const auto start = std::chrono::high_resolution_clock::now();

  createFormulation();
  const auto solverResult = solve();

  const auto end     = std::chrono::high_resolution_clock::now();
  const auto runtime = std::chrono::duration<double>(end - start);

  Results res{};
  res.setRuntime(runtime.count());
  res.setSolverResult(solverResult);

  if (solverResult == Result::SAT) {
    extractResultsFromModel(res);
  }

  cleanup();

  return res;
}

} // namespace cs::encoding
