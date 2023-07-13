//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/SATEncoder.hpp"

#include "LogicUtil/util_logicblock.hpp"
#include "cliffordsynthesis/encoding/MultiGateEncoder.hpp"
#include "cliffordsynthesis/encoding/STQGatesEncoder.hpp"
#include "cliffordsynthesis/encoding/SingleGateEncoder.hpp"
#include "utils/logging.hpp"

#include <chrono>

namespace cs::encoding {

using namespace logicbase;
// TODO: params?
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

  const std::size_t s = config.targetTableau->hasDestabilizers() &&
                          config.initialTableau->hasDestabilizers()
                      ? 2U * N
                      : N;

  if (config.useSTEncoding) {
    T *= 2U;
  }

  tableauEncoder = std::make_shared<TableauEncoder>(N, s, T, lb);
  tableauEncoder->createTableauVariables();
  tableauEncoder->assertTableau(*config.initialTableau, 0U);
  tableauEncoder->assertTableau(*config.targetTableau, T);

  if (config.useMultiGateEncoding) {
    gateEncoder = std::make_shared<MultiGateEncoder>(
        N, s, T, T/2, tableauEncoder->getVariables(), lb);
  } else if (config.useSTEncoding) {
    gateEncoder = std::make_shared<STQGatesEncoder>(
        N, s, T, T/2, tableauEncoder->getVariables(), lb);
  } else {
    gateEncoder = std::make_shared<SingleGateEncoder>(
        N, s, T, T/2, tableauEncoder->getVariables(), lb);
  }
  gateEncoder->createSingleQubitGateVariables();
  gateEncoder->createTwoQubitGateVariables();
  if (config.useSTEncoding) {
    gateEncoder->addIdentityGateToTQGVariables();
  }

  gateEncoder->encodeGates();

  if (config.useSymmetryBreaking) {
    gateEncoder->encodeSymmetryBreakingConstraints();
  }

  objectiveEncoder =
      std::make_shared<ObjectiveEncoder>(N, T, gateEncoder->getVariables(), lb);

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
  gateEncoder->extractCircuitFromModel(res, *model);
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
