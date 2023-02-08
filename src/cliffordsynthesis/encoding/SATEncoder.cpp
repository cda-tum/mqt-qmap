//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/encoding/SATEncoder.hpp"

#include "LogicUtil/util_logicblock.hpp"
#include "cliffordsynthesis/encoding/MultiGateEncoder.hpp"
#include "cliffordsynthesis/encoding/SingleGateEncoder.hpp"
#include "utils/logging.hpp"

#include <chrono>

namespace cs::encoding {

using namespace logicbase;

void SATEncoder::initializeSolver() {
  DEBUG() << "Initializing solver engine.";
  bool success        = false;
  LogicTerm::termType = TermType::BASE;
  logicutil::Params params;
  if (config.useMaxSAT) {
    params.addParam("pb.compile_equality", true);
    params.addParam("maxres.hill_climb", true);
    params.addParam("maxres.pivot_on_correction_set", false);
    lb = logicutil::getZ3LogicOptimizer(success, true, params);
  } else {
    params.addParam("threads", static_cast<std::uint32_t>(config.nThreads / 2));
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

  tableauEncoder = std::make_shared<TableauEncoder>(N, s, T, lb);
  tableauEncoder->createTableauVariables();
  tableauEncoder->assertTableau(*config.initialTableau, 0U);
  tableauEncoder->assertTableau(*config.targetTableau, T);

  if (config.useMultiGateEncoding) {
    gateEncoder = std::make_shared<MultiGateEncoder>(
        N, s, T, tableauEncoder->getVariables(), lb);
  } else {
    gateEncoder = std::make_shared<SingleGateEncoder>(
        N, s, T, tableauEncoder->getVariables(), lb);
  }
  gateEncoder->createSingleQubitGateVariables();
  gateEncoder->createTwoQubitGateVariables();
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

void SATEncoder::produceInstance() const {
  INFO() << "Generating the SAT instance.";
  const auto start = std::chrono::high_resolution_clock::now();
  lb->produceInstance();
  const auto end = std::chrono::high_resolution_clock::now();
  const auto runtime =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  INFO() << "Instance generated in " << runtime << " ms.";

  DEBUG() << "Instance statistics:";
  DEBUG() << "\tClauses: " << TermImpl::getNextId(lb.get());
  DEBUG() << "\tNone terms: " << TermImpl::getNextId();
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
  produceInstance();
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
