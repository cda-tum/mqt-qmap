//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "cliffordsynthesis/SATEncoder.hpp"

#include "Encodings/Encodings.hpp"
#include "LogicUtil/util_logicblock.hpp"
#include "utils/logging.hpp"

#include <chrono>

namespace cs::encoding {

void SATEncoder::initializeSolver(const Configuration& config) {
  DEBUG() << "Initializing solver engine.";
  bool success        = false;
  LogicTerm::termType = TermType::BASE;
  logicutil::Params params;
  if (useMaxSAT) {
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

void SATEncoder::createTableauVariables() {
  DEBUG() << "Creating tableau variables.";
  for (std::size_t t = 0U; t <= T; ++t) {
    auto& x = vars.x.emplace_back();
    auto& z = vars.z.emplace_back();
    for (std::size_t i = 0U; i < N; ++i) {
      const std::string xName =
          "x_" + std::to_string(t) + "_" + std::to_string(i);
      TRACE() << "Creating variable " << xName;
      x.emplace_back(lb->makeVariable(xName, CType::BITVECTOR, N));
      const std::string zName =
          "z_" + std::to_string(t) + "_" + std::to_string(i);
      TRACE() << "Creating variable " << zName;
      z.emplace_back(lb->makeVariable(zName, CType::BITVECTOR, N));
    }
    const std::string rName = "r_" + std::to_string(t);
    TRACE() << "Creating variable " << rName;
    vars.r.emplace_back(lb->makeVariable(rName, CType::BITVECTOR, N));
  }
}

void SATEncoder::createSingleQubitGateVariables() {
  DEBUG() << "Creating single-qubit gate variables.";
  for (std::size_t t = 0U; t < T; ++t) {
    auto& timeStep = vars.gS.emplace_back();
    for (const auto gate : SINGLE_QUBIT_GATES) {
      auto& g = timeStep.emplace_back();
      for (std::size_t q = 0U; q < N; ++q) {
        const std::string gName = "g_" + std::to_string(t) + "_" +
                                  toString(gate) + "_" + std::to_string(q);
        TRACE() << "Creating variable " << gName;
        g.emplace_back(lb->makeVariable(gName));
      }
    }
  }
}

void SATEncoder::createTwoQubitGateVariables() {
  DEBUG() << "Creating two-qubit gate variables.";
  for (std::size_t t = 0U; t < T; ++t) {
    auto& timeStep = vars.gC.emplace_back();
    for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
      auto& control = timeStep.emplace_back();
      for (std::size_t trgt = 0U; trgt < N; ++trgt) {
        const std::string gName = "g_" + std::to_string(t) + "_cx_" +
                                  std::to_string(ctrl) + "_" +
                                  std::to_string(trgt);
        TRACE() << "Creating variable " << gName;
        control.emplace_back(lb->makeVariable(gName));
      }
    }
  }
}

void SATEncoder::asserTableau(const Tableau& tableau, const std::size_t pos) {
  DEBUG() << "Asserting tableau at time step " << pos;
  TRACE() << "Tableau:\n" << tableau;
  for (auto a = 0U; a < N; ++a) {
    const auto targetX = tableau.getBVFrom(a);
    lb->assertFormula(vars.x[pos][a] == LogicTerm(targetX, N));

    const auto targetZ = tableau.getBVFrom(a + N);
    lb->assertFormula(vars.z[pos][a] == LogicTerm(targetZ, N));
  }

  const auto targetR = tableau.getBVFrom(2U * N);
  lb->assertFormula(vars.r[pos] == LogicTerm(targetR, N));
}

void SATEncoder::createFormulation(const Tableau&       initialTableau,
                                   const Tableau&       targetTableau,
                                   const Configuration& config) {
  INFO() << "Creating formulation.";
  const auto start = std::chrono::high_resolution_clock::now();
  initializeSolver(config);

  createTableauVariables();
  createSingleQubitGateVariables();
  createTwoQubitGateVariables();

  asserTableau(initialTableau, 0U);
  asserTableau(targetTableau, T);

  if (!config.useMultiGateEncoding.has_value()) {
    // if not specified otherwise, the appropriate encoding is chosen based on
    // the target metric.
    if (config.target == TargetMetric::GATES) {
      createIndividualGateEncoding();
    } else {
      createMultiGateEncoding();
    }
  } else {
    if (*config.useMultiGateEncoding) {
      createMultiGateEncoding();
    } else {
      createIndividualGateEncoding();
    }
  }

  if (useMaxSAT || config.gateLimit.has_value()) {
    createObjectiveFunction(config);
  }

  const auto end = std::chrono::high_resolution_clock::now();
  const auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  INFO() << "Formulation created in " << duration << " ms.";
}

void SATEncoder::produceInstance() {
  INFO() << "Generating the SAT instance.";
  const auto start = std::chrono::high_resolution_clock::now();
  lb->produceInstance();
  const auto end = std::chrono::high_resolution_clock::now();
  const auto runtime =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  INFO() << "Instance generated in " << runtime << " ms.";

  DEBUG() << "Instance statistics:";
  DEBUG() << "  Clauses: " << TermImpl::getNextId(lb.get());
  DEBUG() << "  None terms: " << TermImpl::getNextId();
}

Result SATEncoder::solve() {
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

void SATEncoder::extractResultsFromModel(Results& res) {
  auto* const model = lb->getModel();
  extractCircuitFromModel(res, *model);
  extractTableauFromModel(res, T, *model);
}

void SATEncoder::extractCircuitFromModel(Results& res, Model& model) {
  std::size_t nSingleQubitGates = 0U;
  std::size_t nTwoQubitGates    = 0U;
  std::size_t depth             = 0U;

  qc::QuantumComputation qc(N);
  for (std::size_t t = 0; t < T; ++t) {
    const auto gateCount = nSingleQubitGates + nTwoQubitGates;
    DEBUG() << "Timestep: " << t << ", gate count: " << gateCount
            << ", depth: " << depth;

    const auto& singleQubitGates = vars.gS[t];
    for (std::size_t q = 0U; q < N; ++q) {
      for (const auto gate : SINGLE_QUBIT_GATES) {
        if (gate == qc::OpType::None) {
          continue;
        }
        if (model.getBoolValue(singleQubitGates[gateToIndex(gate)][q],
                               lb.get())) {
          qc.emplace_back<qc::StandardOperation>(N, q, gate);
          ++nSingleQubitGates;
          DEBUG() << toString(gate) << "(" << q << ")";
        }
      }
    }

    const auto& twoQubitGates = vars.gC[t];
    for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
      for (std::size_t trgt = 0U; trgt < N; ++trgt) {
        if (ctrl == trgt) {
          continue;
        }
        const auto control =
            dd::Control{static_cast<dd::Qubit>(ctrl), dd::Control::Type::pos};
        if (model.getBoolValue(twoQubitGates[ctrl][trgt], lb.get())) {
          qc.emplace_back<qc::StandardOperation>(N, control, trgt,
                                                 qc::OpType::X);
          ++nTwoQubitGates;
          DEBUG() << "CX(" << ctrl << ", " << trgt << ")";
        }
      }
    }
    if ((nSingleQubitGates + nTwoQubitGates) > gateCount) {
      ++depth;
    }
  }

  res.setSingleQubitGates(nSingleQubitGates);
  res.setTwoQubitGates(nTwoQubitGates);
  res.setDepth(depth);
  res.setResultCircuit(qc);
}

void SATEncoder::extractTableauFromModel(Results& res, const std::size_t pos,
                                         Model& model) {
  Tableau tableau(N);
  for (std::size_t i = 0; i < N; ++i) {
    const auto bvx = model.getBitvectorValue(vars.x[pos][i], lb.get());
    tableau.populateTableauFrom(bvx, N, i);
    const auto bvz = model.getBitvectorValue(vars.z[pos][i], lb.get());
    tableau.populateTableauFrom(bvz, N, i + N);
  }
  const auto bvr = model.getBitvectorValue(vars.r[pos], lb.get());
  tableau.populateTableauFrom(bvr, N, 2 * N);

  res.setResultTableau(tableau);
}

void SATEncoder::createObjectiveFunction(const Configuration& config) {
  DEBUG() << "Creating objective function.";

  switch (config.target) {
  case TargetMetric::GATES:
    createGateObjectiveFunction(config);
    break;
  case TargetMetric::DEPTH:
    createDepthObjectiveFunction();
    break;
  }
}

void SATEncoder::createGateObjectiveFunction(const Configuration& config) {
  auto cost = LogicTerm(0);
  for (std::size_t t = 0U; t < T; ++t) {
    const auto& singleQubitGates = vars.gS[t];
    for (std::size_t q = 0U; q < N; ++q) {
      for (const auto gate : SINGLE_QUBIT_GATES) {
        if (gate == qc::OpType::None) {
          continue;
        }
        cost = cost + singleQubitGates[gateToIndex(gate)][q];
      }
    }

    const auto& twoQubitGates = vars.gC[t];
    for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
      for (std::size_t trgt = 0U; trgt < ctrl; ++trgt) {
        cost = cost + (twoQubitGates[ctrl][trgt] + twoQubitGates[trgt][ctrl]);
      }
    }
  }

  if (!useMaxSAT && config.gateLimit.has_value()) {
    const auto gateLimit = config.gateLimit.value();
    DEBUG() << "Limiting number of gates to <= " << gateLimit;
    lb->assertFormula(cost <= LogicTerm(static_cast<int>(gateLimit)));
  } else {
    dynamic_cast<LogicBlockOptimizer*>(lb.get())->minimize(cost);
  }
}

void SATEncoder::createDepthObjectiveFunction() {
  auto cost = LogicTerm(0);
  for (std::size_t t = 0U; t < T; ++t) {
    auto anyGate = LogicTerm(false);

    const auto& singleQubitGates = vars.gS[t];
    for (std::size_t q = 0U; q < N; ++q) {
      for (const auto gate : SINGLE_QUBIT_GATES) {
        if (gate == qc::OpType::None) {
          continue;
        }
        anyGate = anyGate || singleQubitGates[gateToIndex(gate)][q];
      }
    }

    const auto& twoQubitGates = vars.gC[t];
    for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
      for (std::size_t trgt = 0U; trgt < ctrl; ++trgt) {
        anyGate =
            anyGate || (twoQubitGates[ctrl][trgt] || twoQubitGates[trgt][ctrl]);
      }
    }

    cost = cost + LogicTerm::ite(anyGate, LogicTerm(1), LogicTerm(0));
  }

  dynamic_cast<LogicBlockOptimizer*>(lb.get())->minimize(cost);
}

void SATEncoder::assertExactlyOne(const LogicVector& variables) {
  const auto variableGrouping =
      encodings::groupVars(variables, variables.size() / 2U);
  lb->assertFormula(encodings::exactlyOneCmdr(variableGrouping,
                                              LogicTerm::noneTerm(), lb.get()));
}

void SATEncoder::collectSingleQubitGateVariables(const std::size_t pos,
                                                 const std::size_t qubit,
                                                 LogicVector&      variables) {
  const auto& singleQubitGates = vars.gS[pos];
  for (const auto& gate : singleQubitGates) {
    variables.emplace_back(gate[qubit]);
  }
}

void SATEncoder::collectTwoQubitGateVariables(const std::size_t pos,
                                              const std::size_t qubit,
                                              const bool        target,
                                              LogicVector&      variables) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t q = 0; q < N; ++q) {
    if (q == qubit) {
      continue;
    }
    if (target) {
      variables.emplace_back(twoQubitGates[q][qubit]);
    } else {
      variables.emplace_back(twoQubitGates[qubit][q]);
    }
  }
}

LogicTerm SATEncoder::createSingleQubitXChange(const std::size_t pos,
                                               const std::size_t qubit,
                                               const qc::OpType  gate) {
  switch (gate) {
  case qc::OpType::None:
  case qc::OpType::X:
  case qc::OpType::Y:
  case qc::OpType::Z:
  case qc::OpType::S:
  case qc::OpType::Sdag:
    return vars.x[pos][qubit];
  case qc::OpType::H:
    return vars.z[pos][qubit];
  default:
    FATAL() << "Unsupported single-qubit gate: " << toString(gate);
    return LogicTerm::noneTerm();
  }
}

LogicTerm SATEncoder::createSingleQubitZChange(const std::size_t pos,
                                               const std::size_t qubit,
                                               const qc::OpType  gate) {
  switch (gate) {
  case qc::OpType::None:
  case qc::OpType::X:
  case qc::OpType::Y:
  case qc::OpType::Z:
    return vars.z[pos][qubit];
  case qc::OpType::H:
    return vars.x[pos][qubit];
  case qc::OpType::S:
  case qc::OpType::Sdag:
    return (vars.z[pos][qubit] ^ vars.x[pos][qubit]);
  default:
    FATAL() << "Unsupported single-qubit gate: " << toString(gate);
    return LogicTerm::noneTerm();
  }
}

LogicTerm SATEncoder::createSingleQubitRChange(const std::size_t pos,
                                               const std::size_t qubit,
                                               const qc::OpType  gate) {
  switch (gate) {
  case qc::OpType::None:
    return LogicTerm(0, N);
  case qc::OpType::H:
  case qc::OpType::S:
    return vars.x[pos][qubit] & vars.z[pos][qubit];
  case qc::OpType::Sdag:
    return vars.x[pos][qubit] & (vars.x[pos][qubit] ^ vars.z[pos][qubit]);
  case qc::OpType::X:
    return vars.z[pos][qubit];
  case qc::OpType::Y:
    return vars.x[pos][qubit] ^ vars.z[pos][qubit];
  case qc::OpType::Z:
    return vars.x[pos][qubit];
  default:
    FATAL() << "Unsupported single-qubit gate: " << toString(gate);
    return LogicTerm::noneTerm();
  }
}

std::pair<LogicTerm, LogicTerm>
SATEncoder::createTwoQubitXChange(const std::size_t pos, const std::size_t ctrl,
                                  const std::size_t trgt) {
  return {vars.x[pos][ctrl], vars.x[pos][ctrl] ^ vars.x[pos][trgt]};
}

std::pair<LogicTerm, LogicTerm>
SATEncoder::createTwoQubitZChange(const std::size_t pos, const std::size_t ctrl,
                                  const std::size_t trgt) {
  return {vars.z[pos][ctrl] ^ vars.z[pos][trgt], vars.z[pos][trgt]};
}

LogicTerm SATEncoder::createTwoQubitRChange(const std::size_t pos,
                                            const std::size_t ctrl,
                                            const std::size_t trgt) {
  const auto one = LogicTerm((1ULL << N) - 1, N);

  return (vars.x[pos][ctrl] & vars.z[pos][trgt]) &
         ((vars.z[pos][ctrl] ^ vars.x[pos][trgt]) ^ one);
}

LogicTerm SATEncoder::createNoChange(const std::size_t                pos,
                                     const std::size_t                except,
                                     const std::optional<std::size_t> except2) {
  auto changes = LogicTerm(true);
  for (std::size_t q = 0U; q < N; ++q) {
    if (q == except) {
      continue;
    }
    if (except2.has_value() && q == except2.value()) {
      continue;
    }

    changes = changes && (vars.x[pos + 1][q] == vars.x[pos][q]);
    changes = changes && (vars.z[pos + 1][q] == vars.z[pos][q]);
  }
  return changes;
}

void SATEncoder::assertIndividualGateConsistency() {
  DEBUG() << "Asserting individual gate consistency";
  for (std::size_t t = 0U; t < T; ++t) {
    LogicVector gateVariables{};

    for (std::size_t q = 0U; q < N; ++q) {
      collectSingleQubitGateVariables(t, q, gateVariables);
      collectTwoQubitGateVariables(t, q, true, gateVariables);
    }
    IF_PLOG(plog::verbose) {
      TRACE() << "Gate variables at time " << t;
      for (const auto& var : gateVariables) {
        TRACE() << var.getName();
      }
    }
    assertExactlyOne(gateVariables);
  }
}

void SATEncoder::assertIndividualGateSingleQubitGateConstraints(
    const std::size_t pos) {
  const auto& singleQubitGates = vars.gS[pos];
  for (std::size_t q = 0U; q < N; ++q) {
    for (const auto gate : SINGLE_QUBIT_GATES) {
      const auto changes =
          createIndividualGateSingleQubitGateConstraint(pos, q, gate);

      DEBUG() << "Asserting " << toString(gate) << " on " << q;
      IF_PLOG(plog::verbose) {
        std::stringstream ss;
        changes.prettyPrint(ss);
        TRACE() << "\n" << ss.str();
      }

      lb->assertFormula(
          LogicTerm::implies(singleQubitGates[gateToIndex(gate)][q], changes));
    }
  }
}

LogicTerm SATEncoder::createIndividualGateSingleQubitGateConstraint(
    const std::size_t pos, const std::size_t qubit, const qc::OpType gate) {
  auto changes = LogicTerm(true);

  changes = changes && (vars.x[pos + 1][qubit] ==
                        createSingleQubitXChange(pos, qubit, gate));
  changes = changes && (vars.z[pos + 1][qubit] ==
                        createSingleQubitZChange(pos, qubit, gate));
  changes =
      changes && (vars.r[pos + 1] ==
                  (vars.r[pos] ^ createSingleQubitRChange(pos, qubit, gate)));

  return changes && createNoChange(pos, qubit, std::nullopt);
}

void SATEncoder::assertIndividualGateTwoQubitGateConstraints(
    const std::size_t pos) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      const auto changes =
          createIndividualGateTwoQubitGateConstraint(pos, ctrl, trgt);

      DEBUG() << "Asserting CNOT on " << ctrl << " and " << trgt;
      IF_PLOG(plog::verbose) {
        std::stringstream ss;
        changes.prettyPrint(ss);
        TRACE() << "\n" << ss.str();
      }

      lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));
    }
  }
}

LogicTerm SATEncoder::createIndividualGateTwoQubitGateConstraint(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  auto changes              = LogicTerm(true);
  const auto [xCtrl, xTrgt] = createTwoQubitXChange(pos, ctrl, trgt);
  const auto [zCtrl, zTrgt] = createTwoQubitZChange(pos, ctrl, trgt);

  changes = changes && (vars.x[pos + 1][ctrl] == xCtrl);
  changes = changes && (vars.x[pos + 1][trgt] == xTrgt);
  changes = changes && (vars.z[pos + 1][ctrl] == zCtrl);
  changes = changes && (vars.z[pos + 1][trgt] == zTrgt);
  changes = changes && (vars.r[pos + 1] ==
                        (vars.r[pos] ^ createTwoQubitRChange(pos, ctrl, trgt)));

  return changes && createNoChange(pos, ctrl, trgt);
}

void SATEncoder::assertIndividualGateTableauConstraints() {
  DEBUG() << "Asserting individual gate tableau constraints";
  for (std::size_t t = 0U; t < T; ++t) {
    TRACE() << "Asserting tableau constraints at time " << t;
    assertIndividualGateSingleQubitGateConstraints(t);
    assertIndividualGateTwoQubitGateConstraints(t);
  }
}

void SATEncoder::createIndividualGateEncoding() {
  assertIndividualGateConsistency();
  assertIndividualGateTableauConstraints();
}

void SATEncoder::assertMultiGateConsistency() {
  DEBUG() << "Asserting multi-gate consistency";
  for (std::size_t t = 0U; t < T; ++t) {
    // asserting only a single gate is applied on each qubit.
    for (std::size_t q = 0U; q < N; ++q) {
      LogicVector gateVariables{};
      collectSingleQubitGateVariables(t, q, gateVariables);
      collectTwoQubitGateVariables(t, q, true, gateVariables);
      collectTwoQubitGateVariables(t, q, false, gateVariables);

      IF_PLOG(plog::verbose) {
        TRACE() << "Gate variables at time " << t << " and qubit " << q;
        for (const auto& var : gateVariables) {
          TRACE() << var.getName();
        }
      }

      assertExactlyOne(gateVariables);
    }
  }
}

void SATEncoder::assertMultiGateTableauConstraints() {
  DEBUG() << "Asserting multi-gate tableau constraints";
  for (std::size_t t = 0U; t < T; ++t) {
    TRACE() << "Asserting tableau constraints at time " << t;
    auto rChanges = vars.r[t];
    assertMultiGateSingleQubitGateConstraints(t, rChanges);
    assertMultiGateTwoQubitGateConstraints(t, rChanges);
    TRACE() << "Asserting r changes at time " << t;
    lb->assertFormula(vars.r[t + 1] == rChanges);
  }
}

void SATEncoder::assertMultiGateSingleQubitGateConstraints(
    const std::size_t pos, LogicTerm& rChanges) {
  const auto& singleQubitGates = vars.gS[pos];
  for (std::size_t q = 0U; q < N; ++q) {
    for (const auto gate : SINGLE_QUBIT_GATES) {
      const auto changes =
          createMultiGateSingleQubitGateConstraint(pos, q, gate);
      lb->assertFormula(
          LogicTerm::implies(singleQubitGates[gateToIndex(gate)][q], changes));

      rChanges =
          rChanges ^ LogicTerm::ite(singleQubitGates[gateToIndex(gate)][q],
                                    createSingleQubitRChange(pos, q, gate),
                                    LogicTerm(0, N));

      DEBUG() << "Asserting " << toString(gate) << " on " << q;
      IF_PLOG(plog::verbose) {
        std::stringstream ss;
        changes.prettyPrint(ss);
        TRACE() << "\n" << ss.str();
      }
    }
  }
}

LogicTerm SATEncoder::createMultiGateSingleQubitGateConstraint(
    const std::size_t pos, const std::size_t qubit, const qc::OpType gate) {
  auto changes = LogicTerm(true);

  changes = changes && (vars.x[pos + 1][qubit] ==
                        createSingleQubitXChange(pos, qubit, gate));
  changes = changes && (vars.z[pos + 1][qubit] ==
                        createSingleQubitZChange(pos, qubit, gate));

  return changes;
}

void SATEncoder::assertMultiGateTwoQubitGateConstraints(const std::size_t pos,
                                                        LogicTerm& rChanges) {
  const auto& twoQubitGates = vars.gC[pos];
  for (std::size_t ctrl = 0U; ctrl < N; ++ctrl) {
    for (std::size_t trgt = 0U; trgt < N; ++trgt) {
      if (ctrl == trgt) {
        continue;
      }
      const auto changes =
          createMultiGateTwoQubitGateConstraint(pos, ctrl, trgt);
      lb->assertFormula(LogicTerm::implies(twoQubitGates[ctrl][trgt], changes));

      rChanges =
          rChanges ^ LogicTerm::ite(twoQubitGates[ctrl][trgt],
                                    createTwoQubitRChange(pos, ctrl, trgt),
                                    LogicTerm(0, N));

      DEBUG() << "Asserting CNOT on " << ctrl << " and " << trgt;
      IF_PLOG(plog::verbose) {
        std::stringstream ss;
        changes.prettyPrint(ss);
        TRACE() << "\n" << ss.str();
      }
    }
  }
}

LogicTerm SATEncoder::createMultiGateTwoQubitGateConstraint(
    const std::size_t pos, const std::size_t ctrl, const std::size_t trgt) {
  auto changes              = LogicTerm(true);
  const auto [xCtrl, xTrgt] = createTwoQubitXChange(pos, ctrl, trgt);
  const auto [zCtrl, zTrgt] = createTwoQubitZChange(pos, ctrl, trgt);

  changes = changes && (vars.x[pos + 1][ctrl] == xCtrl);
  changes = changes && (vars.x[pos + 1][trgt] == xTrgt);
  changes = changes && (vars.z[pos + 1][ctrl] == zCtrl);
  changes = changes && (vars.z[pos + 1][trgt] == zTrgt);

  return changes;
}

void SATEncoder::createMultiGateEncoding() {
  assertMultiGateConsistency();
  assertMultiGateTableauConstraints();
}

void SATEncoder::cleanup() {
  if (lb) {
    lb->reset();
  }
}

} // namespace cs::encoding
