/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "LogicBlock/LogicBlock.hpp"
#include "LogicTerm/LogicTerm.hpp"
#include "LogicUtil/util_logicblock.hpp"
#include "cliffordsynthesis/ExactStrategy.hpp"
#include "cliffordsynthesis/GateEncoding.hpp"
#include "cliffordsynthesis/HeuristicStrategy.hpp"
#include "cliffordsynthesis/OptimizationStrategy.hpp"
#include "cliffordsynthesis/TargetMetricHandler.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"
#include "utils.hpp"
#include "utils/logging.hpp"

namespace cs {

void CliffordSynthesizer::initConfiguration(Configuration& configuration) {
  // we dont already have a tableau
  if (configuration.initialTableau == nullptr ||
      configuration.targetTableau == nullptr ||
      configuration.initialTableau->empty() ||
      configuration.targetTableau->empty()) {
    configuration.nqubits = configuration.targetCircuit->getNqubits();
    configuration.initialTableau =
        std::make_shared<Tableau>(configuration.nqubits);
    configuration.targetTableau =
        std::make_shared<Tableau>(*configuration.targetCircuit);
  }
}

void CliffordSynthesizer::synthesize(Configuration& configuration) {
  TRACE() << "OptimizationStrategy: " << toString(configuration.strategy)
          << std::endl;
  TRACE() << "Target: " << toString(configuration.target) << std::endl;

  initConfiguration(configuration);

  initResults();

  auto totalStart = std::chrono::high_resolution_clock::now();
  initCouplingMaps(configuration);

  for (const auto& subset : couplingMaps) {
    auto qubitSet{Architecture::getQubitSet(subset)};

    DEBUG() << "Reduced Coupling Map"
            << (configuration.chooseBest ? " (best)" : "") << ": ";
    std::stringstream strings;
    Architecture::printCouplingMap(subset, strings);
    DEBUG() << strings.str();
    DEBUG() << "Qubit Set: " << qubitSet;
    DEBUG() << "Coupling Map Fidelity: "
            << Architecture::getAverageArchitectureFidelity(
                   configuration.architecture.getCouplingMap(), qubitSet,
                   configuration.architecture.getProperties());
    std::size_t timesteps = configuration.initialTimestep;
    if (timesteps == 0U) {
      timesteps = configuration.nqubits * configuration.nqubits;
    }
    if (isExact(configuration.strategy)) {
      ExactStrategy::runExactStrategy(timesteps, subset, qubitSet,
                                      configuration, *this);
    } else {
      HeuristicStrategy::runHeuristicStrategy(subset, qubitSet, configuration,
                                              *this);
    }

    if (configuration.chooseBest && optimalResults.sat) {
      break;
    }
  }

  auto totalEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = totalEnd - totalStart;
  INFO() << "Time: " << diff.count() << std::endl;
  optimalResults.totalSeconds = diff.count();
}

Results CliffordSynthesizer::mainOptimization(
    const std::size_t timesteps, const CouplingMap& reducedCM,
    const QubitSubset& qubitChoice, const Tableau& targetTableau,
    const Tableau& initialTableau, const Configuration& configuration) {
  using namespace logicbase;
  std::unique_ptr<LogicBlock> lb;
  bool                        success = false;
  LogicTerm::termType                 = TermType::BASE;
  if (configuration.strategy == OptimizationStrategy::UseMinimizer ||
      configuration.strategy == OptimizationStrategy::SplitIter) {
    logicutil::Params params;
    params.addParam("pb.compile_equality", true);
    params.addParam("maxres.hill_climb", true);
    params.addParam("maxres.pivot_on_correction_set", false);
    lb = logicutil::getZ3LogicOptimizer(success, true, params);
  } else {
    logicutil::Params params;
    params.addParam("threads",
                    static_cast<unsigned>(configuration.nThreads / 2));
    lb = logicutil::getZ3LogicBlock(success, true, params);
  }
  if (!success) {
    throw QMAPException("Could not initialize Z3 logic block optimizer");
  }
  DEBUG() << "lb 1: " << lb.get() << std::endl;
  logicbase::LogicMatrix   x{};
  logicbase::LogicMatrix   z{};
  logicbase::LogicVector   r{};
  logicbase::LogicMatrix3D gS{};
  logicbase::LogicMatrix3D gC{};

  auto changes = logicbase::LogicTerm(true);

  auto start = std::chrono::high_resolution_clock::now();
  /*
   * Tableau Variables x/z
   * k before gate k is applied
   * i column
   */
  std::stringstream xName{};
  std::stringstream zName{};
  std::stringstream rName{};
  for (std::size_t k = 0U; k < timesteps + 1U; ++k) {
    x.emplace_back();
    z.emplace_back();
    for (std::size_t i = 0U; i < configuration.nqubits; ++i) {
      xName.str("");
      zName.str("");
      xName << "x_" << k << "_" << i;
      zName << "z_" << k << "_" << i;
      x.back().emplace_back(lb->makeVariable(xName.str(), CType::BITVECTOR,
                                             configuration.nqubits));
      z.back().emplace_back(lb->makeVariable(zName.str(), CType::BITVECTOR,
                                             configuration.nqubits));
    }
    rName.str("");
    rName << "r_" << k;
    r.emplace_back(
        lb->makeVariable(rName.str(), CType::BITVECTOR, configuration.nqubits));
  }

  /*
   * Gate Variables
   * k before gates k are applied
   * i qubit 1
   * j qubit
   */
  std::stringstream gName{};
  for (std::size_t gateStep = 0U; gateStep < timesteps + 1U; ++gateStep) {
    gS.emplace_back();
    for (const auto gate : Gates::SINGLE_QUBIT) {
      gS.back().emplace_back();
      for (std::size_t j = 0U; j < configuration.nqubits; ++j) {
        gName.str("");
        gName << "g_" << gateStep << "_" << Gates::gateName(gate) << "_" << j;
        gS.back().back().emplace_back(lb->makeVariable(gName.str()));
      }
    }
  }
  for (std::size_t gateStep = 0U; gateStep < timesteps + 1U; ++gateStep) {
    gC.emplace_back();
    for (std::size_t j = 0U; j < configuration.nqubits; ++j) {
      gC.back().emplace_back();
      for (std::size_t l = 0U; l < configuration.nqubits; ++l) {
        gName.str("");
        gName << "g_" << gateStep << "_CNOT_" << j << "_" << l;
        gC.back().back().emplace_back(lb->makeVariable(gName.str()));
      }
    }
  }

  assertTableau(SynthesisData{configuration.nqubits, timesteps, reducedCM,
                              qubitChoice, lb, x, z, r, gS, gC},
                initialTableau, 0U);
  assertTableau(SynthesisData{configuration.nqubits, timesteps, reducedCM,
                              qubitChoice, lb, x, z, r, gS, gC},
                targetTableau, timesteps);

  // assert gate limits
  GateEncoding::makeGateEncoding(SynthesisData{configuration.nqubits, timesteps,
                                               reducedCM, qubitChoice, lb, x, z,
                                               r, gS, gC},
                                 configuration);

  // assert cost functions for respective target metrics
  TargetMetricHandler::makeTargetMetric(
      SynthesisData{configuration.nqubits, timesteps, reducedCM, qubitChoice,
                    lb, x, z, r, gS, gC},
      configuration);

  auto formulation = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = formulation - start;
  INFO() << "Time to produce Formulation: " << diff.count() << std::endl;

  lb->produceInstance();

  auto modGen = std::chrono::high_resolution_clock::now();
  diff        = modGen - formulation;
  INFO() << "Time to generate Model: " << diff.count() << std::endl;

  TRACE() << "Clauses: " << TermImpl::getNextId(lb.get()) << std::endl;
  TRACE() << "None Terms: " << TermImpl::getNextId() << std::endl;

  Result result = lb->solve();

  auto end = std::chrono::high_resolution_clock::now();
  diff     = end - modGen;
  INFO() << "Time to solve Model: " << diff.count() << std::endl;
  std::chrono::duration<double, std::milli> elapsedMilliseconds = end - start;
  Results                                   results{};
  results.verbose          = configuration.verbosity;
  results.chooseBest       = configuration.chooseBest;
  results.nqubits          = configuration.nqubits;
  results.initialTimesteps = timesteps;
  results.strategy         = configuration.strategy;
  results.target           = configuration.target;
  results.totalSeconds     = elapsedMilliseconds.count();
  results.sat              = result == Result::SAT;
  results.doubleFidelity   = configuration.architecture.getFidelityTable();
  results.singleFidelity =
      configuration.architecture.getSingleQubitFidelities();
  results.resultCM = configuration.architecture.getCouplingMap();
  results.resultTableaus.clear();

  if (result == Result::SAT) {
    results.result               = logicbase::Result::SAT;
    auto*                  model = lb->getModel();
    qc::QuantumComputation localResultCircuit{};
    if (configuration.architecture.isArchitectureAvailable()) {
      localResultCircuit.addQubitRegister(
          configuration.architecture.getNqubits());
    } else {
      localResultCircuit.addQubitRegister(configuration.nqubits);
    }
    results.singleQubitGates = 0U;
    results.twoQubitGates    = 0U;
    results.depth            = 0U;
    results.fidelity         = 1.;
    assert(configuration.nqubits == qubitChoice.size());
    for (std::size_t gateStep = 0U; gateStep < timesteps + 1U; ++gateStep) {
      auto oldGateCount = results.singleQubitGates + results.twoQubitGates;
      TRACE() << "Gate Step: " << gateStep << std::endl
              << " Actual gate count: " << oldGateCount << std::endl
              << " Depth: " << results.depth << std::endl
              << " Fidelity: " << results.fidelity << std::endl;
      if (gateStep > 0U) {
        auto aIt = qubitChoice.cbegin();
        for (std::size_t a = 0U; a < configuration.nqubits; ++a) {
          const auto q0 = *aIt;
          for (auto gate : Gates::SINGLE_QUBIT_WITHOUT_NOP) {
            if (model->getBoolValue(gS[gateStep][Gates::toIndex(gate)][a],
                                    lb.get())) {
              localResultCircuit.emplace_back<qc::StandardOperation>(
                  configuration.nqubits, q0, Gates::toOpType(gate));
              if (configuration.architecture.isCalibrationDataAvailable()) {
                results.fidelity *=
                    (configuration.architecture.getSingleQubitFidelities()[q0]);
              }
              TRACE() << Gates::gateName(gate) << "(" << q0 << ")" << std::endl;
              if (configuration.architecture.isCalibrationDataAvailable()) {
                TRACE()
                    << " Fidelity: "
                    << configuration.architecture.getSingleQubitFidelities()[q0]
                    << std::endl;
              }
              ++results.singleQubitGates;
            }
          }
          auto bIt = qubitChoice.cbegin();
          for (std::size_t b = 0U; b < configuration.nqubits; ++b) {
            const auto q1 = *bIt;
            if (model->getBoolValue(gC[gateStep][a][b], lb.get())) {
              localResultCircuit.emplace_back<qc::StandardOperation>(
                  configuration.nqubits,
                  dd::Control{static_cast<dd::Qubit>(q0)}, q1, qc::X);
              if (configuration.architecture.isCalibrationDataAvailable()) {
                results.fidelity *=
                    configuration.architecture.getFidelityTable()[q0][q1];
              }
              TRACE() << "X(" << q0 << "," << q1 << ")" << std::endl;
              if (configuration.architecture.isCalibrationDataAvailable()) {
                TRACE() << "Fidelity: "
                        << configuration.architecture.getFidelityTable()[q0][q1]
                        << std::endl;
              }
              ++results.twoQubitGates;
            }
            ++bIt;
          }
          ++aIt;
        }
      }
      if (oldGateCount < results.singleQubitGates + results.twoQubitGates) {
        results.depth++;
      }
      auto tableau                  = results.resultTableaus.emplace_back();
      tableau                       = Tableau(localResultCircuit);
      results.resultTableaus.back() = tableau;
      modelTableau                  = Tableau(configuration.nqubits);
      for (int i = 0; i < configuration.nqubits; ++i) {
        modelTableau.populateTableauFrom(
            model->getBitvectorValue(x[gateStep][i], lb.get()),
            configuration.nqubits, i);
        modelTableau.populateTableauFrom(
            model->getBitvectorValue(z[gateStep][i], lb.get()),
            configuration.nqubits, i + configuration.nqubits);
      }
      modelTableau.populateTableauFrom(
          model->getBitvectorValue(r[gateStep], lb.get()),
          configuration.nqubits, 2 * configuration.nqubits);
      if (configuration.verbosity >= 5) {
        TRACE() << modelTableau;
      }
    }
    std::stringstream ss;
    localResultCircuit.dumpOpenQASM(ss);
    results.resultStringCircuit = ss.str();
    resultCircuit               = localResultCircuit.clone();
  }
  lb->reset();
  if (result == Result::SAT) {
    DEBUG() << "SAT" << std::endl;
    return results;
  }

  results.result = logicbase::Result::UNSAT;
  DEBUG() << "UNSAT" << std::endl;
  return results;
}

void CliffordSynthesizer::assertTableau(const SynthesisData& data,
                                        const Tableau&       tableau,
                                        std::size_t          position) {
  for (auto a = 0U; a < data.nqubits; ++a) {
    data.lb->assertFormula(
        data.x[position][a] ==
        logicbase::LogicTerm(tableau.getBVFrom(a), data.nqubits));
    data.lb->assertFormula(
        data.z[position][a] ==
        logicbase::LogicTerm(tableau.getBVFrom(a + data.nqubits),
                             data.nqubits));
  }
  data.lb->assertFormula(
      data.r[position] ==
      logicbase::LogicTerm(tableau.getBVFrom(2 * data.nqubits), data.nqubits));
}

void CliffordSynthesizer::initCouplingMaps(const Configuration& configuration) {
  if (configuration.chooseBest) {
    auto& cm = couplingMaps.emplace_back();
    configuration.architecture.getHighestFidelityCouplingMap(
        configuration.nqubits, cm);
  } else {
    configuration.architecture.getReducedCouplingMaps(configuration.nqubits,
                                                      couplingMaps);
  }
}
void CliffordSynthesizer::initResults() { optimalResults = Results(); }
} // namespace cs
