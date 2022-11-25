//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "LogicBlock/LogicBlock.hpp"
#include "cliffordsynthesis/Configuration.hpp"
#include "cliffordsynthesis/Results.hpp"
#include "cliffordsynthesis/Tableau.hpp"
#include "operations/OpType.hpp"

#include <cstddef>
#include <memory>

namespace cs::encoding {

using namespace logicbase;

class SATEncoder {
public:
  SATEncoder() = default;
  SATEncoder(const std::size_t nQubits, const std::size_t timestepLimit)
      : N(nQubits), T(timestepLimit) {}

  void   createFormulation(const Tableau&       initialTableau,
                           const Tableau&       targetTableau,
                           const Configuration& config);
  void   produceInstance();
  Result solve();
  void   extractResultsFromModel(Results& res);
  void   cleanup();

protected:
  void initializeSolver(const Configuration& config);

  void createTableauVariables();
  void createSingleQubitGateVariables();
  void createTwoQubitGateVariables();

  void asserTableau(const Tableau& tableau, std::size_t pos);

  void assertExactlyOne(const LogicVector& variables);
  void collectSingleQubitGateVariables(std::size_t pos, std::size_t qubit,
                                       LogicVector& variables);
  void collectTwoQubitGateVariables(std::size_t pos, std::size_t qubit,
                                    bool target, LogicVector& variables);

  LogicTerm createSingleQubitXChange(std::size_t pos, std::size_t qubit,
                                     qc::OpType gate);
  LogicTerm createSingleQubitZChange(std::size_t pos, std::size_t qubit,
                                     qc::OpType gate);
  LogicTerm createSingleQubitRChange(std::size_t pos, std::size_t qubit,
                                     qc::OpType gate);
  std::pair<LogicTerm, LogicTerm>
  createTwoQubitXChange(std::size_t pos, std::size_t ctrl, std::size_t trgt);
  std::pair<LogicTerm, LogicTerm>
  createTwoQubitZChange(std::size_t pos, std::size_t ctrl, std::size_t trgt);
  LogicTerm createTwoQubitRChange(std::size_t pos, std::size_t ctrl,
                                  std::size_t trgt);
  LogicTerm createNoChange(std::size_t pos, std::size_t except,
                           std::optional<std::size_t> except2);

  void      createIndividualGateEncoding();
  void      assertIndividualGateConsistency();
  void      assertIndividualGateTableauConstraints();
  void      assertIndividualGateSingleQubitGateConstraints(std::size_t pos);
  LogicTerm createIndividualGateSingleQubitGateConstraint(std::size_t pos,
                                                          std::size_t qubit,
                                                          qc::OpType  gate);
  void      assertIndividualGateTwoQubitGateConstraints(std::size_t pos);
  LogicTerm createIndividualGateTwoQubitGateConstraint(std::size_t pos,
                                                       std::size_t ctrl,
                                                       std::size_t trgt);

  void      createMultiGateEncoding();
  void      assertMultiGateConsistency();
  void      assertMultiGateTableauConstraints();
  void      assertMultiGateSingleQubitGateConstraints(std::size_t pos,
                                                      LogicTerm&  rChanges);
  LogicTerm createMultiGateSingleQubitGateConstraint(std::size_t pos,
                                                     std::size_t qubit,
                                                     qc::OpType  gate);
  void      assertMultiGateTwoQubitGateConstraints(std::size_t pos,
                                                   LogicTerm&  rChanges);
  LogicTerm createMultiGateTwoQubitGateConstraint(std::size_t pos,
                                                  std::size_t ctrl,
                                                  std::size_t trgt);

  void createObjectiveFunction(const Configuration& config);
  void createGateObjectiveFunction();
  void createDepthObjectiveFunction();

  void extractCircuitFromModel(Results& res, Model& model);
  void extractTableauFromModel(Results& res, std::size_t pos, Model& model);

  static constexpr std::array<qc::OpType, 7> SINGLE_QUBIT_GATES = {
      qc::OpType::None, qc::OpType::X, qc::OpType::Y,   qc::OpType::Z,
      qc::OpType::H,    qc::OpType::S, qc::OpType::Sdag};
  [[nodiscard]] static constexpr std::size_t
  gateToIndex(const qc::OpType type) {
    switch (type) {
    case qc::OpType::None:
      return 0;
    case qc::OpType::X:
      return 1;
    case qc::OpType::Y:
      return 2;
    case qc::OpType::Z:
      return 3;
    case qc::OpType::H:
      return 4;
    case qc::OpType::S:
      return 5;
    case qc::OpType::Sdag:
      return 6;
    default:
      return 0;
    }
  }

  struct Variables {
    // variables for the X parts of the tableaus
    LogicMatrix x{};
    // variables for the Z parts of the tableaus
    LogicMatrix z{};
    // variables for the phase parts of the tableaus
    LogicVector r{};
    // variables for the single-qubit gates
    LogicMatrix3D gS{};
    // variables for the two-qubit gates
    LogicMatrix3D gC{};
  };

  std::unique_ptr<LogicBlock> lb;
  Variables                   vars{};

  // number of qubits N
  const std::size_t N{};
  // timestep limit T
  const std::size_t T{};
};

} // namespace cs::encoding
