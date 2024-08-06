#include "nasp/Solver.hpp"

#include "Definitions.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// NOLINTNEXTLINE (misc-include-cleaner)
#include <yaml-cpp/yaml.h>
#include <z3++.h>

namespace na {
context NASolver::ctx{};

auto NASolver::minBitsToRepresentUInt(std::int32_t num) -> std::uint32_t {
  std::uint32_t bits = 0;
  while (num > 0) {
    num >>= 1;
    ++bits;
  }
  return bits;
}

auto NASolver::minBitsToRepresentInt(const std::int32_t num) -> std::uint32_t {
  return minBitsToRepresentUInt(num) + 1;
}

auto NASolver::initVariables() -> void {
  stages.clear();
  stages.reserve(numStages);
  for (std::uint16_t t = 0; t < numStages; ++t) {
    stages.emplace_back(t, numQubits, maxX, maxY, maxC, maxR, maxHOffset,
                        maxVOffset);
  }
  transfers.clear();
  if (numTransfers.has_value()) {
    transfers.reserve(numTransfers.value());
    for (std::uint16_t t = 0; t < numTransfers.value(); ++t) {
      transfers.emplace_back(
          ctx.bv_const(("transfer_" + std::to_string(t)).c_str(),
                       minBitsToRepresentUInt(numStages)));
    }
  } else {
    transfers.reserve(numStages);
    for (std::uint16_t t = 0; t < numStages; ++t) {
      transfers.emplace_back(
          ctx.bool_const(("transfer_" + std::to_string(t)).c_str()));
    }
  }
}

auto NASolver::getExactNumTransfersConstraints() const -> std::vector<expr> {
  std::vector<expr> constraints;
  if (numTransfers.has_value() && numTransfers.value() != 0) {
    constraints.reserve(numTransfers.value() + 1);
    // constraints.emplace_back(0 <= transfers[0]);
    for (std::uint16_t t = 1; t < numTransfers.value(); ++t) {
      constraints.emplace_back(ult(transfers[t - 1], transfers[t]));
    }
    constraints.emplace_back(
        ult(transfers[numTransfers.value() - 1],
            ctx.bv_val(numStages, minBitsToRepresentUInt(numStages))));
  }
  return constraints;
}

auto NASolver::getHaveSamePositionConstraint(
    const std::uint16_t q0, const std::uint16_t q1,
    const std::uint16_t t) const -> expr {
  return stages[t].getQubit(q0).getX() == stages[t].getQubit(q1).getX() &&
         stages[t].getQubit(q0).getY() == stages[t].getQubit(q1).getY();
}

auto NASolver::getHaveDifferentPositionConstraint(
    const std::uint16_t q0, const std::uint16_t q1,
    const std::uint16_t t) const -> expr {
  return !getHaveSamePositionConstraint(q0, q1, t);
}

// NOLINTNEXTLINE (bugprone-switch-missing-default-case)
auto NASolver::getAffectedByRydbergBeamConstraint(
    const std::uint16_t q, const std::uint16_t t) const -> expr {
  switch (storage) {
  case Storage::None:
    return ctx.bool_val(true);
  case Storage::Bottom:
    return ule(stages[t].getQubit(q).getY(),
               ctx.bv_val(maxEntanglingY, minBitsToRepresentUInt(maxY)));
  case Storage::TwoSided:
    return ule(minEntanglingY, stages[t].getQubit(q).getY()) &&
           ule(stages[t].getQubit(q).getY(),
               ctx.bv_val(maxEntanglingY, minBitsToRepresentUInt(maxY)));
  }
}

auto NASolver::getShieldedFromRydbergBeamConstraint(
    const std::uint16_t q, const std::uint16_t t) const -> expr {
  return !getAffectedByRydbergBeamConstraint(q, t);
}

auto NASolver::getValidRydbergTransitionConstraints(const std::uint16_t t) const
    -> std::vector<expr> {
  if (t == numStages - 1) {
    std::stringstream msg;
    msg << "There is no next stage after the last stage " << t;
    throw std::invalid_argument(msg.str());
  }
  std::vector<expr> constraints;
  constraints.reserve(3UL * numQubits);
  for (std::uint16_t i = 0; i < numQubits; ++i) {
    // For all: AOD/SLM state is preserved
    constraints.emplace_back(implies(getRydbergStageConstraint(t),
                                     stages[t].getQubit(i).getA() ==
                                         stages[t + 1].getQubit(i).getA()));
    // For AOD: column and row are preserved
    constraints.emplace_back(implies(
        getRydbergStageConstraint(t) && stages[t].getQubit(i).getA(),
        stages[t].getQubit(i).getC() == stages[t + 1].getQubit(i).getC() &&
            stages[t].getQubit(i).getR() == stages[t + 1].getQubit(i).getR()));
    // For SLM: stay at the same position
    constraints.emplace_back(implies(
        getRydbergStageConstraint(t) && !stages[t].getQubit(i).getA(),
        stages[t].getQubit(i).getX() == stages[t + 1].getQubit(i).getX() &&
            stages[t].getQubit(i).getY() == stages[t + 1].getQubit(i).getY()));
  }
  // we set load and store variables to false because they do not carray any
  // meaning in the rydberg stage anyways
  for (std::uint16_t i = 0; i <= static_cast<std::uint16_t>(maxC); ++i) {
    constraints.emplace_back(
        implies(getRydbergStageConstraint(t), !stages[t].getLoadCol(i)));
    constraints.emplace_back(
        implies(getRydbergStageConstraint(t), !stages[t].getStoreCol(i)));
  }
  for (std::uint16_t i = 0; i <= static_cast<std::uint16_t>(maxR); ++i) {
    constraints.emplace_back(
        implies(getRydbergStageConstraint(t), !stages[t].getLoadRow(i)));
    constraints.emplace_back(
        implies(getRydbergStageConstraint(t), !stages[t].getStoreRow(i)));
  }
  return constraints;
}

auto NASolver::getValidTransferTransitionConstraints(
    const std::uint16_t t) const -> std::vector<expr> {
  if (t == numStages - 1) {
    std::stringstream msg;
    msg << "There is no next stage after the last stage " << t;
    throw std::invalid_argument(msg.str());
  }
  std::vector<expr> constraints;
  // TODO: check all reserves for the correct capacity
  constraints.reserve(3UL * numQubits);
  for (std::uint16_t i = 0; i < numQubits; ++i) {
    // For SLM and Stored: stay at the same position
    constraints.emplace_back(implies(
        getTransferStageConstraint(t) && !stages[t + 1].getQubit(i).getA(),
        stages[t].getQubit(i).getX() == stages[t + 1].getQubit(i).getX() &&
            stages[t].getQubit(i).getY() == stages[t + 1].getQubit(i).getY()));
    // For Stored and Loaded: Must be at zero horizontal and vertical offset
    constraints.emplace_back(implies(getTransferStageConstraint(t) &&
                                         stages[t].getQubit(i).getA() !=
                                             stages[t + 1].getQubit(i).getA(),
                                     0 == stages[t].getQubit(i).getH() &&
                                         0 == stages[t].getQubit(i).getV()));
    // For Loaded: Entire row or column must be loaded
    {
      expr colClauses = ctx.bool_val(true);
      for (std::uint16_t c = 0; c <= static_cast<std::uint16_t>(maxC); ++c) {
        colClauses = colClauses &&
                     implies(stages[t].getQubit(i).getC() ==
                                 ctx.bv_val(c, minBitsToRepresentUInt(maxC)),
                             stages[t].getLoadCol(c));
      }
      expr rowClauses = ctx.bool_val(true);
      for (std::uint16_t r = 0; r <= static_cast<std::uint16_t>(maxR); ++r) {
        rowClauses = rowClauses &&
                     implies(stages[t].getQubit(i).getR() ==
                                 ctx.bv_val(r, minBitsToRepresentUInt(maxR)),
                             stages[t].getLoadRow(r));
      }
      constraints.emplace_back(implies(
          getTransferStageConstraint(t) && !stages[t].getQubit(i).getA(),
          stages[t + 1].getQubit(i).getA() == (colClauses || rowClauses)));
    }
    // For Stored: Entire row or column must be stored
    {
      expr colClauses = ctx.bool_val(true);
      for (std::uint16_t c = 0; c <= static_cast<std::uint16_t>(maxC); ++c) {
        colClauses = colClauses &&
                     implies(stages[t].getQubit(i).getC() ==
                                 ctx.bv_val(c, minBitsToRepresentUInt(maxC)),
                             stages[t].getStoreCol(c));
      }
      expr rowClauses = ctx.bool_val(true);
      for (std::uint16_t r = 0; r <= static_cast<std::uint16_t>(maxR); ++r) {
        rowClauses = rowClauses &&
                     implies(stages[t].getQubit(i).getR() ==
                                 ctx.bv_val(r, minBitsToRepresentUInt(maxR)),
                             stages[t].getStoreRow(r));
      }
      constraints.emplace_back(implies(
          getTransferStageConstraint(t) && stages[t].getQubit(i).getA(),
          !stages[t + 1].getQubit(i).getA() == (colClauses || rowClauses)));
    }
    for (std::uint16_t j = 0; j < numQubits; ++j) {
      // For Loaded: Loaded atoms remain in the relative position with respect
      // to other AOD/loaded atoms
      constraints.emplace_back(implies(
          getTransferStageConstraint(t) && stages[t + 1].getQubit(i).getA() &&
              stages[t + 1].getQubit(j).getA(),
          (ult(stages[t].getQubit(i).getX(), stages[t].getQubit(j).getX()) ||
           ((stages[t].getQubit(i).getX() == stages[t].getQubit(j).getX()) &&
            slt(stages[t].getQubit(i).getH(), stages[t].getQubit(j).getH()))) ==
              ult(stages[t + 1].getQubit(i).getC(),
                  stages[t + 1].getQubit(j).getC())));
      constraints.emplace_back(implies(
          getTransferStageConstraint(t) && stages[t + 1].getQubit(i).getA() &&
              stages[t + 1].getQubit(j).getA(),
          (ult(stages[t].getQubit(i).getY(), stages[t].getQubit(j).getY()) ||
           ((stages[t].getQubit(i).getY() == stages[t].getQubit(j).getY()) &&
            slt(stages[t].getQubit(i).getV(), stages[t].getQubit(j).getV()))) ==
              ult(stages[t + 1].getQubit(i).getR(),
                  stages[t + 1].getQubit(j).getR())));
    }
  }
  return constraints;
}

auto NASolver::getCircuitExecutionConstraints(
    const std::vector<std::pair<qc::Qubit, qc::Qubit>>& ops,
    const bool mindOpsOrder, const bool shieldIdleAtoms) -> std::vector<expr> {
  const auto numGates = static_cast<std::uint16_t>(ops.size());
  std::unordered_map<std::pair<qc::Qubit, qc::Qubit>, std::vector<expr>,
                     qc::PairHash<qc::Qubit, qc::Qubit>>
      pairToGates;
  gates.clear();
  gates.reserve(numGates);
  std::vector<std::vector<expr>> gatesForQubit(numQubits);
  for (std::uint16_t i = 0; i < numGates; ++i) {
    gates.emplace_back(ctx.bv_const(("gate_" + std::to_string(i)).c_str(),
                                    minBitsToRepresentUInt(numStages)));
    const auto key = std::make_pair(std::min(ops[i].first, ops[i].second),
                                    std::max(ops[i].first, ops[i].second));
    if (pairToGates.find(key) == pairToGates.end()) {
      pairToGates.emplace(key, std::vector<expr>());
    }
    pairToGates.at(key).emplace_back(gates.back());
    gatesForQubit[ops[i].first].emplace_back(gates.back());
    gatesForQubit[ops[i].second].emplace_back(gates.back());
  }
  std::vector<expr> constraints;
  if (mindOpsOrder) {
    std::vector<std::optional<const expr>> lastGateOnQubit(numQubits,
                                                           std::nullopt);
    for (std::uint16_t i = 0; i < numGates; ++i) {
      const auto& gate = gates[i];
      // if (!lastGateOnQubit[ops[i].first].has_value() &&
      // !lastGateOnQubit[ops[i].
      //       second].has_value()) {
      //   constraints.emplace_back(
      //       0 <= gate);
      // }
      if (lastGateOnQubit[ops[i].first].has_value()) {
        constraints.emplace_back(
            // NOLINTNEXTLINE (bugprone-unchecked-optional-access)
            ult(lastGateOnQubit[ops[i].first].value(), gate));
      }
      if (lastGateOnQubit[ops[i].second].has_value()) {
        constraints.emplace_back(
            // NOLINTNEXTLINE (bugprone-unchecked-optional-access)
            ult(lastGateOnQubit[ops[i].second].value(), gate));
      }
      lastGateOnQubit[ops[i].first].emplace(gate);
      lastGateOnQubit[ops[i].second].emplace(gate);
    }
    std::unordered_set<expr, ExprHash> lastGates;
    for (const auto& lastGate : lastGateOnQubit) {
      if (lastGate.has_value()) {
        lastGates.emplace(lastGate.value());
      }
    }
    for (const auto& lastGate : lastGates) {
      constraints.emplace_back(ult(lastGate, numStages));
    }
    constraints.reserve(3 * numGates +
                        numStages * numQubits * (numQubits - 1) / 2);
  } else {
    constraints.reserve(numGates * (numGates - 1) / 2 + 4 * numGates +
                        numStages * numQubits * (numQubits - 1) / 2);
    for (std::uint16_t g = 0; g < numGates; ++g) {
      constraints.emplace_back(
          // 0 <= gates[g] &&
          ult(gates[g], numStages));
      for (std::uint16_t h = g + 1; h < numGates; ++h) {
        if (ops[g].first == ops[h].first || ops[g].first == ops[h].second ||
            ops[g].second == ops[h].first || ops[g].second == ops[h].second) {
          constraints.emplace_back(gates[g] != gates[h]);
        }
      }
    }
  }
  for (std::uint16_t t = 0; t < numStages; ++t) {
    for (std::uint16_t i = 0; i < numQubits; ++i) {
      for (std::uint16_t j = i + 1; j < numQubits; ++j) {
        if (const auto& it = pairToGates.find({i, j});
            it != pairToGates.end() && !it->second.empty()) {
          const auto& gatesForPair = it->second;
          for (const auto& gate : gatesForPair) {
            auto hDiff =
                stages[t].getQubit(j).getH() - stages[t].getQubit(i).getH();
            auto hDiffSign = ashr(
                hDiff,
                static_cast<signed int>(minBitsToRepresentInt(maxHOffset)) - 1);
            auto absHDiff = (hDiff ^ hDiffSign) - hDiffSign;
            auto vDiff =
                stages[t].getQubit(j).getV() - stages[t].getQubit(i).getV();
            auto vDiffSign = ashr(
                vDiff,
                static_cast<signed int>(minBitsToRepresentInt(maxVOffset)) - 1);
            auto absVDiff = (vDiff ^ vDiffSign) - vDiffSign;

            constraints.emplace_back(implies(
                gate == ctx.bv_val(t, minBitsToRepresentUInt(numStages)),
                getRydbergStageConstraint(t) &&
                    getHaveSamePositionConstraint(i, j, t) &&
                    getAffectedByRydbergBeamConstraint(i, t) &&
                    getAffectedByRydbergBeamConstraint(j, t) &&
                    ult(absHDiff, maxHDist) &&
                    ult(absVDiff, ctx.bv_val(maxVDist, minBitsToRepresentInt(
                                                           maxVOffset)))));
          }
          expr premisses = getRydbergStageConstraint(t) &&
                           getAffectedByRydbergBeamConstraint(i, t) &&
                           getAffectedByRydbergBeamConstraint(j, t);
          for (const auto& gate : gatesForPair) {
            premisses =
                premisses &&
                gate != ctx.bv_val(t, minBitsToRepresentUInt(numStages));
          }
          constraints.emplace_back(
              implies(premisses, getHaveDifferentPositionConstraint(i, j, t)));
        } else {
          constraints.emplace_back(
              implies(getRydbergStageConstraint(t) &&
                          getAffectedByRydbergBeamConstraint(i, t) &&
                          getAffectedByRydbergBeamConstraint(j, t),
                      getHaveDifferentPositionConstraint(i, j, t)));
        }
      }
      if (shieldIdleAtoms) {
        if (gatesForQubit[i].empty()) {
          constraints.emplace_back(
              implies(getRydbergStageConstraint(t),
                      getShieldedFromRydbergBeamConstraint(i, t)));
        } else {
          expr premisses = getRydbergStageConstraint(t);
          for (const auto& gate : gatesForQubit[i]) {
            premisses =
                premisses &&
                gate != ctx.bv_val(t, minBitsToRepresentUInt(numStages));
          }
          constraints.emplace_back(
              implies(premisses, getShieldedFromRydbergBeamConstraint(i, t)));
        }
      }
    }
  }
  return constraints;
}

auto NASolver::getRydbergStageConstraint(const std::uint16_t t) const -> expr {
  return !getTransferStageConstraint(t);
}

auto NASolver::getTransferStageConstraint(const std::uint16_t t) const -> expr {
  if (numTransfers.has_value()) {
    expr clauses = ctx.bool_val(false);
    for (const auto& transfer : transfers) {
      clauses = clauses ||
                transfer == ctx.bv_val(t, minBitsToRepresentUInt(numStages));
    }
    return clauses;
  }
  return transfers[t];
}

auto NASolver::getValidStageConstraints(const std::uint16_t t) const
    -> std::vector<expr> {
  std::vector<expr> constraints;
  constraints.reserve(numQubits * (5UL * numQubits + 3UL));
  for (std::uint16_t i = 0; i < numQubits; ++i) {
    // 0 <= x <= maxX
    constraints.emplace_back(
        // ule(0, stages[t].getQubit(i).
        // getX()) &&
        ule(stages[t].getQubit(i).getX(),
            ctx.bv_val(maxX, minBitsToRepresentUInt(maxX))));
    // 0 <= y <= maxY
    constraints.emplace_back(
        // 0 <= stages[t].getQubit(i).
        // getY() &&
        ule(stages[t].getQubit(i).getY(),
            ctx.bv_val(maxY, minBitsToRepresentUInt(maxY))));
    // For AOD: 0 <= c <= maxC
    constraints.emplace_back(
        implies(stages[t].getQubit(i).getA(),
                // 0 <= stages[t].getQubit(i).
                // getC() &&
                ule(stages[t].getQubit(i).getC(),
                    ctx.bv_val(maxC, minBitsToRepresentUInt(maxC)))));
    // For AOD: 0 <= r <= maxR
    constraints.emplace_back(
        implies(stages[t].getQubit(i).getA(),
                // 0 <= stages[t].getQubit(i).
                // getR() &&
                ule(stages[t].getQubit(i).getR(),
                    ctx.bv_val(maxR, minBitsToRepresentUInt(maxR)))));
    // For AOD: - maxHOffset <= h <= maxHOffset
    constraints.emplace_back(
        implies(stages[t].getQubit(i).getA(),
                sle(ctx.bv_val(-static_cast<signed int>(maxHOffset),
                               minBitsToRepresentInt(maxHOffset)),
                    stages[t].getQubit(i).getH()) &&
                    sle(stages[t].getQubit(i).getH(), maxHOffset)));
    // For AOD: - maxVOffset <= v <= maxVOffset
    constraints.emplace_back(
        implies(stages[t].getQubit(i).getA(),
                sle(ctx.bv_val(-static_cast<signed int>(maxVOffset),
                               minBitsToRepresentInt(maxVOffset)),
                    stages[t].getQubit(i).getV()) &&
                    sle(stages[t].getQubit(i).getV(), maxVOffset)));
    // For SLM: c = 0, r = 0, h = 0, v = 0
    constraints.emplace_back(
        implies(!stages[t].getQubit(i).getA(),
                stages[t].getQubit(i).getC() ==
                        ctx.bv_val(0, minBitsToRepresentUInt(maxC)) &&
                    stages[t].getQubit(i).getR() ==
                        ctx.bv_val(0, minBitsToRepresentUInt(maxR)) &&
                    stages[t].getQubit(i).getH() ==
                        ctx.bv_val(0, minBitsToRepresentInt(maxHOffset)) &&
                    stages[t].getQubit(i).getV() ==
                        ctx.bv_val(0, minBitsToRepresentInt(maxVOffset))));
    for (std::uint16_t j = 0; j < numQubits; ++j) {
      // For AOD: (x < x' ∨ x == x' ∧ v < v') <-> r < r'
      constraints.emplace_back(implies(
          stages[t].getQubit(i).getA() && stages[t].getQubit(j).getA(),
          (ult(stages[t].getQubit(i).getX(), stages[t].getQubit(j).getX()) ||
           ((stages[t].getQubit(i).getX() == stages[t].getQubit(j).getX()) &&
            slt(stages[t].getQubit(i).getH(), stages[t].getQubit(j).getH()))) ==
              ult(stages[t].getQubit(i).getC(), stages[t].getQubit(j).getC())));
      // For AOD: (y < y' ∨ y == y' ∧ h < h') <-> c < c'
      constraints.emplace_back(implies(
          stages[t].getQubit(i).getA() && stages[t].getQubit(j).getA(),
          (ult(stages[t].getQubit(i).getY(), stages[t].getQubit(j).getY()) ||
           ((stages[t].getQubit(i).getY() == stages[t].getQubit(j).getY()) &&
            slt(stages[t].getQubit(i).getV(), stages[t].getQubit(j).getV()))) ==
              ult(stages[t].getQubit(i).getR(), stages[t].getQubit(j).getR())));
    }
    for (std::uint16_t j = i + 1; j < numQubits; ++j) {
      // For all: h == h' -> v == v' -> x != x' or y != y'
      constraints.emplace_back(implies(
          stages[t].getQubit(i).getH() == stages[t].getQubit(j).getH() &&
              stages[t].getQubit(i).getV() == stages[t].getQubit(j).getV(),
          getHaveDifferentPositionConstraint(i, j, t)));
    }
  }
  return constraints;
}

auto NASolver::init(const std::uint16_t newMaxX, const std::uint16_t newMaxY,
                    const std::uint16_t newMaxC, const std::uint16_t newMaxR,
                    const std::uint16_t newMaxHOffset,
                    const std::uint16_t newMaxVOffset,
                    const std::uint16_t newMaxHDist,
                    const std::uint16_t newMaxVDist,
                    const std::uint16_t newMinEntanglingY,
                    const std::uint16_t newMaxEntanglingY) -> void {
  maxX           = newMaxX;
  maxY           = newMaxY;
  minEntanglingY = newMinEntanglingY;
  maxEntanglingY = newMaxEntanglingY;
  maxC           = newMaxC;
  maxR           = newMaxR;
  maxHOffset     = newMaxHOffset;
  maxVOffset     = newMaxVOffset;
  maxHDist       = newMaxHDist;
  maxVDist       = newMaxVDist;
  if (minEntanglingY == 0 && maxEntanglingY < maxY) {
    storage = Storage::Bottom;
  } else if (minEntanglingY > 0 && maxEntanglingY < maxY) {
    storage = Storage::TwoSided;
  } else if (minEntanglingY == 0 && maxEntanglingY == maxY) {
    storage = Storage::None;
  } else {
    throw std::invalid_argument("One sided storage zone is only supported "
                                "below the entangling zone (higher Y).");
  }
}

auto NASolver::solve(const std::vector<std::pair<qc::Qubit, qc::Qubit>>& ops,
                     const std::uint16_t newNumQubits,
                     const std::uint16_t newNumStages,
                     const std::uint16_t newNumTransfers,
                     const bool          mindOpsOrder,
                     const bool          shieldIdleQubits) -> Result {
  if (shieldIdleQubits) {
    if (storage == Storage::None) {
      throw std::invalid_argument("No storage zone is available.");
    }
  }

  numQubits    = newNumQubits;
  numStages    = newNumStages;
  numTransfers = newNumTransfers;

  // solver solver(ctx, "QF_LIA");
  solver solver(ctx, "QF_BV");

  initVariables();

  // Now we assert the constraints to the solver.
  for (const auto& c : getExactNumTransfersConstraints()) {
    solver.add(c);
  }
  for (const auto& c :
       getCircuitExecutionConstraints(ops, mindOpsOrder, shieldIdleQubits)) {
    solver.add(c);
  }
  for (std::uint16_t t = 0; t < numStages; ++t) {
    for (const auto& c : getValidStageConstraints(t)) {
      solver.add(c);
    }
    if (t < numStages - 1) {
      for (const auto& c : getValidRydbergTransitionConstraints(t)) {
        solver.add(c);
      }
      for (const auto& c : getValidTransferTransitionConstraints(t)) {
        solver.add(c);
      }
    }
  }

  // Check satisfiability
  if (solver.check() == unsat) {
    return Result(false);
  }
  const auto                 model  = solver.get_model();
  std::uint16_t              nTrans = 0;
  std::vector<Result::Stage> resultStages;
  resultStages.reserve(numStages);
  for (const auto& stage : stages) {
    const bool rydberg =
        nTrans == numTransfers ||
        model.eval(transfers[nTrans]).as_uint64() != stage.getT();
    if (!rydberg) {
      ++nTrans;
    }
    std::vector<Result::Qubit> resultQubits;
    resultQubits.reserve(numQubits);
    for (std::uint16_t i = 0; i < numQubits; ++i) {
      resultQubits.emplace_back(
          model.eval(stage.getQubit(i).getX()).as_uint64(),
          model.eval(stage.getQubit(i).getY()).as_uint64(),
          model.eval(stage.getQubit(i).getA()).is_true(),
          model.eval(stage.getQubit(i).getC()).as_uint64(),
          model.eval(stage.getQubit(i).getR()).as_uint64(),
          model.eval(bv2int(stage.getQubit(i).getH(), true)).get_numeral_int(),
          model.eval(bv2int(stage.getQubit(i).getV(), true)).get_numeral_int());
    }
    std::vector<Result::Gate> resultGates;
    for (std::uint16_t i = 0; i < static_cast<std::uint16_t>(gates.size());
         ++i) {
      if (model.eval(gates[i]).as_uint64() == stage.getT()) {
        resultGates.emplace_back(stage.getT(), ops[i]);
      }
    }
    resultStages.emplace_back(rydberg, resultQubits, resultGates);
  }
  return {true, resultStages};
}

auto NASolver::solve(const std::vector<std::pair<qc::Qubit, qc::Qubit>>& ops,
                     const std::uint16_t newNumQubits,
                     const std::uint16_t newNumStages, const bool mindOpsOrder,
                     const bool shieldIdleQubits) -> Result {
  // this is almost the same as the method above just that it does not fix the
  // number of transfers which implies few minor changes marked in the following
  if (shieldIdleQubits) {
    if (storage == Storage::None) {
      throw std::invalid_argument("No storage zone is available.");
    }
  }

  numQubits    = newNumQubits;
  numStages    = newNumStages;
  numTransfers = std::nullopt;
  // CHANGE: instead of a fixed number of transfers

  // solver solver(ctx, "QF_LIA");
  solver solver(ctx, "QF_BV");

  initVariables();

  // Now we assert the constraints to the solver.
  // CHANGE: no exact number of transfers constraint
  for (const auto& c :
       getCircuitExecutionConstraints(ops, mindOpsOrder, shieldIdleQubits)) {
    solver.add(c);
  }
  for (std::uint16_t t = 0; t < numStages; ++t) {
    for (const auto& c : getValidStageConstraints(t)) {
      solver.add(c);
    }
    if (t < numStages - 1) {
      for (const auto& c : getValidRydbergTransitionConstraints(t)) {
        solver.add(c);
      }
      for (const auto& c : getValidTransferTransitionConstraints(t)) {
        solver.add(c);
      }
    }
  }

  // Check satisfiability
  if (solver.check() == unsat) {
    return Result(false);
  }
  const auto                 model = solver.get_model();
  std::vector<Result::Stage> resultStages;
  resultStages.reserve(numStages);
  for (const auto& stage : stages) {
    const bool rydberg = model.eval(transfers[stage.getT()]).is_false();
    // CHANGE: Read boolan variable directly
    std::vector<Result::Qubit> resultQubits;
    resultQubits.reserve(numQubits);
    for (std::uint16_t i = 0; i < numQubits; ++i) {
      resultQubits.emplace_back(
          model.eval(stage.getQubit(i).getX()).get_numeral_uint(),
          model.eval(stage.getQubit(i).getY()).get_numeral_uint(),
          model.eval(stage.getQubit(i).getA()).is_true(),
          model.eval(stage.getQubit(i).getC()).get_numeral_uint(),
          model.eval(stage.getQubit(i).getR()).get_numeral_uint(),
          model.eval(bv2int(stage.getQubit(i).getH(), true)).get_numeral_int(),
          model.eval(bv2int(stage.getQubit(i).getV(), true)).get_numeral_int());
    }
    std::vector<Result::Gate> resultGates;
    for (std::uint16_t i = 0; i < static_cast<std::uint16_t>(gates.size());
         ++i) {
      if (model.eval(gates[i]).as_uint64() == stage.getT()) {
        resultGates.emplace_back(stage.getT(), ops[i]);
      }
    }
    resultStages.emplace_back(rydberg, resultQubits, resultGates);
  }
  return {true, resultStages};
}

/// Initialize a Qubit from a YAML string.
// NOLINTNEXTLINE (misc-include-cleaner)
auto NASolver::Result::Qubit::fromYAML(const YAML::Node& yaml) -> Qubit {
  Qubit qubit{};
  qubit.a = yaml["a"].as<bool>();
  qubit.c = yaml["c"].as<std::uint16_t>();
  qubit.h = yaml["h"].as<int>();
  qubit.r = yaml["r"].as<std::uint16_t>();
  qubit.v = yaml["v"].as<int>();
  qubit.x = yaml["x"].as<std::uint16_t>();
  qubit.y = yaml["y"].as<std::uint16_t>();
  return qubit;
}

auto NASolver::Result::Qubit::yaml(std::size_t indent, const bool item,
                                   const bool compact) const -> std::string {
  std::stringstream ss;
  ss << std::boolalpha;
  ss << std::string(indent, ' ');
  if (item) {
    ss << "- ";
    indent += 2;
  }
  if (compact) {
    ss << "{x: " << x << ", y: " << y << ", a: " << a << ", c: " << c
       << ", r: " << r << ", h: " << h << ", v: " << v << "}\n";
    return ss.str();
  }
  ss << "x: " << x << "\n";
  ss << std::string(indent, ' ') << "y: " << y << "\n";
  ss << std::string(indent, ' ') << "a: " << a << "\n";
  ss << std::string(indent, ' ') << "c: " << c << "\n";
  ss << std::string(indent, ' ') << "r: " << r << "\n";
  ss << std::string(indent, ' ') << "h: " << h << "\n";
  ss << std::string(indent, ' ') << "v: " << v << "\n";
  return ss.str();
}
auto NASolver::Result::Qubit::operator==(const Qubit& other) const -> bool {
  return x == other.x && y == other.y && a == other.a && c == other.c &&
         r == other.r && h == other.h && v == other.v;
}

auto NASolver::Result::Stage::yaml(std::size_t indent, const bool item,
                                   const bool compact) const -> std::string {
  std::stringstream ss;
  ss << std::boolalpha;
  ss << std::string(indent, ' ');
  if (item) {
    ss << "- ";
    indent += 2;
  }
  ss << "rydberg: " << rydberg << "\n";
  ss << std::string(indent, ' ') << "gates:\n";
  for (const auto& gate : gates) {
    ss << gate.yaml(indent + 2, true, compact);
  }
  ss << std::string(indent, ' ') << "qubits:\n";
  for (const auto& qubit : qubits) {
    ss << qubit.yaml(indent + 2, true, compact);
  }
  return ss.str();
}
auto NASolver::Result::Stage::operator==(const Stage& other) const -> bool {
  return rydberg == other.rydberg && gates == other.gates &&
         qubits == other.qubits;
}

auto NASolver::Result::fromYAML(const YAML::Node& yaml) -> Result {
  Result result{};
  result.sat = yaml["sat"].as<bool>();
  if (result.sat) {
    for (const auto& stage : yaml["stages"]) {
      result.stages.emplace_back(Stage::fromYAML(stage));
    }
  }
  return result;
}

auto NASolver::Result::Gate::fromYAML(const YAML::Node& yaml) -> Gate {
  Gate gate{};
  gate.qubits.first  = yaml["qubits"][0].as<qc::Qubit>();
  gate.qubits.second = yaml["qubits"][1].as<qc::Qubit>();
  return gate;
}

auto NASolver::Result::Gate::yaml(std::size_t indent, const bool item,
                                  const bool compact) const -> std::string {
  std::stringstream ss;
  ss << std::string(indent, ' ');
  if (item) {
    ss << "- ";
    indent += 2;
  }
  if (compact) {
    ss << "qubits: [" << qubits.first << ", " << qubits.second << "]\n";
    return ss.str();
  }
  ss << "qubits:\n";
  ss << std::string(indent, ' ') << "- " << qubits.first << "\n";
  ss << std::string(indent, ' ') << "- " << qubits.second << "\n";
  return ss.str();
}
auto NASolver::Result::Gate::operator==(const Gate& other) const -> bool {
  return qubits == other.qubits;
}

auto NASolver::Result::Stage::fromYAML(const YAML::Node& yaml) -> Stage {
  Stage stage{};
  stage.rydberg = yaml["rydberg"].as<bool>();
  for (const auto& gate : yaml["gates"]) {
    stage.gates.emplace_back(Gate::fromYAML(gate));
  }
  for (const auto& qubit : yaml["qubits"]) {
    stage.qubits.emplace_back(Qubit::fromYAML(qubit));
  }
  return stage;
}

auto NASolver::Result::yaml(const std::size_t indent,
                            const bool        compact) const -> std::string {
  std::stringstream ss;
  ss << std::boolalpha;
  ss << std::string(indent, ' ') << "sat: " << sat << "\n";
  if (sat) {
    ss << std::string(indent, ' ') << "stages:\n";
    for (const auto& stage : stages) {
      ss << stage.yaml(indent + 2, true, compact);
    }
  }
  return ss.str();
}
auto NASolver::Result::operator==(const Result& other) const -> bool {
  return sat == other.sat && stages == other.stages;
}
} // namespace na
