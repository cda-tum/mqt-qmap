/*
 * Copyright (c) 2025 Chair for Design Automation, TUM
 * All rights reserved.
 *
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License
 */

#pragma once

#include "../NADefinitions.hpp"
#include "Definitions.hpp"
#include "NAOperation.hpp"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace na::nalac {
class NALocalOperation : public NAOperation {
protected:
  qc::OpType opType;
  std::size_t ctrls;
  std::vector<qc::fp> params;
  std::vector<std::shared_ptr<Point>> positions;

public:
  NALocalOperation(const qc::OpType opType, const std::size_t ctrls,
                   const std::vector<qc::fp>& parameter,
                   const std::vector<std::shared_ptr<Point>>& pos)
      : opType(opType), ctrls(ctrls), params(parameter), positions(pos) {
    if (!isSingleQubitGate(opType)) {
      throw std::invalid_argument("Operation is not single qubit.");
    }
    if (ctrls > 0) {
      throw std::logic_error("Control qubits are not supported.");
    }
  }
  explicit NALocalOperation(const qc::OpType opType, const std::size_t ctrls,
                            const std::vector<std::shared_ptr<Point>>& pos)
      : NALocalOperation(opType, ctrls, {}, pos) {}
  explicit NALocalOperation(const qc::OpType opType, const std::size_t ctrls,
                            const std::vector<qc::fp>& parameters,
                            std::shared_ptr<Point> pos)
      : NALocalOperation(opType, ctrls, parameters,
                         std::vector<std::shared_ptr<Point>>{std::move(pos)}) {}
  explicit NALocalOperation(const qc::OpType opType, const std::size_t ctrls,
                            std::shared_ptr<Point> pos)
      : NALocalOperation(opType, ctrls, {}, std::move(pos)) {}
  [[nodiscard]] auto getPositions() const
      -> const std::vector<std::shared_ptr<Point>>& {
    return positions;
  }
  [[nodiscard]] auto getParams() const -> const std::vector<qc::fp>& {
    return params;
  }
  [[nodiscard]] auto getType() const -> std::pair<qc::OpType, std::size_t> {
    return {opType, ctrls};
  }
  [[nodiscard]] auto isLocalOperation() const -> bool override { return true; }
  [[nodiscard]] auto toString() const -> std::string override;
  [[nodiscard]] auto clone() const -> std::unique_ptr<NAOperation> override {
    return std::make_unique<NALocalOperation>(*this);
  }
};
} // namespace na::nalac
