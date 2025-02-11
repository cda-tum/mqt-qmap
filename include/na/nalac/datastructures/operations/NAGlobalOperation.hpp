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
#include <vector>

namespace na::nalac {
class NAGlobalOperation : public NAOperation {
protected:
  qc::OpType opType;
  std::size_t ctrls;
  std::vector<qc::fp> params;

public:
  explicit NAGlobalOperation(const qc::OpType opType, const std::size_t ctrls,
                             const std::vector<qc::fp>& parameters)
      : opType(opType), ctrls(ctrls), params(parameters) {
    if (!isSingleQubitGate(opType)) {
      throw std::invalid_argument("Operation is not single qubit.");
    }
  }
  explicit NAGlobalOperation(const qc::OpType opType, const std::size_t ctrls)
      : NAGlobalOperation(opType, ctrls, {}) {}
  [[nodiscard]] auto getParams() const -> const std::vector<qc::fp>& {
    return params;
  }
  [[nodiscard]] auto getType() const -> std::pair<qc::OpType, std::size_t> {
    return {opType, ctrls};
  }
  [[nodiscard]] auto isGlobalOperation() const -> bool override { return true; }
  [[nodiscard]] auto toString() const -> std::string override;
  [[nodiscard]] auto clone() const -> std::unique_ptr<NAOperation> override {
    return std::make_unique<NAGlobalOperation>(*this);
  }
};
} // namespace na::nalac
