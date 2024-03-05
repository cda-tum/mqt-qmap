#pragma once

#include "NAOperation.hpp"
#include "na/Definitions.hpp"
#include "operations/OpType.hpp"

#include <cmath>
#include <utility>
namespace na {
class NALocalOperation : public NAOperation {
protected:
  OpType                               type;
  std::vector<qc::fp>                  params;
  std::vector<std::shared_ptr<Point>> positions;

public:
  explicit NALocalOperation(
      const OpType& type, const std::vector<qc::fp>& params,
      const std::vector<std::shared_ptr<Point>>& positions)
      : type(type), params(params), positions(positions) {
    if (!isSingleQubitGate(type.type)) {
      throw std::invalid_argument("Operation is not single qubit.");
    }
    if (type.nctrl != 0) {
      throw std::logic_error("Operation control qubits are not supported.");
    }
  }
  explicit NALocalOperation(const OpType&                    type,
                            const std::vector<std::shared_ptr<Point>>& positions)
      : NALocalOperation(type, {}, positions) {}
  explicit NALocalOperation(const OpType&              type,
                            const std::vector<qc::fp>& params,
                            std::shared_ptr<Point>               position)
      : NALocalOperation(type, params, std::vector<std::shared_ptr<Point>>{std::move(position)}) {}
  explicit NALocalOperation(const OpType& type, std::shared_ptr<Point> position)
      : NALocalOperation(type, {}, std::move(position)) {}
  [[nodiscard]] auto getPositions() const -> std::vector<std::shared_ptr<Point>> {
    return positions;
  }
  [[nodiscard]] auto getParams() const -> std::vector<qc::fp> {
    return params;
  }
  auto        isLocalOperation() -> bool override { return true; }
  [[nodiscard]] auto toString() const -> std::string override {
    std::stringstream ss;
    ss << type;
    if (!params.empty()) {
      ss << "(";
      for (const auto& p : params) {
        ss << p << ", ";
      }
      ss.seekp(-2, std::ios_base::end);
      ss << ")";
    }
    ss << " at ";
    for (const auto& p : positions) {
      ss << *p << ", ";
    }
    ss.seekp(-2, std::ios_base::end);
    ss << ";" << std::endl;
    return ss.str();
  }
  [[nodiscard]] auto clone() const -> std::unique_ptr<NAOperation> override {
    return std::make_unique<NALocalOperation>(*this);
  }
};
} // namespace na