#pragma once

#include "NAOperation.hpp"
#include "na/Definitions.hpp"
#include "operations/OpType.hpp"

#include <cmath>
#include <utility>
#include <vector>
namespace na {

enum ShuttleType { LOAD, MOVE, STORE };

class NAShuttlingOperation : public NAOperation {
protected:
  ShuttleType                         type;
  std::vector<std::shared_ptr<Point>> start;
  std::vector<std::shared_ptr<Point>> end;

public:
  explicit NAShuttlingOperation(
      const ShuttleType type, const std::vector<std::shared_ptr<Point>>& start,
      const std::vector<std::shared_ptr<Point>>& end)
      : type(type), start(start), end(end) {
    if (start.size() != end.size()) {
      throw std::logic_error("Shuttling operation must have the same number of "
                             "start and end qubits.");
    }
  }
  explicit NAShuttlingOperation(const ShuttleType      type,
                                std::shared_ptr<Point> start,
                                std::shared_ptr<Point> end)
      : NAShuttlingOperation(
            type, std::vector<std::shared_ptr<Point>>{std::move(start)},
            std::vector<std::shared_ptr<Point>>{std::move(end)}) {}
  [[nodiscard]] auto getStart() const
      -> const std::vector<std::shared_ptr<Point>>& {
    return start;
  }
  [[nodiscard]] auto getEnd() const
      -> const std::vector<std::shared_ptr<Point>>& {
    return end;
  }
  static auto        isShuttlingOperation() -> bool { return true; }
  [[nodiscard]] auto toString() const -> std::string override {
    std::stringstream ss;
    switch (type) {
    case LOAD:
      ss << "load";
      break;
    case MOVE:
      ss << "move";
      break;
    case STORE:
      ss << "store";
      break;
    }
    ss << " ";
    for (const auto& p : start) {
      ss << *p << ", ";
    }
    ss.seekp(-2, std::ios_base::end);
    ss << " to ";
    for (const auto& p : end) {
      ss << *p << ", ";
    }
    ss.seekp(-2, std::ios_base::end);
    ss << ";" << std::endl;
    return ss.str();
  }
};
} // namespace na