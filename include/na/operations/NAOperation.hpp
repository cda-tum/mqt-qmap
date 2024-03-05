#pragma once

#include <string>
namespace na {
class NAOperation {
public:
  NAOperation() = default;
  NAOperation(const NAOperation& op) = default;
  NAOperation(NAOperation&& op) noexcept = default;
  NAOperation& operator=(const NAOperation& op) = default;
  NAOperation& operator=(NAOperation&& op) noexcept = default;
  static auto                isShuttlingOperation() -> bool { return false; }
  static auto                isLocalOperation() -> bool { return false; }
  static auto                isGlobalOperation() -> bool { return false; }
  [[nodiscard]] virtual auto toString() const -> std::string = 0;
  virtual ~NAOperation()                                     = default;
};
} // namespace na