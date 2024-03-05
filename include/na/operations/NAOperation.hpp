#pragma once

#include <memory>
#include <string>
namespace na {
class NAOperation {
public:
  NAOperation() = default;
  NAOperation(const NAOperation& op) = default;
  NAOperation(NAOperation&& op) noexcept = default;
  NAOperation& operator=(const NAOperation& op) = default;
  NAOperation& operator=(NAOperation&& op) noexcept = default;
  virtual auto                isShuttlingOperation() -> bool { return false; }
  virtual auto                isLocalOperation() -> bool { return false; }
  virtual auto                isGlobalOperation() -> bool { return false; }
  [[nodiscard]] virtual auto toString() const -> std::string = 0;
  virtual ~NAOperation()                                     = default;
  [[nodiscard]] virtual auto clone() const -> std::unique_ptr<NAOperation> = 0;
};
} // namespace na