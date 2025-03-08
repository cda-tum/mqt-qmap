#pragma once

namespace na {
class CodeGenerator {
  std::reference_wrapper<const na::Architecture> architecture_;

protected:
  CodeGenerator(const na::Architecture& architecture,
                const nlohmann::json& /* unused */)
      : architecture_(architecture) {}
  [[nodiscard]] auto generateCode(
      const std::vector<std::vector<const qc::Operation*>>& oneQubitGateLayers,
      const std::vector<std::vector<Site>>& placement,
      const std::vector<std::vector<std::vector<qc::Qubit>>>& routing) const
      -> na::NAComputation {
    std::cout << "Generating code with architecture \""
              << architecture_.get().name << "\"...\n";
    return {};
  }
};
} // namespace na
