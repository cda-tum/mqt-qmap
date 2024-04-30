//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include <cstddef>
#include <cstdint>
#include <istream>
#include <string>
#include <unordered_map>
namespace na {

enum class NaMappingMethod : std::uint8_t { NAIVE, SMART };
static const std::unordered_map<std::string, NaMappingMethod> STRING_TO_METHOD =
    {{"Naive", NaMappingMethod::NAIVE}, {"Smart", NaMappingMethod::SMART},
     {"naive", NaMappingMethod::NAIVE}, {"smart", NaMappingMethod::SMART},
     {"NAIVE", NaMappingMethod::NAIVE}, {"SMART", NaMappingMethod::SMART}};
[[nodiscard]] inline auto
getMethodOfString(const std::string& method) -> NaMappingMethod;
class Configuration {
private:
  std::size_t     patchRows = 1;
  std::size_t     patchCols = 1;
  NaMappingMethod method    = NaMappingMethod::SMART;

public:
  Configuration() = default;
  explicit Configuration(const NaMappingMethod mappingMethod)
      : method(mappingMethod) {};
  explicit Configuration(const std::size_t rows, const std::size_t cols)
      : patchRows(rows), patchCols(cols) {};
  explicit Configuration(const std::size_t rows, const std::size_t cols,
                         const NaMappingMethod mappingMethod)
      : patchRows(rows), patchCols(cols), method(mappingMethod) {};
  explicit Configuration(const std::string& filename);
  explicit Configuration(std::istream& fs);
  virtual ~Configuration() = default;

  [[nodiscard]] auto getPatchRows() const -> std::size_t { return patchRows; }
  [[nodiscard]] auto getPatchCols() const -> std::size_t { return patchCols; }
  [[nodiscard]] auto getMethod() const -> NaMappingMethod { return method; }
};
} // namespace na
