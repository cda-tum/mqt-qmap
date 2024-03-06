//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/mqt-qmap for more
// information.
//

#pragma once

#include <string>

namespace na {
class Configuration {
private:
  std::size_t  patchRows             = 1;
  std::size_t  patchCols             = 1;

public:
  Configuration() = default;
  explicit Configuration(std::size_t rows, std::size_t cols)
      : patchRows(rows), patchCols(cols){};
  explicit Configuration(const std::string& filename);
  explicit Configuration(std::istream& fs);
  virtual ~Configuration() = default;

  [[nodiscard]] auto getPatchRows() const -> std::size_t { return patchRows; }
  [[nodiscard]] auto getPatchCols() const -> std::size_t { return patchCols; }
};
} // namespace na
