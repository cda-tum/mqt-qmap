//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

enum class Method : std::uint8_t { None, Exact, Heuristic };

[[maybe_unused]] static inline std::string toString(const Method method) {
  switch (method) {
  case Method::None:
    return "none";
  case Method::Exact:
    return "exact";
  case Method::Heuristic:
    return "heuristic";
  }
  return " ";
}

[[maybe_unused]] static Method methodFromString(const std::string& method) {
  if (method == "none" || method == "0") {
    return Method::None;
  }
  if (method == "exact" || method == "1") {
    return Method::Exact;
  }
  if (method == "heuristic" || method == "2") {
    return Method::Heuristic;
  }
  throw std::invalid_argument("Invalid method value: " + method);
}
