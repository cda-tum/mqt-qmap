//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <iostream>

/// Identity: q_i -> Q_i
/// Static: first layer is mapped q_c -> Q_c and q_t -> Q_t
/// Dynamic: Layout is generated on demand upon encountering a specific gate
enum class InitialLayout { None, Identity, Static, Dynamic };

[[maybe_unused]] static inline std::string
toString(const InitialLayout strategy) {
  switch (strategy) {
  case InitialLayout::Identity:
    return "identity";
  case InitialLayout::Static:
    return "static";
  case InitialLayout::Dynamic:
    return "dynamic";
  case InitialLayout::None:
    return "none";
  }
  return " ";
}

[[maybe_unused]] static InitialLayout
initialLayoutFromString(const std::string& initialLayout) {
  if (initialLayout == "none" || initialLayout == "0") {
    return InitialLayout::None;
  }
  if (initialLayout == "identity" || initialLayout == "1") {
    return InitialLayout::Identity;
  }
  if (initialLayout == "static" || initialLayout == "2") {
    return InitialLayout::Static;
  }
  if (initialLayout == "dynamic" || initialLayout == "3") {
    return InitialLayout::Dynamic;
  }
  throw std::invalid_argument("Invalid initial layout value: " + initialLayout);
}
