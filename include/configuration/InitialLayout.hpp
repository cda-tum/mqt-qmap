//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#ifndef QMAP_INITIALLAYOUT_HPP
#define QMAP_INITIALLAYOUT_HPP

#include <iostream>

/// Identity: q_i -> Q_i
/// Static: first layer is mapped q_c -> Q_c and q_t -> Q_t
/// Dynamic: Layout is generated on demand upon encountering a specific gate
enum class InitialLayout { None, Identity, Static, Dynamic };

static std::string toString(const InitialLayout strategy) {
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
  } else if (initialLayout == "identity" || initialLayout == "1") {
    return InitialLayout::Identity;
  } else if (initialLayout == "static" || initialLayout == "2") {
    return InitialLayout::Static;
  } else if (initialLayout == "dynamic" || initialLayout == "3") {
    return InitialLayout::Dynamic;
  } else {
    throw std::invalid_argument("Invalid initial layout value: " +
                                initialLayout);
  }
}

#endif // QMAP_INITIALLAYOUT_HPP
