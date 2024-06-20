//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

enum class SwapReduction : std::uint8_t {
  None,
  CouplingLimit,
  Custom,
  Increasing
};

[[maybe_unused]] static inline std::string
toString(const SwapReduction strategy) {
  switch (strategy) {
  case SwapReduction::CouplingLimit:
    return "coupling_limit";
  case SwapReduction::Custom:
    return "custom";
  case SwapReduction::None:
    return "none";
  case SwapReduction::Increasing:
    return "increasing";
  }
  return " ";
}

[[maybe_unused]] static SwapReduction
swapReductionFromString(const std::string& reduction) {
  if (reduction == "none" || reduction == "0") {
    return SwapReduction::None;
  }
  if (reduction == "coupling_limit" || reduction == "1") {
    return SwapReduction::CouplingLimit;
  }
  if (reduction == "custom" || reduction == "2") {
    return SwapReduction::Custom;
  }
  if (reduction == "increasing" || reduction == "3") {
    return SwapReduction::Increasing;
  }
  throw std::invalid_argument("Invalid swap reduction value: " + reduction);
}
