//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#ifndef QMAP_SWAPREDUCTION_HPP
#define QMAP_SWAPREDUCTION_HPP

#include <iostream>

enum class SwapReduction { None, CouplingLimit, Custom, Increasing };

static std::string toString(const SwapReduction strategy) {
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
  } else if (reduction == "coupling_limit" || reduction == "1") {
    return SwapReduction::CouplingLimit;
  } else if (reduction == "custom" || reduction == "2") {
    return SwapReduction::Custom;
  } else if (reduction == "increasing" || reduction == "3") {
    return SwapReduction::Increasing;
  } else {
    throw std::invalid_argument("Invalid swap reduction value: " + reduction);
  }
}

#endif // QMAP_SWAPREDUCTION_HPP
