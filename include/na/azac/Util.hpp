#pragma once

#include <cmath>
#include <utility>

namespace na {

template <typename T1, typename T2>
auto distance(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b)
    -> double {
  return std::sqrt(
      std::pow(static_cast<double>(a.first) - static_cast<double>(b.first), 2) +
      std::pow(static_cast<double>(a.second) - static_cast<double>(b.second),
               2));
}

} // namespace na