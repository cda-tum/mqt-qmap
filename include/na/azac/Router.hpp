#pragma once

#include "na/azac/CompilerBase.hpp"

#include <type_traits>

namespace na {

template <typename T> class Router {
  static_assert(std::is_base_of_v<CompilerBase, T>,
                "T must be a subclass of CompilerBase");
};

} // namespace na
