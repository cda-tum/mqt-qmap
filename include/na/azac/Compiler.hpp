#pragma once

#include "na/azac/CompilerBase.hpp"
#include "na/azac/Placer.hpp"
#include "na/azac/Router.hpp"
#include "na/azac/Scheduler.hpp"

namespace na {

class Compiler : public CompilerBase,
                 private Scheduler<Compiler>,
                 private Placer<Compiler>,
                 private Router<Compiler> {

};

} // namespace na
