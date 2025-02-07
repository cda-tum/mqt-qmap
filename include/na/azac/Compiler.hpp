#pragma once

#include "na/azac/CompilerBase.hpp"
#include "na/azac/Placer.hpp"
#include "na/azac/Router.hpp"
#include "na/azac/Scheduler.hpp"
#include <iostream>
#include <chrono>
#include <stdexcept>
#include <fstream>
#include <nlohmann/json_fwd.hpp>

namespace na {

template <typename T, template <typename> class... Mixins>
class Compiler : public CompilerBase, public Mixins<T>... {
private:
  Compiler() = default;
  friend T;
public:
auto solve(bool saveFile = true) -> Result {
    // member to hold intermediate results
    gateScheduling.clear();
    gate1QScheduling.clear();
    qubitMapping.clear();
    reuseQubits.clear();
    qubitMapping.clear();

    std::cout << "[INFO] ZAC: A compiler for neutral atom-based compute-store "
                 "architecture\n";
    std::cout << *this;
    // todo: check if the program input is valid, i.e., #q < #p
    const auto tS = std::chrono::system_clock::now();
    // gate scheduling with graph coloring
    std::cout << "[INFO] ZAC: Run scheduling\n";
    static_cast<T*>(this)->schedule();

    if (reuse) {
      collectReuseQubit();
    } else {
      reuseQubits.reserve (gateScheduling.size());
      for (std::size_t i = 0; i < gateScheduling.size(); ++i) {
        reuseQubits.emplace_back();
      }
    }

    static_cast<T*>(this)->placeQubitInitial();
    std::cout << "[INFO]               Time for initial placement: "
              << runtimeAnalysis.initialPlacement.count() << "s\n";
    static_cast<T*>(this)->placeQubitIntermediate();
    std::cout << "[INFO]               Time for intermediate placement: ";
    std::cout << runtimeAnalysis.intermediatePlacement.count() << "s\n";
    static_cast<T*>(this)->routeQubit();
    runtimeAnalysis.total = std::chrono::system_clock::now() - tS;
    std::cout << "[INFO]               Time for routing: "
              << runtimeAnalysis.routing.count() << "s\n";
    std::cout << "[INFO] ZAC: Toal Time: " << runtimeAnalysis.total.count() << "s\n";
    if  (saveFile) {
      if (dir.empty()) {
        dir = "./result/";
      }
      const auto codeFilename = dir / "code" / (result.name + "_code.na");
      std::ofstream codeOfs (codeFilename);
      codeOfs << result.instructions;

      const auto resultJsonFilename =
          dir / "time" / (result.name + "_time.json");
      std::ofstream resultJsonOfs (resultJsonFilename);
      nlohmann::json resultJson{};
      resultJson["scheduling"] = runtimeAnalysis.scheduling.count();
      resultJson["initial_placement"] = runtimeAnalysis.initialPlacement.count();
      resultJson["intermediate_placement"] =
          runtimeAnalysis.intermediatePlacement.count();
      resultJson["routing"] = runtimeAnalysis.routing.count();
      resultJson["total"] = runtimeAnalysis.total.count();
      resultJsonOfs << resultJson;
    }

    if  (toVerify) {
      std::cout << "[INFO] ZAC: Start Verification\n";
      throw std::invalid_argument("Verification is not implemented yet");
    }

    return result;
  }
};

class ZACompiler final : public Compiler<ZACompiler, Placer, Router, Scheduler> {};

} // namespace na
