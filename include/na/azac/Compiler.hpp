#pragma once

#include "na/azac/CompilerBase.hpp"
#include "na/azac/Placer.hpp"
#include "na/azac/Router.hpp"
#include "na/azac/Scheduler.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <nlohmann/json_fwd.hpp>
#include <stdexcept>

namespace na {

template <typename T, template <typename> class... Mixins>
class Compiler : public CompilerBase, public Mixins<T>... {
private:
  Compiler() = default;
  friend T;

public:
  auto solve(bool saveFile = true) -> const Result& {
    // member to hold intermediate results
    gateScheduling.clear();
    gate1QScheduling.clear();
    qubitMapping.clear();
    reuseQubits.clear();
    qubitMapping.clear();

    std::cout << "[INFO] AZAC: An advanced compiler for neutral atom-based "
                 "compute-store "
                 "architecture\n";
    std::cout << *this;
    // todo: check if the program input is valid, i.e., #q < #p
    const auto tS = std::chrono::system_clock::now();
    // gate scheduling with graph coloring
    std::cout << "[INFO] AZAC: Run scheduling\n";
    static_cast<T*>(this)->schedule();

    if (reuse) {
      collectReuseQubit();
    } else {
      reuseQubits.reserve(gateScheduling.size());
      for (std::size_t i = 0; i < gateScheduling.size(); ++i) {
        reuseQubits.emplace_back();
      }
    }

    static_cast<T*>(this)->placeQubitInitial();
    std::cout << "[INFO]           Time for initial placement: "
              << runtimeAnalysis.initialPlacement.count() << "µs\n";
    static_cast<T*>(this)->placeQubitIntermediate();
    std::cout << "[INFO]           Time for intermediate placement: ";
    std::cout << runtimeAnalysis.intermediatePlacement.count() << "µs\n";
    static_cast<T*>(this)->routeQubit();
    runtimeAnalysis.total = std::chrono::system_clock::now() - tS;
    std::cout << "[INFO]           Time for routing: "
              << runtimeAnalysis.routing.count() << "µs\n";
    std::cout << "[INFO] AZAC: Total Time: " << runtimeAnalysis.total.count()
              << "µs\n";
    if (saveFile) {
      if (dir.empty()) {
        dir = "./result/";
      }
      const auto codeFilename = dir / "code" / (result.name + "_code.json");
      create_directories(codeFilename.parent_path());
      std::ofstream codeOfs(codeFilename);
      if (!codeOfs) {
        std::stringstream ss{};
        ss << "Cannot open file " << absolute(codeFilename);
        throw std::runtime_error(ss.str());
      }
      nlohmann::json resultJson;
      resultJson["instructions"] = result.instructions;
      resultJson["runtime"] = result.runtime;
      resultJson["name"] = result.name;
      resultJson["architecture_spec_path"] = result.architectureSpecPath;
      codeOfs << resultJson;
      std::cout << "[INFO]           Saved code to " << absolute(codeFilename)
                << "\n";

      const auto timingJsonFilename =
          dir / "time" / (result.name + "_time.json");
      create_directories(timingJsonFilename.parent_path());
      std::ofstream timingtJsonOfs(timingJsonFilename);
      if (!timingtJsonOfs) {
        std::stringstream ss{};
        ss << "Cannot open file " << absolute(timingJsonFilename);
        throw std::runtime_error(ss.str());
      }
      nlohmann::json timingJson{};
      timingJson["scheduling"] = runtimeAnalysis.scheduling.count();
      timingJson["initial_placement"] =
          runtimeAnalysis.initialPlacement.count();
      timingJson["intermediate_placement"] =
          runtimeAnalysis.intermediatePlacement.count();
      timingJson["routing"] = runtimeAnalysis.routing.count();
      timingJson["total"] = runtimeAnalysis.total.count();
      timingtJsonOfs << timingJson;
      std::cout << "[INFO]           Saved results to "
                << absolute(timingJsonFilename) << "\n";

      //===----------------------------------------------------------------===//
      // NAComputation: save NAComputation to file
      const auto naFilename = dir / "na" / (result.name + "_code.naviz");
      create_directories(naFilename.parent_path());
      std::ofstream naOfs(naFilename);
      if (!naOfs) {
        std::stringstream ss{};
        ss << "Cannot open file " << absolute(naFilename);
        throw std::runtime_error(ss.str());
      }
      naOfs << result.naComputation;
      std::cout << "[INFO]           Saved NAComputation to "
                << absolute(naFilename) << "\n";
      //===----------------------------------------------------------------===//
    }

    if (toVerify) {
      std::cout << "[INFO] AZAC: Start Verification\n";
      throw std::invalid_argument("Verification is not implemented yet");
    }

    return result;
  }
};

class AZACompiler final
    : public Compiler<AZACompiler, Placer, Router, Scheduler> {};

} // namespace na
