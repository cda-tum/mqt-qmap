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
auto solve(bool save_file = true) -> Result {
    // member to hold intermediate results
    gate_scheduling.clear();
    gate_1q_scheduling.clear();
    qubit_mapping.clear();
    reuse_qubits.clear();
    qubit_mapping.clear();

    std::cout << "[INFO] ZAC: A compiler for neutral atom-based compute-store "
                 "architecture\n";
    std::cout << *this;
    // todo: check if the program input is valid, i.e., #q < #p
    const auto tS = std::chrono::system_clock::now();
    // gate scheduling with graph coloring
    std::cout << "[INFO] ZAC: Run scheduling\n";
    static_cast<T*>(this)->schedule();

    if (reuse) {
      collect_reuse_qubit();
    } else {
      reuse_qubits.reserve(gate_scheduling.size());
      for (std::size_t i = 0; i < gate_scheduling.size(); ++i) {
        reuse_qubits.emplace_back();
      }
    }

    static_cast<T*>(this)->place_qubit_initial();
    std::cout << "[INFO]               Time for initial placement: "
              << runtime_analysis.initial_placement.count() << "s\n";
    static_cast<T*>(this)->place_qubit_intermediate();
    std::cout << "[INFO]               Time for intermediate placement: ";
    std::cout << runtime_analysis.intermediate_placement.count() << "s\n";
    static_cast<T*>(this)->route_qubit();
    runtime_analysis.total = std::chrono::system_clock::now() - tS;
    std::cout << "[INFO]               Time for routing: "
              << runtime_analysis.routing.count() << "s\n";
    std::cout << "[INFO] ZAC: Toal Time: " << runtime_analysis.total.count() << "s\n";
    if (save_file) {
      if (dir.empty()) {
        dir = "./result/";
      }
      const auto code_filename = dir / "code" / (result.name + "_code.na");
      std::ofstream code_ofs(code_filename);
      code_ofs << result.instructions;

      const auto result_json_filename =
          dir / "time" / (result.name + "_time.json");
      std::ofstream result_json_ofs(result_json_filename);
      nlohmann::json result_json{};
      result_json["scheduling"] = runtime_analysis.scheduling.count();
      result_json["initial_placement"] = runtime_analysis.initial_placement.count();
      result_json["intermediate_placement"] =
          runtime_analysis.intermediate_placement.count();
      result_json["routing"] = runtime_analysis.routing.count();
      result_json["total"] = runtime_analysis.total.count();
      result_json_ofs << result_json;
    }

    if (to_verify) {
      std::cout << "[INFO] ZAC: Start Verification\n";
      throw std::invalid_argument("Verification is not implemented yet");
    }

    return result;
  }
};

class ZACompiler final : public Compiler<ZACompiler, Placer, Router, Scheduler> {};

} // namespace na
