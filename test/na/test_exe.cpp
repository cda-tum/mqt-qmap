#include "Architecture.hpp"
#include "Configuration.hpp"
#include "Definitions.hpp"
#include "Layer.hpp"
#include "NeutralAtomMapper.hpp"
#include "QuantumComputation.hpp"
#include "na/Definitions.hpp"
#include "operations/CompoundOperation.hpp"
#include "operations/NAGlobalOperation.hpp"
#include "operations/NALocalOperation.hpp"
#include "operations/NAQuantumComputation.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <ostream>
#include <string>
#include <vector>

namespace na {

auto test(const std::string& inputFile, const std::string& architecture,
          const std::string& layout, const std::string& configuration,
          std::ostream& outputStream = std::cout) -> void {
  /*
                     ┌─────────┐┌─────────┐┌──────────┐
q_0: ────────────────┤         ├┤ Rz(π/4) ├┤          ├──────
                     |         |├─────────┤|          |
q_1: ──■──■──────────┤         ├┤ Rz(π/4) ├┤          ├─■────
       │  │          |         |├─────────┤|          | │
q_2: ──■──┼──■───────┤         ├┤ Rz(π/4) ├┤          ├─■────
          │  │       |         |└─────────┘|          |
q_3: ──■──┼──┼───────┤ Ry(π/2) ├───────────┤ Ry(-π/2) ├─■────
       │  │  │       |         |           |          | │
q_4: ──■──┼──┼──■────┤         ├───────────┤          ├─┼──■─
          │  │  │    |         |           |          | │  │
q_5: ─────┼──┼──┼──■─┤         ├───────────┤          ├─■──■─
          │  │  │  │ |         |           |          |
q_6: ──■──■──┼──┼──┼─┤         ├───────────┤          ├──────
       │     │  │  │ |         |           |          |
q_7: ──■─────■──■──■─┤         ├───────────┤          ├──────
                     └─────────┘           └──────────┘
  */
  //   const std::vector<qc::Qubit> qubits{0, 1, 2, 3, 4, 5, 6, 7};
  //   auto                         qc = qc::QuantumComputation(8);
  //   qc.cz(1, 2);
  //   qc.cz(1, 6);
  //   qc.cz(2, 7);
  //   qc.cz(3, 4);
  //   qc.cz(4, 7);
  //   qc.cz(5, 7);
  //   qc.cz(6, 7);
  //   qc.emplace_back<qc::CompoundOperation>(
  //       std::accumulate(qubits.begin(), qubits.end(),
  //       qc::CompoundOperation(8),
  //                       [](qc::CompoundOperation& op, const auto q) {
  //                         op.emplace_back<qc::StandardOperation>(
  //                             8, q, qc::RY, std::vector<qc::fp>{qc::PI_2});
  //                         return op;
  //                       }));
  //   qc.rz(qc::PI_4, 0);
  //   qc.rz(qc::PI_4, 1);
  //   qc.rz(qc::PI_4, 2);
  //   qc.emplace_back<qc::CompoundOperation>(
  //       std::accumulate(qubits.begin(), qubits.end(),
  //       qc::CompoundOperation(8),
  //                       [](qc::CompoundOperation& op, const auto q) {
  //                         op.emplace_back<qc::StandardOperation>(
  //                             8, q, qc::RY, std::vector<qc::fp>{-qc::PI_2});
  //                         return op;
  //                       }));
  //   qc.cz(3, 5);
  //   qc.cz(4, 5);
  //   qc.cz(1, 2);
  const qc::QuantumComputation qc(inputFile);
  //   qc.dumpOpenQASM2(outputStream);
  auto layer  = Layer(qc);
  auto graph  = layer.constructInteractionGraph({qc::OpType::Z, 1});
  auto arch   = Architecture(architecture, layout);
  auto mapper = NeutralAtomMapper(arch, Configuration(configuration));
  mapper.map(qc);
  const auto& mappedQc = mapper.getResult();
  outputStream << mappedQc << std::endl;
}

} // namespace na

struct Options {
  std::string architecture;
  std::string layout;
  std::string configuration;
  std::string inputFilename;
  std::string outputFilename;
};

Options parseCommandLine(int argc, char* argv[]) {
  Options options;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];

    if (arg == "-a" || arg == "--architecture") {
      if (i + 1 < argc) {
        options.architecture = argv[++i];
      } else {
        std::cerr << "Error: Missing argument for architecture option"
                  << std::endl;
        exit(1);
      }
    } else if (arg == "-l" || arg == "--layout") {
      if (i + 1 < argc) {
        options.layout = argv[++i];
      } else {
        std::cerr << "Error: Missing argument for layout option" << std::endl;
        exit(1);
      }
    } else if (arg == "-c" || arg == "--configuration") {
      if (i + 1 < argc) {
        options.configuration = argv[++i];
      } else {
        std::cerr << "Error: Missing argument for configuration option"
                  << std::endl;
        exit(1);
      }
    } else if (arg == "-i" || arg == "--input") {
      if (i + 1 < argc) {
        options.inputFilename = argv[++i];
      } else {
        std::cerr << "Error: Missing argument for input filename" << std::endl;
        exit(1);
      }
    } else if (options.inputFilename.empty()) {
      options.inputFilename = argv[i];
    } else if (arg == "-o" || arg == "--output") {
      if (i + 1 < argc) {
        options.outputFilename = argv[++i];
      } else {
        std::cerr << "Error: Missing argument for output filename" << std::endl;
        exit(1);
      }
    } else if (options.outputFilename.empty()) {
      options.outputFilename = argv[i];
    } else {
      std::cerr << "Error: Unknown option or argument: " << arg << std::endl;
      exit(1);
    }
  }
  if (options.inputFilename.empty()) {
    std::cerr << "Error: Missing input filename" << std::endl;
    exit(1);
  }
  return options;
}

auto main(int argc, char* argv[]) -> int {
  Options options = parseCommandLine(argc, argv);
  if (options.architecture.empty()) {
    options.architecture = "examples/na/nature.json";
  }
  if (options.layout.empty()) {
    options.layout = "examples/na/nature.csv";
  }
  if (options.configuration.empty()) {
    options.configuration = "examples/na/config.json";
  }
  if (!options.outputFilename.empty()) {
    std::ofstream outputStream(options.outputFilename);
    if (!outputStream) {
      std::cerr << "Error: Could not open output file: "
                << options.outputFilename << std::endl;
      return 1;
    }
    na::test(options.inputFilename, options.architecture, options.layout,
             options.configuration, outputStream);
  } else {
    na::test(options.inputFilename, options.architecture, options.layout,
             options.configuration);
  }
}
