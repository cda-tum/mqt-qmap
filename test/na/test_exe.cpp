#include "Architecture.hpp"
#include "Configuration.hpp"
#include "Definitions.hpp"
#include "Layer.hpp"
#include "NeutralAtomMapper.hpp"
#include "QuantumComputation.hpp"
#include "na/Definitions.hpp"
#include "operations/NAGlobalOperation.hpp"
#include "operations/NALocalOperation.hpp"
#include "operations/NAQuantumComputation.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"

#include <iostream>
#include <vector>
namespace na {

auto test() -> void {
  auto qc = qc::QuantumComputation(8);
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
  qc.cz(1, 2);
  qc.cz(1, 6);
  qc.cz(2, 7);
  qc.cz(3, 4);
  qc.cz(4, 7);
  qc.cz(5, 7);
  qc.cz(6, 7);
  qc.emplace_back<qc::StandardOperation>(8, qc::Targets{0, 1, 2, 3, 4, 5, 6, 7},
                                         qc::OpType::RY,
                                         std::vector<qc::fp>{qc::PI_2});
  qc.rz(qc::PI_4, 0);
  qc.rz(qc::PI_4, 1);
  qc.rz(qc::PI_4, 2);
  qc.emplace_back<qc::StandardOperation>(8, qc::Targets{0, 1, 2, 3, 4, 5, 6, 7},
                                         qc::OpType::RY,
                                         std::vector<qc::fp>{-qc::PI_2});
  qc.cz(3, 5);
  qc.cz(4, 5);
  qc.cz(1, 2);
//   qc.dumpOpenQASM2(std::cout);
  // std::cout << "-------------------------------------------------" << std::endl;
  auto layer = Layer(qc);
  auto graph = layer.constructInteractionGraph({qc::OpType::Z, 1});
  auto arch = Architecture("examples/na/nature.json", "examples/na/nature.csv");
  auto mapper = NeutralAtomMapper(arch, Configuration{3, 2});
  mapper.map(qc);
  const auto& mappedQc = mapper.getResult();
  std::cout << mappedQc << std::endl;
}

} // namespace na

auto main() -> int {
  na::test();
  return 0;
}
