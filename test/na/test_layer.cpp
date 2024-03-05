//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Definitions.hpp"
#include "Layer.hpp"
#include "QuantumComputation.hpp"
#include "operations/OpType.hpp"
#include "operations/StandardOperation.hpp"

#include "gtest/gtest.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

TEST(TestNALayer, ExecutableSet) {
  auto qc = qc::QuantumComputation(3);
  /* construct the following circuit
┌─────────┐┌─────────┐┌──────────┐      ┌─────────┐┌─────────┐┌──────────┐
┤         ├┤ Rz(π/4) ├┤          ├─■──■─┤         ├┤ Rz(π/4) ├┤          ├───
│         │├─────────┤│          │ │  | │         │└─────────┘│          │
┤ Ry(π/2) ├┤ Rz(π/4) ├┤ Ry(-π/2) ├─■──┼─┤ Ry(π/2) ├───────────┤ Ry(-π/2) ├─■─
│         │├─────────┤│          │    │ │         │           │          │ │
┤         ├┤ Rz(π/4) ├┤          ├────■─┤         ├───────────┤          ├─■─
└─────────┘└─────────┘└──────────┘      └─────────┘           └──────────┘
    (1)        (2)        (3)     (4)(5)    (6)        (7)        (8)     (9)
  */
  qc.emplace_back<qc::StandardOperation>(
      3, qc::Targets{0, 1, 2}, qc::OpType::RY, std::vector<qc::fp>{qc::PI_2});
  qc.rz(qc::PI_4, 0);
  qc.rz(qc::PI_4, 1);
  qc.rz(qc::PI_4, 2);
  qc.emplace_back<qc::StandardOperation>(
      3, qc::Targets{0, 1, 2}, qc::OpType::RY, std::vector<qc::fp>{-qc::PI_2});
  qc.cz(0, 1);
  qc.cz(0, 2);
  qc.emplace_back<qc::StandardOperation>(
      3, qc::Targets{0, 1, 2}, qc::OpType::RY, std::vector<qc::fp>{qc::PI_2});
  qc.rz(qc::PI_4, 0);
  qc.emplace_back<qc::StandardOperation>(
      3, qc::Targets{0, 1, 2}, qc::OpType::RY, std::vector<qc::fp>{-qc::PI_2});
  qc.cz(1, 2);

  na::Layer const layer(qc);
  EXPECT_EQ((*layer.getExecutableSet())->size(), 1); // layer (1)
  std::shared_ptr<na::Layer::DAGVertex> v =
      *(*layer.getExecutableSet())->begin();
  na::Layer::execute(v);
  EXPECT_EQ((*layer.getExecutableSet())->size(), 3); // layer (2)
  v = *(*layer.getExecutableSet())->begin();
  na::Layer::execute(v);
  EXPECT_EQ((*layer.getExecutableSet())->size(), 2); // rest of layer (2)
  v = *(*layer.getExecutableSet())->begin();
  na::Layer::execute(v);
  EXPECT_EQ((*layer.getExecutableSet())->size(), 1); // rest of layer (2)
  v = *(*layer.getExecutableSet())->begin();
  na::Layer::execute(v);
  EXPECT_EQ((*layer.getExecutableSet())->size(), 1); // layer (3)
  v = *(*layer.getExecutableSet())->begin();
  na::Layer::execute(v);
  EXPECT_EQ((*layer.getExecutableSet())->size(), 3); // layer (4), (5), (9)
  // execute layer (4) and (5)
  for (const auto& u : **layer.getExecutableSet()) {
    if (const auto& it = (*u->getOperation())->getUsedQubits();
        it.find(0) != it.end()) {
      na::Layer::execute(u);
    }
  }
  EXPECT_EQ((*layer.getExecutableSet())->size(), 2); // layer (6), (9)
}

TEST(TestNALayer, AllExecutable) {
  qc::QuantumComputation qc{};
  na::Layer              layer{};
  qc = qc::QuantumComputation(8);
  qc.cz(1, 2);
  qc.cz(1, 6);
  qc.cz(2, 7);
  qc.cz(3, 4);
  qc.cz(3, 5);
  qc.cz(4, 5);
  qc.cz(4, 6);
  qc.cz(4, 7);
  qc.cz(5, 7);
  qc.cz(6, 7);
  layer = na::Layer(qc);
  EXPECT_EQ((*layer.getExecutableSet())->size(), 10);
}
