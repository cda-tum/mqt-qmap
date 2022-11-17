/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "Architecture.hpp"
#include "QuantumComputation.hpp"
#include "cliffordsynthesis/Tableau.hpp"

#include "gtest/gtest.h"

class TestTableau : public testing::TestWithParam<std::string> {
protected:
  Tableau tableau{"['+ZI', '+IZ']"};
  Tableau tableau2{};
};

TEST_F(TestTableau, InitialTableau) {
  tableau2.init(2);

  EXPECT_EQ(tableau, tableau2);

  EXPECT_EQ(tableau[0][0], 0);
  EXPECT_EQ(tableau[0][1], 0);
  EXPECT_EQ(tableau[0][2], 1);
  EXPECT_EQ(tableau[0][3], 0);
  EXPECT_EQ(tableau[0][4], 0);
  EXPECT_EQ(tableau[1][0], 0);
  EXPECT_EQ(tableau[1][1], 0);
  EXPECT_EQ(tableau[1][2], 0);
  EXPECT_EQ(tableau[1][3], 1);
  EXPECT_EQ(tableau[1][4], 0);

  const std::string result_string = "0;0;1;0;0;\n0;0;0;1;0;\n";
  EXPECT_EQ(tableau.toString(), result_string);
}

TEST_F(TestTableau, TableauIO) {
  const auto filename = "tableau.txt";
  tableau.dump(filename);
  tableau2.import(filename);
  EXPECT_EQ(tableau, tableau2);
}

TEST_F(TestTableau, BellCircuit) {
  using namespace dd::literals;

  auto qc = qc::QuantumComputation(2U);
  qc.h(0);
  qc.x(1, 0_pc);

  tableau = Tableau(qc);

  tableau2.fromString("[+XX, +ZZ]");
  EXPECT_EQ(tableau, tableau2);
}

TEST_F(TestTableau, TestOperations) {
  using namespace dd::literals;

  qc::QuantumComputation qc1(2U);
  qc1.x(0);
  qc1.y(0);
  qc1.z(0);
  qc1.h(0);
  qc1.s(0);
  qc1.sdag(0);
  qc1.x(1, 0_pc);
  qc1.y(1, 0_pc);
  qc1.z(1, 0_pc);
  qc1.swap(0, 1);

  EXPECT_NO_THROW(tableau = Tableau(qc1););
}

TEST_F(TestTableau, TestCompoundOperation) {
  using namespace dd::literals;

  auto qc = qc::QuantumComputation(2U);

  auto compOP = std::make_unique<qc::CompoundOperation>(2);
  compOP->emplace_back<qc::StandardOperation>(2, 0, qc::H);
  compOP->emplace_back<qc::StandardOperation>(2, 0_pc, 1, qc::X);
  qc.emplace_back(compOP);

  tableau = Tableau(qc);

  tableau2.fromString("[+XX, +ZZ]");
  EXPECT_EQ(tableau, tableau2);
}

TEST_F(TestTableau, BVTableau) {
  const auto bitvector1 = 0b10;
  const auto bitvector2 = 0b01;

  tableau.populateTableauFrom(bitvector1, 2, 0);
  tableau.populateTableauFrom(bitvector2, 2, 1);

  EXPECT_EQ(tableau[0][0], 0);
  EXPECT_EQ(tableau[1][0], 1);
  EXPECT_EQ(tableau[0][1], 1);
  EXPECT_EQ(tableau[1][1], 0);

  const auto bitvector3 = tableau.getBVFrom(0);
  const auto bitvector4 = tableau.getBVFrom(1);

  EXPECT_EQ(bitvector3, bitvector1);
  EXPECT_EQ(bitvector4, bitvector2);
}
