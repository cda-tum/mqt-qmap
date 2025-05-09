#include "ir/operations/CompoundOperation.hpp"
#include "ir/operations/NonUnitaryOperation.hpp"
#include "ir/operations/StandardOperation.hpp"
#include "ir/operations/SymbolicOperation.hpp"
#include "na/azac/CodeGenerator.hpp"

#include <cstddef>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>
#include <map>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace na::azac {
constexpr std::string_view architectureJson = R"({
  "name": "code_generator_architecture",
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{"id": 0, "site_separation": [3, 3], "r": 20, "c": 20, "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [60, 60]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {"id": 1, "site_separation": [12, 10], "r": 4, "c": 4, "location": [5, 70]},
      {"id": 2, "site_separation": [12, 10], "r": 4, "c": 4, "location": [7, 70]}
    ],
    "offset": [5, 70],
    "dimension": [50, 40]
  }],
  "aods":[{"id": 0, "site_separation": 2, "r": 20, "c": 20}],
  "arch_range": [[0, 0], [60, 110]],
  "rydberg_range": [[[5, 70], [55, 110]]]
})";
constexpr std::string_view configJson = R"({
  "code_generator" : {
    "parking_offset" : 1,
    "warn_unsupported_gates" : true
  }
})";
class CodeGeneratorGenerateTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  CodeGenerator codeGenerator;
  CodeGeneratorGenerateTest()
      : architecture(nlohmann::json::parse(architectureJson)),
        config(nlohmann::json::parse(configJson)),
        codeGenerator(architecture, config) {}
};
TEST_F(CodeGeneratorGenerateTest, Empty) {
  const auto& slm = *architecture.storageZones.front();
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, GlobalRYGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto ry = qc::StandardOperation(0, qc::RY, {0.1});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {ry}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ ry 0.10000 global\n");
}
TEST_F(CodeGeneratorGenerateTest, GlobalYGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto y = qc::StandardOperation(0, qc::Y);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {y}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ ry 3.14159 global\n");
}
TEST_F(CodeGeneratorGenerateTest, GlobalCompoundRYGate) {
  const auto& slm = *architecture.storageZones.front();
  // Create a compound operation with a single RY gate
  qc::CompoundOperation cry;
  cry.emplace_back<qc::StandardOperation>(0, qc::RY, std::vector{0.1});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {cry}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ ry 0.10000 global\n");
}
TEST_F(CodeGeneratorGenerateTest, GlobalCompoundYGate) {
  const auto& slm = *architecture.storageZones.front();
  qc::CompoundOperation cy;
  cy.emplace_back<qc::StandardOperation>(0, qc::Y);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {cy}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ ry 3.14159 global\n");
}
TEST_F(CodeGeneratorGenerateTest, RZGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto rz = qc::StandardOperation(0, qc::RZ, {0.1});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {rz}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ rz 0.10000 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, PGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto p = qc::StandardOperation(0, qc::P, {0.1});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {p}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ rz 0.10000 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, ZGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto z = qc::StandardOperation(0, qc::Z);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {z}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ rz 3.14159 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, SGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto s = qc::StandardOperation(0, qc::S);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {s}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ rz 1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, SdgGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto sdg = qc::StandardOperation(0, qc::Sdg);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {sdg}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ rz -1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, TGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto t = qc::StandardOperation(0, qc::T);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {t}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ rz 0.78540 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, TdgGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto tdg = qc::StandardOperation(0, qc::Tdg);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {tdg}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ rz -0.78540 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, U3Gate) {
  const auto& slm = *architecture.storageZones.front();
  const auto u = qc::StandardOperation(0, qc::U, {0.1, 0.2, 0.3});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {u}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u 0.10000 0.20000 0.30000 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, U2Gate) {
  const auto& slm = *architecture.storageZones.front();
  const auto u2 = qc::StandardOperation(0, qc::U2, {0.1, 0.2});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {u2}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u 1.57080 0.10000 0.20000 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, RXGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto rx = qc::StandardOperation(0, qc::RX, {0.1});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {rx}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u 0.10000 -1.57080 1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, RYGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto ry = qc::StandardOperation(0, qc::RY, {0.1});
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {ry}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}, {slm, 0, 1}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "atom (3.000, 0.000) atom1\n"
      "@+ u 0.10000 0.00000 0.00000 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, YGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto ry = qc::StandardOperation(0, qc::Y);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {ry}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}, {slm, 0, 1}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "atom (3.000, 0.000) atom1\n"
      "@+ u 3.14159 1.57080 1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, HGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto h = qc::StandardOperation(0, qc::H);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {h}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u 1.57080 0.00000 3.14159 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, XGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto x = qc::StandardOperation(0, qc::X);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {x}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u 3.14159 0.00000 3.14159 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, VGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto v = qc::StandardOperation(0, qc::V);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {v}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u -1.57080 -1.57080 1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, VdgGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto vdg = qc::StandardOperation(0, qc::Vdg);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {vdg}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u -1.57080 1.57080 -1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, SXGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto sx = qc::StandardOperation(0, qc::SX);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {sx}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u 1.57080 -1.57080 1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, SXdgGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto sxdg = qc::StandardOperation(0, qc::SXdg);
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{
                  {sxdg}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{slm, 0, 0}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{})
          .toString(),
      "atom (0.000, 0.000) atom0\n"
      "@+ u -1.57080 -1.57080 1.57080 atom0\n");
}
TEST_F(CodeGeneratorGenerateTest, UnsupportedGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto unsupported = qc::NonUnitaryOperation(0, 0);
  EXPECT_THROW(
      std::ignore = codeGenerator.generate(
          std::vector<std::vector<std::reference_wrapper<const qc::Operation>>>{
              {unsupported}},
          std::vector<std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>{
              {{slm, 0, 0}}},
          std::vector<std::vector<std::vector<qc::Qubit>>>{}),
      std::invalid_argument);
}
TEST_F(CodeGeneratorGenerateTest, TwoQubitGate) {
  const auto& storage = *architecture.storageZones.front();
  const auto& entanglementLeft =
      architecture.entanglementZones.front()->front();
  const auto& entanglementRight =
      architecture.entanglementZones.front()->back();
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{{},
                                                                            {}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{storage, 19, 0}, {storage, 19, 1}},
                  {{entanglementLeft, 0, 0}, {entanglementRight, 0, 0}},
                  {{storage, 19, 0}, {storage, 19, 1}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{{{0U, 1U}},
                                                               {{0U, 1U}}})
          .toString(),
      "atom (0.000, 57.000) atom0\n"
      "atom (3.000, 57.000) atom1\n"
      "@+ load [\n"
      "    atom0\n"
      "    atom1\n"
      "]\n"
      "@+ move [\n"
      "    (5.000, 70.000) atom0\n"
      "    (7.000, 70.000) atom1\n"
      "]\n"
      "@+ store [\n"
      "    atom0\n"
      "    atom1\n"
      "]\n"
      "@+ cz zone_cz0\n"
      "@+ load [\n"
      "    atom0\n"
      "    atom1\n"
      "]\n"
      "@+ move [\n"
      "    (0.000, 57.000) atom0\n"
      "    (3.000, 57.000) atom1\n"
      "]\n"
      "@+ store [\n"
      "    atom0\n"
      "    atom1\n"
      "]\n");
}
TEST_F(CodeGeneratorGenerateTest, Offset) {
  // STORAGE     ...         │ ...         │ ...
  //         18  0 1 o o ... │ o o o o ... │ 0 1 o o ...
  //         19  2 3 o o ... │ o o o o ... │ 2 3 o o ...
  //                         │  ╲╲         │ ↑ ↑
  // ENTANGLEMENT            │   ↓↓        │  ╲╲
  //          0    oo    ... │   01    ... │   oo    ...
  //          1    oo    ... │   23    ... │   oo    ...
  //               ...       │   ...       │   ...
  const auto& storage = *architecture.storageZones.front();
  const auto& entanglementLeft =
      architecture.entanglementZones.front()->front();
  const auto& entanglementRight =
      architecture.entanglementZones.front()->back();
  EXPECT_EQ(
      codeGenerator
          .generate(
              std::vector<
                  std::vector<std::reference_wrapper<const qc::Operation>>>{{},
                                                                            {}},
              std::vector<std::vector<std::tuple<
                  std::reference_wrapper<const SLM>, size_t, size_t>>>{
                  {{storage, 18, 0},
                   {storage, 18, 1},
                   {storage, 19, 0},
                   {storage, 19, 1}},
                  {{entanglementLeft, 0, 0},
                   {entanglementRight, 0, 0},
                   {entanglementLeft, 1, 0},
                   {entanglementRight, 1, 0}},
                  {{storage, 18, 0},
                   {storage, 18, 1},
                   {storage, 19, 0},
                   {storage, 19, 1}}},
              std::vector<std::vector<std::vector<qc::Qubit>>>{
                  {{0U, 1U, 2U, 3U}}, {{0U, 1U, 2U, 3U}}})
          .toString(),
      "atom (0.000, 54.000) atom0\n"
      "atom (0.000, 57.000) atom2\n"
      "atom (3.000, 54.000) atom1\n"
      "atom (3.000, 57.000) atom3\n"
      "@+ load [\n"
      "    atom0\n"
      "    atom1\n"
      "]\n"
      "@+ move [\n"
      "    (0.000, 55.000) atom0\n"
      "    (3.000, 55.000) atom1\n"
      "]\n"
      "@+ load [\n"
      "    atom2\n"
      "    atom3\n"
      "]\n"
      "@+ move [\n"
      "    (5.000, 70.000) atom0\n"
      "    (7.000, 70.000) atom1\n"
      "    (5.000, 80.000) atom2\n"
      "    (7.000, 80.000) atom3\n"
      "]\n"
      "@+ store [\n"
      "    atom0\n"
      "    atom1\n"
      "    atom2\n"
      "    atom3\n"
      "]\n"
      "@+ cz zone_cz0\n"
      "@+ load [\n"
      "    atom0\n"
      "    atom1\n"
      "]\n"
      "@+ move [\n"
      "    (5.000, 71.000) atom0\n"
      "    (7.000, 71.000) atom1\n"
      "]\n"
      "@+ load [\n"
      "    atom2\n"
      "    atom3\n"
      "]\n"
      "@+ move [\n"
      "    (0.000, 54.000) atom0\n"
      "    (3.000, 54.000) atom1\n"
      "    (0.000, 57.000) atom2\n"
      "    (3.000, 57.000) atom3\n"
      "]\n"
      "@+ store [\n"
      "    atom0\n"
      "    atom1\n"
      "    atom2\n"
      "    atom3\n"
      "]\n");
}
TEST(CodeGeneratorTest, InvalidConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({
  "code_generator": {
    "parking_offset": "invalid",
    "warn_unsupported_gates" : "invalid",
    "unknown_key": 42
  }
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = CodeGenerator(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_THAT(
      buffer.str(),
      ::testing::AllOf(
          ::testing::HasSubstr("\033[1;35m[WARN]\033[0m Configuration for "
                               "CodeGenerator contains an invalid value for "
                               "parking_offset. Using default (1)."),
          ::testing::HasSubstr("\033[1;35m[WARN]\033[0m Configuration for "
                               "CodeGenerator contains an invalid value for "
                               "warn_unsupported_gates. Using default (true)."),
          ::testing::HasSubstr("\033[1;35m[WARN]\033[0m Configuration for "
                               "CodeGenerator contains an unknown key: "
                               "unknown_key. Ignoring.")));
  size_t warnings = 0;
  size_t pos = 0;
  std::string target = "\033[1;35m[WARN]\033[0m";
  while ((pos = buffer.str().find(target, pos)) != std::string::npos) {
    ++warnings;
    pos += target.length();
  }
  EXPECT_EQ(warnings, 3);
}
TEST(CodeGeneratorTest, EmptyConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({
  "code_generator": {}
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = CodeGenerator(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_THAT(
      buffer.str(),
      ::testing::AllOf(
          ::testing::HasSubstr("\033[1;35m[WARN]\033[0m Configuration for "
                               "CodeGenerator does not contain a value for "
                               "parking_offset. Using default (1).\n"),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for CodeGenerator "
              "does not contain a value for warn_unsupported_gates. Using "
              "default (true).\n")));
  size_t warnings = 0;
  size_t pos = 0;
  std::string target = "\033[1;35m[WARN]\033[0m";
  while ((pos = buffer.str().find(target, pos)) != std::string::npos) {
    ++warnings;
    pos += target.length();
  }
  EXPECT_EQ(warnings, 2);
}
TEST(CodeGeneratorTest, NoConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = CodeGenerator(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_EQ(buffer.str(),
            "\033[1;35m[WARN]\033[0m Configuration does not contain settings "
            "for CodeGenerator or is malformed. Using default settings.\n");
}
} // namespace na::azac
