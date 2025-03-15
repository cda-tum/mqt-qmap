#include "ir/operations/StandardOperation.hpp"
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

namespace na {
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
    "parking_offset" : 1
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
TEST_F(CodeGeneratorGenerateTest, OneQubitGate) {
  const auto& slm = *architecture.storageZones.front();
  const auto rz = qc::StandardOperation(0, qc::Z, {qc::PI});
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
      "@+ rz 3.14159 atom0\n");
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
          ::testing::MatchesRegex(
              ".*\\[WARN\\].*\n.*\\[WARN\\].*\n.*\\[WARN\\].*\n"),
          ::testing::HasSubstr("\033[1;35m[WARN]\033[0m Configuration for "
                               "CodeGenerator contains an invalid "
                               "value for parking_offset. Using default."),
          ::testing::HasSubstr(
              "\033[1;35m[WARN]\033[0m Configuration for CodeGenerator does "
              "not contain a value for parking_offset. Using default."),
          ::testing::HasSubstr("\033[1;35m[WARN]\033[0m Configuration for "
                               "CodeGenerator contains an "
                               "unknown key: unknown_key. Ignoring.")));
}
} // namespace na
