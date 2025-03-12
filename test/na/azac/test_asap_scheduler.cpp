#include "ir/operations/StandardOperation.hpp"
#include "na/azac/ASAPScheduler.hpp"

#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>
#include <gtest/gtest.h>
#include <utility>
#include <vector>

namespace testing {
MATCHER_P(RefEq, value, "") { return arg.get() == value; }
} // namespace testing
namespace na {
constexpr std::string_view architectureJson = R"({
  "name": "asap_scheduler_architecture",
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
class ASAPSchedulerScheduleTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  ASAPScheduler scheduler;
  ASAPSchedulerScheduleTest()
      : architecture(nlohmann::json::parse(architectureJson)),
        scheduler(architecture, config) {}
};
TEST_F(ASAPSchedulerScheduleTest, NoGate) {
  qc::QuantumComputation qc;
  const auto& [oneQubitGateLayers, twoQubitGateLayers] = scheduler.schedule(qc);
  EXPECT_THAT(oneQubitGateLayers, ::testing::IsEmpty());
  EXPECT_THAT(twoQubitGateLayers, ::testing::IsEmpty());
}
TEST_F(ASAPSchedulerScheduleTest, OneQubitGate) {
  //    ┌───────┐
  // q: ┤ Rz(π) ├
  //    └───────┘
  qc::QuantumComputation qc(1);
  qc.rz(qc::PI, 0);
  const auto& [oneQubitGateLayers, twoQubitGateLayers] = scheduler.schedule(qc);
  EXPECT_THAT(oneQubitGateLayers,
              ::testing::ElementsAre(::testing::ElementsAre(::testing::RefEq(
                  static_cast<qc::StandardOperation&>(*qc.at(0))))));
  EXPECT_THAT(twoQubitGateLayers, ::testing::IsEmpty());
}
TEST_F(ASAPSchedulerScheduleTest, TwoQubitGate) {
  // q_0: ─■─
  //       │
  // q_1: ─■─
  qc::QuantumComputation qc(2);
  qc.cz(0, 1);
  const auto& [oneQubitGateLayers, twoQubitGateLayers] = scheduler.schedule(qc);
  EXPECT_THAT(oneQubitGateLayers, ::testing::ElementsAre(::testing::IsEmpty(),
                                                         ::testing::IsEmpty()));
  EXPECT_THAT(twoQubitGateLayers, ::testing::ElementsAre(::testing::ElementsAre(
                                      ::testing::Pair(0U, 1U))));
}
TEST_F(ASAPSchedulerScheduleTest, OneQubitSandwich) {
  // q_0: ──────────■──────────
  //      ┌───────┐ │ ┌───────┐
  // q_1: ┤ Rz(π) ├─■─┤ Rz(π) ├
  //      └───────┘   └───────┘
  qc::QuantumComputation qc(2);
  qc.rz(qc::PI, 1);
  qc.cz(0, 1);
  qc.rz(qc::PI, 1);
  const auto& [oneQubitGateLayers, twoQubitGateLayers] = scheduler.schedule(qc);
  EXPECT_THAT(oneQubitGateLayers,
              ::testing::ElementsAre(
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(0)))),
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(2))))));
  EXPECT_THAT(twoQubitGateLayers, ::testing::ElementsAre(::testing::ElementsAre(
                                      ::testing::Pair(0U, 1U))));
}
TEST_F(ASAPSchedulerScheduleTest, TwoQubitSequence) {
  // q_0: ─■───────
  //       │
  // q_1: ─■──■────
  //          │
  // q_2: ────■──■─
  //             │
  // q_3: ───────■─
  qc::QuantumComputation qc(4);
  qc.cz(0, 1);
  qc.cz(1, 2);
  qc.cz(2, 3);
  const auto& [oneQubitGateLayers, twoQubitGateLayers] = scheduler.schedule(qc);
  EXPECT_THAT(oneQubitGateLayers, ::testing::SizeIs(4));
  EXPECT_THAT(oneQubitGateLayers, ::testing::Each(::testing::IsEmpty()));
  EXPECT_THAT(
      twoQubitGateLayers,
      ::testing::ElementsAre(::testing::ElementsAre(::testing::Pair(0U, 1U)),
                             ::testing::ElementsAre(::testing::Pair(1U, 2U)),
                             ::testing::ElementsAre(::testing::Pair(2U, 3U))));
}
TEST_F(ASAPSchedulerScheduleTest, Mixed) {
  //            INPUT ORDER                         SCHEDULED ORDER
  // q_0: ─■─────────────────────────  >>>  ─────────░─■─░─────────░───░─
  //       │ ┌───────┐                 >>>           ░ │ ░┌───────┐░   ░
  // q_1: ─■─┤ Rz(π) ├─────────────■─  >>>  ─────────░─■─░┤ Rz(π) ├░─■─░─
  //         └───────┘┌───────┐    │   >>>  ┌───────┐░   ░└───────┘░ │ ░
  // q_2: ────────────┤ Rz(π) ├─■──■─  >>>  ┤ Rz(π) ├░─■─░─────────░─■─░─
  //                  └───────┘ │      >>>  └───────┘░ │ ░         ░   ░
  // q_3: ──────────────────────■────  >>>  ─────────░─■─░─────────░───░─
  qc::QuantumComputation qc(4);
  qc.cz(0, 1);
  qc.rz(qc::PI, 1);
  qc.rz(qc::PI, 2);
  qc.cz(2, 3);
  qc.cz(1, 2);
  const auto& [oneQubitGateLayers, twoQubitGateLayers] = scheduler.schedule(qc);
  EXPECT_THAT(oneQubitGateLayers,
              ::testing::ElementsAre(
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(2)))),
                  ::testing::ElementsAre(::testing::RefEq(
                      static_cast<qc::StandardOperation&>(*qc.at(1)))),
                  ::testing::IsEmpty()));
  EXPECT_THAT(
      twoQubitGateLayers,
      ::testing::ElementsAre(::testing::ElementsAre(::testing::Pair(0U, 1U),
                                                    ::testing::Pair(2U, 3U)),
                             ::testing::ElementsAre(::testing::Pair(1U, 2U))));
}
TEST(ASAPSchedulerTest, Config) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  const auto config = R"({
  "asap_scheduler": {
    "unknown_key": 42
  }
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = ASAPScheduler(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_EQ(buffer.str(), "[WARN] Configuration for ASAPScheduler contains an "
                          "unknown key: unknown_key. Ignoring.\n");
}
} // namespace na
