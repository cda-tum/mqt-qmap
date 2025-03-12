#include "na/azac/ASAPScheduler.hpp"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <optional>
#include <utility>
#include <vector>

namespace na {
class ASAPSchedulerScheduleTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  ASAPScheduler scheduler;

public:
  ASAPSchedulerScheduleTest() : scheduler{architecture, config} {}
};
// TEST_F(ASAPSchedulerScheduleTest, NoGates) {
//   std::vector<std::vector<std::pair<qc::Qubit, qc::Qubit>>>
//   twoQubitGateLayers; EXPECT_NO_THROW(
//       { EXPECT_TRUE(scheduler.analyzeReuse(twoQubitGateLayers).empty()); });
// }
// TEST(VMReuseAnalyzerTest, Config) {
//   Architecture architecture;
//   nlohmann::json config;
//   std::istringstream iss(R"({
//   "asap_scheduler": {
//     "unknown_key": 42
//   }
// })");
//   iss >> config;
//   std::stringstream buffer;
//   std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
//   ASAPScheduler scheduler{architecture, config};
//   std::cout.rdbuf(oldCout);
//   EXPECT_EQ(buffer.str(), "[WARN] Configuration for ASAPScheduler contains an
//   unknown "
//                           "key: unknown_key. Ignoring.\n");
//   // silence unused variable warning
//   (void)analyzer;
// }
} // namespace na
