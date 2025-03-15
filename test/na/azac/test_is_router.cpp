#include "na/azac/ISRouter.hpp"

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
  "name": "is_router_architecture",
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
class ISRouterPlaceTest : public ::testing::Test {
protected:
  Architecture architecture;
  nlohmann::json config;
  ISRouter router;
  ISRouterPlaceTest()
      : architecture(nlohmann::json::parse(architectureJson)),
        router(architecture, config) {}
};
TEST_F(ISRouterPlaceTest, Empty) {
  EXPECT_THAT(
      router.route(
          std::vector<std::vector<std::tuple<std::reference_wrapper<const SLM>,
                                             size_t, size_t>>>{}),
      ::testing::IsEmpty());
}
TEST_F(ISRouterPlaceTest, Initial) {
  const auto& slm = *architecture.storageZones.front();
  EXPECT_THAT(
      router.route(
          std::vector<std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>{
              {{slm, 0, 0}}}),
      ::testing::IsEmpty());
}
TEST_F(ISRouterPlaceTest, OneLayer) {
  // STORAGE     ...         │ ...         │ ...
  //         18  o o o o ... │ o o o o ... │ o o o o ...
  //         19  0 1 o o ... │ o o o o ... │ 0 1 o o ...
  //                         │  ╲╲         │ ↑ ↑
  // ENTANGLEMENT            │   ↓↓        │  ╲ ╲
  //          0    oo    ... │   01    ... │   oo    ...
  //          1    oo    ... │   oo    ... │   oo    ...
  //               ...       │   ...       │   ...
  const auto& storage = *architecture.storageZones.front();
  const auto& entanglementLeft =
      architecture.entanglementZones.front()->front();
  const auto& entanglementRight =
      architecture.entanglementZones.front()->back();
  EXPECT_THAT(
      router.route(
          std::vector<std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>{
              {{storage, 19, 0}, {storage, 19, 1}},
              {{entanglementLeft, 0, 0}, {entanglementRight, 0, 0}},
              {{storage, 19, 0}, {storage, 19, 1}}}),
      ::testing::ElementsAre(::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(0U, 1U)),
                             ::testing::UnorderedElementsAre(
                                 ::testing::UnorderedElementsAre(0U, 1U))));
}
TEST_F(ISRouterPlaceTest, Cross) {
  // STORAGE     ...         │ ...
  //         18  o o o o ... │ o o o o ...
  //         19  0 1 o o ... │ o o o o ...
  //                         │  ╲|
  // ENTANGLEMENT            │   ↓↘
  //          0    oo    ... │   10    ...
  //          1    oo    ... │   oo    ...
  //               ...       │   ...
  const auto& storage = *architecture.storageZones.front();
  const auto& entanglementLeft =
      architecture.entanglementZones.front()->front();
  const auto& entanglementRight =
      architecture.entanglementZones.front()->back();
  EXPECT_THAT(
      router.route(
          std::vector<std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>{
              {{storage, 19, 0}, {storage, 19, 1}},
              {{entanglementRight, 0, 0}, {entanglementLeft, 0, 0}}}),
      ::testing::ElementsAre(::testing::UnorderedElementsAre(
          ::testing::UnorderedElementsAre(0U),
          ::testing::UnorderedElementsAre(1U))));
}
TEST_F(ISRouterPlaceTest, Overtake) {
  // STORAGE     ...         │ ...
  //         18  0 1 o o ... │ o o o o ...
  //         19  2 3 o o ... │ o o o o ...
  //                         │  ╲╲
  // ENTANGLEMENT            │   ↓↓
  //          0    oo    ... │   23    ...
  //          1    oo    ... │   01    ...
  //               ...       │   ...
  const auto& storage = *architecture.storageZones.front();
  const auto& entanglementLeft =
      architecture.entanglementZones.front()->front();
  const auto& entanglementRight =
      architecture.entanglementZones.front()->back();
  EXPECT_THAT(
      router.route(
          std::vector<std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>{
              {{storage, 18, 0},
               {storage, 18, 1},
               {storage, 19, 0},
               {storage, 19, 1}},
              {{entanglementLeft, 1, 0},
               {entanglementRight, 1, 0},
               {entanglementLeft, 0, 0},
               {entanglementRight, 0, 0}}}),
      ::testing::ElementsAre(::testing::UnorderedElementsAre(
          ::testing::UnorderedElementsAre(0U, 1U),
          ::testing::UnorderedElementsAre(2U, 3U))));
}
TEST_F(ISRouterPlaceTest, Array) {
  // STORAGE     ...             │ ...
  //         18  0 1 2 3 o o ... │ o o o o o o ...
  //         19  4 5 6 7 o o ... │ o o o o o o ...
  //                             │  ╲╲   ╲╲
  // ENTANGLEMENT                │   ↓↓    ↘↘
  //          0    oo     oo ... │   01     23 ...
  //          1    oo     oo ... │   45     67 ...
  //               ...           │   ...
  const auto& storage = *architecture.storageZones.front();
  const auto& entanglementLeft =
      architecture.entanglementZones.front()->front();
  const auto& entanglementRight =
      architecture.entanglementZones.front()->back();
  EXPECT_THAT(
      router.route(
          std::vector<std::vector<
              std::tuple<std::reference_wrapper<const SLM>, size_t, size_t>>>{
              {{storage, 18, 0},
               {storage, 18, 1},
               {storage, 18, 2},
               {storage, 18, 3},
               {storage, 19, 0},
               {storage, 19, 1},
               {storage, 19, 2},
               {storage, 19, 3}},
              {{entanglementLeft, 0, 0},
               {entanglementRight, 0, 0},
               {entanglementLeft, 0, 1},
               {entanglementRight, 0, 1},
               {entanglementLeft, 1, 0},
               {entanglementRight, 1, 0},
               {entanglementLeft, 1, 1},
               {entanglementRight, 1, 1}}}),
      ::testing::ElementsAre(::testing::UnorderedElementsAre(
          ::testing::UnorderedElementsAre(0U, 1U, 2U, 3U, 4U, 5U, 6U, 7U))));
}
TEST(ISRouterTest, InvalidConfig) {
  Architecture architecture(nlohmann::json::parse(architectureJson));
  nlohmann::json config = R"({
  "is_router": {
    "unknown_key": 42
  }
})"_json;
  std::stringstream buffer;
  std::streambuf* oldCout = std::cout.rdbuf(buffer.rdbuf());
  std::ignore = ISRouter(architecture, config);
  std::cout.rdbuf(oldCout);
  EXPECT_EQ(buffer.str(), "[WARN] Configuration for ISRouter contains an "
                          "unknown key: unknown_key. Ignoring.\n");
}
} // namespace na
