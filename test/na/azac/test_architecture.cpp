#include "na/azac/Architecture.hpp"

#include <gtest/gtest.h>
#include <sstream>

namespace na {

class TestArchitecture : public testing::Test {
protected:
  Architecture arch;
  void SetUp() override {
    std::istringstream archIS(R"({
  "name": "full_compute_store_architecture",
  "operation_duration": {"rydberg": 0.36, "1qGate": 52, "atom_transfer": 15},
  "operation_fidelity": {
    "two_qubit_gate": 0.995,
    "single_qubit_gate": 0.9997,
    "atom_transfer": 0.999
  },
  "qubit_spec": {"T": 1.5e6},
  "storage_zones": [{
    "zone_id": 0,
    "slms": [{
      "id": 0,
      "site_seperation": [3, 3],
      "r": 100,
      "c": 100,
      "location": [0, 0]}],
    "offset": [0, 0],
    "dimenstion": [300, 300]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {
        "id": 1,
        "site_seperation": [12, 10],
        "r": 7,
        "c": 20,
        "location": [35, 307]
      },
      {
        "id": 2,
        "site_seperation": [12, 10],
        "r": 7,
        "c": 20,
        "location": [37, 307]
      }],
    "offset": [35, 307],
    "dimension": [240, 70]
  }],
  "aods":[{"id": 0, "site_seperation": 2, "r": 100, "c": 100}],
  "arch_range": [[0, 0], [297, 402]],
  "rydberg_range": [[[5, 305], [292, 402]]]
})");
    ASSERT_NO_THROW(arch.load(archIS));
    ASSERT_NO_THROW(arch.preprocessing());
  }
};

TEST_F(TestArchitecture, Load) {}

TEST_F(TestArchitecture, Storage) {
  EXPECT_EQ(arch.storage_zone.size(), 1);
  EXPECT_EQ(arch.storage_zone.front()->n_r, 100);
  EXPECT_EQ(arch.storage_zone.front()->n_c, 100);
}

} // namespace na
