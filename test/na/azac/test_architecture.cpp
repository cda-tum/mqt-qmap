#include "na/azac/Architecture.hpp"
#include "na/azac/Utils.hpp"

#include <cstddef>
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
      "site_separation": [3, 3],
      "r": 100,
      "c": 100,
      "location": [0, 0]}],
    "offset": [0, 0],
    "dimension": [300, 300]
  }],
  "entanglement_zones": [{
    "zone_id": 0,
    "slms": [
      {
        "id": 1,
        "site_separation": [12, 10],
        "r": 7,
        "c": 20,
        "location": [35, 307]
      },
      {
        "id": 2,
        "site_separation": [12, 10],
        "r": 7,
        "c": 20,
        "location": [37, 307]
      }],
    "offset": [35, 307],
    "dimension": [240, 70]
  }],
  "aods":[{"id": 0, "site_separation": 2, "r": 100, "c": 100}],
  "arch_range": [[0, 0], [297, 402]],
  "rydberg_range": [[[5, 305], [292, 402]]]
})");
    ASSERT_NO_THROW(arch.load(archIS));
    ASSERT_NO_THROW(arch.preprocessing());
  }
};

TEST_F(TestArchitecture, Load) {}

TEST_F(TestArchitecture, Storage) {
  EXPECT_EQ(arch.storageZones.size(), 1);
  EXPECT_EQ(arch.storageZones.front()->nRows, 100);
  EXPECT_EQ(arch.storageZones.front()->nCols, 100);
}

TEST_F(TestArchitecture, Distance) {
  const auto& slm1 = *arch.storageZones.front();
  EXPECT_EQ(arch.distance(slm1, 0, 0, slm1, 0, 1), slm1.siteSeparation.first);
  EXPECT_EQ(arch.distance(slm1, 0, 0, slm1, 1, 0), slm1.siteSeparation.second);

  const auto& slm2 = *arch.entanglementZones.front().front();
  EXPECT_EQ(arch.distance(slm1, 0, 0, slm2, 0, 0),
            distance(slm1.location, slm2.location));
}

TEST_F(TestArchitecture, NearestStorageSite) {
  const auto& entanglementSLM = *arch.entanglementZones.front().front();
  const auto nearestStorageSite =
      arch.nearestStorageSite(entanglementSLM, 0, 0);
  const auto minDistance =
      arch.distance({&entanglementSLM, 0, 0}, nearestStorageSite);
  for (const auto& slm : arch.storageZones) {
    for (std::size_t r = 0; r < slm->nRows; ++r) {
      for (std::size_t c = 0; c < slm->nCols; ++c) {
        const auto distance =
            arch.distance({&entanglementSLM, 0, 0}, {slm.get(), r, c});
        EXPECT_GE(distance, minDistance);
      }
    }
  }
}
TEST_F(TestArchitecture, NearestEntanglementSite) {
  const auto& storageSlm = *arch.storageZones.front();
  const auto nearestEntanglementSite =
      arch.nearestEntanglementSite(storageSlm, 0, 0, storageSlm, 0, 1);
  const auto minDistance =
      arch.distance({&storageSlm, 0, 0}, nearestEntanglementSite) +
      arch.distance({&storageSlm, 0, 1}, nearestEntanglementSite);
  for (const auto& slms : arch.entanglementZones) {
    for (const auto& slm : slms) {
      for (std::size_t r = 0; r < slm->nRows; ++r) {
        for (std::size_t c = 0; c < slm->nCols; ++c) {
          const auto distance =
              arch.distance({&storageSlm, 0, 0}, {slm.get(), r, c}) +
              arch.distance({&storageSlm, 0, 1}, {slm.get(), r, c});
          EXPECT_GE(distance, minDistance);
        }
      }
    }
  }
}

TEST_F(TestArchitecture, ExportNoThrow) {
  ASSERT_NO_THROW(arch.exportNAVizMachine(arch.name + ".namachine"));
}

} // namespace na
