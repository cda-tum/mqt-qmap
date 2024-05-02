//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"
#include "Configuration.hpp"
#include "operations/OpType.hpp"

#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <optional>
#include <string>

class NAArchitecture : public testing::Test {
protected:
  na::Architecture arch;
  void             SetUp() override {
    // write content to a file
    std::istringstream archIS(R"({
      "name": "Nature",
      "initialZones": [
          "storage"
      ],
      "zones": [
          {
              "name": "entangling",
              "xmin": -300,
              "xmax": 656,
              "ymin": -10,
              "ymax": 46,
              "fidelity": 0.9959
          },
          {
              "name": "storage",
              "xmin": -300,
              "xmax": 656,
              "ymin": 47,
              "ymax": 121,
              "fidelity": 1
          },
          {
              "name": "readout",
              "xmin": -300,
              "xmax": 656,
              "ymin": 122,
              "ymax": 156,
              "fidelity": 0.99
          }
      ],
      "operations": [
          {
              "name": "rz",
              "type": "local",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "ry",
              "type": "global",
              "zones": [
                  "entangling",
                  "storage",
                  "readout"
              ],
              "time": 0.5,
              "fidelity": 0.999
          },
          {
              "name": "cz",
              "type": "global",
              "zones": [
                  "entangling"
              ],
              "time": 0.2,
              "fidelity": 0.9959
          },
          {
              "name": "measure",
              "type": "global",
              "zones": [
                  "readout"
              ],
              "time": 0.2,
              "fidelity": 0.95
          }
      ],
      "decoherence": {
          "t1": 100000000,
          "t2": 1500000
      },
      "interactionRadius": 2,
      "noInteractionRadius": 5,
      "minAtomDistance": 1,
      "shuttling": [
          {
              "rows": 5,
              "columns": 5,
              "xmin": -2.5,
              "xmax": 2.5,
              "ymin": -2.5,
              "ymax": 2.5,
              "move": {
                  "speed": 0.55,
                  "fidelity": 1
              },
              "load": {
                  "time": 20,
                  "fidelity": 1
              },
              "store": {
                  "time": 20,
                  "fidelity": 1
              }
          }
      ]
  })");
    std::stringstream  gridSS;
    gridSS << "x,y\n";
    // entangling zone (4 x 36 = 144 sites)
    for (std::size_t y = 0; y <= 36; y += 12) {
      for (std::size_t x = 3; x <= 353; x += 10) {
        gridSS << x << "," << y << "\n";
      }
    }
    // storage zone (12 x 72 = 864 sites)
    for (std::size_t y = 56; y <= 111; y += 5) {
      for (std::size_t x = 0; x <= 355; x += 5) {
        gridSS << x << "," << y << "\n";
      }
    }
    // readout zone (4 x 72 = 288 sites)
    for (std::size_t y = 131; y <= 146; y += 5) {
      for (std::size_t x = 0; x <= 355; x += 5) {
        gridSS << x << "," << y << "\n";
      }
    }
    // total: 1296 sites
    arch.fromFileStream(archIS, gridSS);
  }
};

TEST_F(NAArchitecture, ScopeString) {
  EXPECT_EQ(na::getScopeOfString("local"), na::Scope::Local);
  EXPECT_EQ(na::getScopeOfString("gLoBaL"), na::Scope::Global);
  EXPECT_THROW(na::getScopeOfString(""), std::invalid_argument);
}

TEST_F(NAArchitecture, Import) {
  EXPECT_EQ(arch.getNZones(), 3);
  EXPECT_EQ(arch.getNSites(), 1296);
  EXPECT_EQ(arch.getName(), "Nature");
  EXPECT_EQ(arch.getZoneLabel(arch.getZoneOfSite(0)), "entangling");
  EXPECT_ANY_THROW(
      na::Architecture("file_does_not_Exist.json", "file_does_not_Exist.csv"));
  std::ofstream tmp("temp.json");
  tmp.close();
  EXPECT_ANY_THROW(na::Architecture("temp.json", "file_does_not_Exist.csv"));
  std::remove("temp.json");
  {
    std::istringstream archIS("{}");
    std::istringstream gridIS("x,y\n0;0");
    EXPECT_ANY_THROW(na::Architecture(archIS, gridIS));
  }
  {
    std::istringstream archIS("{");
    std::istringstream gridIS("x,y\n0,0");
    EXPECT_ANY_THROW(na::Architecture(archIS, gridIS));
  }
}

TEST_F(NAArchitecture, GateApplicability) {
  EXPECT_TRUE(arch.isAllowedGlobally({qc::OpType::RY, 0}, 1));
  EXPECT_TRUE(arch.isAllowedGlobally({qc::OpType::Z, 1}, 0));
  EXPECT_TRUE(arch.isAllowedLocally({qc::OpType::RZ, 0}, 1));
}

TEST_F(NAArchitecture, GateProperty) {
  EXPECT_EQ(arch.getPropertiesOfOperation({qc::OpType::RY, 0}).scope,
            na::Scope::Global);
  EXPECT_EQ(arch.getPropertiesOfOperation({qc::OpType::Z, 1}).fidelity, 0.9959);
  EXPECT_THROW(std::ignore = arch.getPropertiesOfOperation({
                   qc::OpType::RX,
               }),
               std::invalid_argument);
}

TEST_F(NAArchitecture, WithConfiguration) {
  na::Configuration const config(2, 3);
  const auto              modArch = arch.withConfig(config);
  EXPECT_EQ(modArch.getNSites(), 216);
}

TEST_F(NAArchitecture, SiteUp) {
  EXPECT_TRUE(arch.hasSiteUp({3, 3}, false, true).second);
  EXPECT_FALSE(arch.hasSiteUp({3, 0}, true, true).second);
  EXPECT_EQ(arch.getPositionOfSite(*arch.getNearestSiteUp({3, 3}, true, true)),
            (na::Point{3, 0}));
  EXPECT_EQ(arch.getNearestSiteUp({3, 0}, true, true), std::nullopt);
}

TEST_F(NAArchitecture, SiteDown) {
  EXPECT_FALSE(arch.hasSiteDown({0, 3}, false, true).second);
  EXPECT_TRUE(arch.hasSiteDown({3, 0}, true, true).second);
  EXPECT_EQ(arch.getNearestSiteDown({0, 3}, false, true), std::nullopt);
  EXPECT_EQ(
      arch.getPositionOfSite(*arch.getNearestSiteDown({3, 0}, true, true)),
      (na::Point{3, 12}));
}

TEST_F(NAArchitecture, SiteLeft) {
  EXPECT_TRUE(arch.hasSiteLeft({3, 0}, false, true).second);
  EXPECT_FALSE(arch.hasSiteLeft({3, 0}, true, true).second);
  EXPECT_EQ(
      arch.getPositionOfSite(*arch.getNearestSiteLeft({3, 0}, false, true)),
      (na::Point{3, 0}));
  EXPECT_EQ(arch.getNearestSiteLeft({3, 0}, true, true), std::nullopt);
}

TEST_F(NAArchitecture, SiteRight) {
  EXPECT_TRUE(arch.hasSiteRight({3, 0}, false, true).second);
  EXPECT_FALSE(arch.hasSiteRight({3, 3}, true, true).second);
  EXPECT_EQ(
      arch.getPositionOfSite(*arch.getNearestSiteRight({3, 0}, true, true)),
      (na::Point{13, 0}));
  EXPECT_EQ(arch.getNearestSiteRight({3, 3}, true, true), std::nullopt);
}
