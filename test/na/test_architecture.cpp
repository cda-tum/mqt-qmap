//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"
#include "operations/OpType.hpp"

#include "gtest/gtest.h"
#include <fstream>
#include <random>
#include <sstream>
#include <string>

constexpr const char* ARCH_FN = "arch.json";
constexpr const char* GRID_FN = "grid.csv";

class TestNAArchitecture : public testing::Test {
protected:
  void SetUp() override {
    // write content to a file
    std::ofstream archF(ARCH_FN);
    archF << "{\n"
             "  \"name\": \"Nature\",\n"
             "  \"initialZones\": [\n"
             "    \"storage\"\n"
             "  ],\n"
             "  \"zones\": [\n"
             "    {\n"
             "      \"name\": \"entangling\",\n"
             "      \"xmin\": -12,\n"
             "      \"xmax\": 126,\n"
             "      \"ymin\": -12,\n"
             "      \"ymax\": 31,\n"
             "      \"fidelity\": 0.999\n"
             "    },\n"
             "    {\n"
             "      \"name\": \"storage\",\n"
             "      \"xmin\": -12,\n"
             "      \"xmax\": 126,\n"
             "      \"ymin\": 32,\n"
             "      \"ymax\": 123,\n"
             "      \"fidelity\": 0.99\n"
             "    },\n"
             "    {\n"
             "      \"name\": \"readout\",\n"
             "      \"xmin\": -12,\n"
             "      \"xmax\": 126,\n"
             "      \"ymin\": 124,\n"
             "      \"ymax\": 160,\n"
             "      \"fidelity\": 0.95\n"
             "    }\n"
             "  ],\n"
             "  \"operations\": [\n"
             "    {\n"
             "      \"name\": \"rz\",\n"
             "      \"type\": \"local\",\n"
             "      \"zones\": [\n"
             "        \"storage\",\n"
             "        \"readout\"\n"
             "      ],\n"
             "      \"time\": 0.5,\n"
             "      \"fidelity\": 0.999\n"
             "    },\n"
             "    {\n"
             "      \"name\": \"ry\",\n"
             "      \"type\": \"global\",\n"
             "      \"zones\": [\n"
             "        \"storage\",\n"
             "        \"readout\"\n"
             "      ],\n"
             "      \"time\": 0.5,\n"
             "      \"fidelity\": 0.999\n"
             "    },\n"
             "    {\n"
             "      \"name\": \"cz\",\n"
             "      \"type\": \"global\",\n"
             "      \"zones\": [\n"
             "        \"entangling\"\n"
             "      ],\n"
             "      \"time\": 0.2,\n"
             "      \"fidelity\": 0.9959\n"
             "    },\n"
             "    {\n"
             "      \"name\": \"measure\",\n"
             "      \"type\": \"global\",\n"
             "      \"zones\": [\n"
             "        \"readout\"\n"
             "      ],\n"
             "      \"time\": 0.2,\n"
             "      \"fidelity\": 0.95\n"
             "    }\n"
             "  ],\n"
             "  \"decoherence\": {\n"
             "    \"t1\": 100000000,\n"
             "    \"t2\": 1500000\n"
             "  },\n"
             "  \"interactionRadius\": 3,\n"
             "  \"minAtomDistance\": 0.1,\n"
             "  \"AOD\": [\n"
             "    {\n"
             "      \"rows\": 5,\n"
             "      \"columns\": 5,\n"
             "      \"xmin\": -2.5,\n"
             "      \"xmax\": 2.5,\n"
             "      \"ymin\": -2.5,\n"
             "      \"ymax\": 2.5,\n"
             "      \"move\": {\n"
             "        \"speed\": 0.55,\n"
             "        \"fidelity\": 1\n"
             "      },\n"
             "      \"load\": {\n"
             "        \"time\": 20,\n"
             "        \"fidelity\": 1\n"
             "      },\n"
             "      \"store\": {\n"
             "        \"time\": 20,\n"
             "        \"fidelity\": 1\n"
             "      }\n"
             "    }\n"
             "  ]\n"
             "}\n";
    archF.close();
    std::ofstream gridF(GRID_FN);
    gridF << "x,y\n"
             "3,0\n"
             "15,0\n"
             "27,0\n"
             "39,0\n"
             "51,0\n"
             "63,0\n"
             "75,0\n"
             "87,0\n"
             "99,0\n"
             "111,0\n"
             "3,6\n"
             "15,6\n"
             "27,6\n"
             "39,6\n"
             "51,6\n"
             "63,6\n"
             "75,6\n"
             "87,6\n"
             "99,6\n"
             "111,6\n"
             "3,12\n"
             "15,12\n"
             "27,12\n"
             "39,12\n"
             "51,12\n"
             "63,12\n"
             "75,12\n"
             "87,12\n"
             "99,12\n"
             "111,12\n"
             "3,18\n"
             "15,18\n"
             "27,18\n"
             "39,18\n"
             "51,18\n"
             "63,18\n"
             "75,18\n"
             "87,18\n"
             "99,18\n"
             "111,18\n"
             "0,44\n"
             "6,44\n"
             "12,44\n"
             "18,44\n"
             "24,44\n"
             "30,44\n"
             "36,44\n"
             "42,44\n"
             "48,44\n"
             "54,44\n"
             "60,44\n"
             "66,44\n"
             "72,44\n"
             "78,44\n"
             "84,44\n"
             "90,44\n"
             "96,44\n"
             "102,44\n"
             "108,44\n"
             "114,44\n"
             "0,50\n"
             "6,50\n"
             "12,50\n"
             "18,50\n"
             "24,50\n"
             "30,50\n"
             "36,50\n"
             "42,50\n"
             "48,50\n"
             "54,50\n"
             "60,50\n"
             "66,50\n"
             "72,50\n"
             "78,50\n"
             "84,50\n"
             "90,50\n"
             "96,50\n"
             "102,50\n"
             "108,50\n"
             "114,50\n"
             "0,56\n"
             "6,56\n"
             "12,56\n"
             "18,56\n"
             "24,56\n"
             "30,56\n"
             "36,56\n"
             "42,56\n"
             "48,56\n"
             "54,56\n"
             "60,56\n"
             "66,56\n"
             "72,56\n"
             "78,56\n"
             "84,56\n"
             "90,56\n"
             "96,56\n"
             "102,56\n"
             "108,56\n"
             "114,56\n"
             "0,62\n"
             "6,62\n"
             "12,62\n"
             "18,62\n"
             "24,62\n"
             "30,62\n"
             "36,62\n"
             "42,62\n"
             "48,62\n"
             "54,62\n"
             "60,62\n"
             "66,62\n"
             "72,62\n"
             "78,62\n"
             "84,62\n"
             "90,62\n"
             "96,62\n"
             "102,62\n"
             "108,62\n"
             "114,62\n"
             "0,68\n"
             "6,68\n"
             "12,68\n"
             "18,68\n"
             "24,68\n"
             "30,68\n"
             "36,68\n"
             "42,68\n"
             "48,68\n"
             "54,68\n"
             "60,68\n"
             "66,68\n"
             "72,68\n"
             "78,68\n"
             "84,68\n"
             "90,68\n"
             "96,68\n"
             "102,68\n"
             "108,68\n"
             "114,68\n"
             "0,74\n"
             "6,74\n"
             "12,74\n"
             "18,74\n"
             "24,74\n"
             "30,74\n"
             "36,74\n"
             "42,74\n"
             "48,74\n"
             "54,74\n"
             "60,74\n"
             "66,74\n"
             "72,74\n"
             "78,74\n"
             "84,74\n"
             "90,74\n"
             "96,74\n"
             "102,74\n"
             "108,74\n"
             "114,74\n"
             "0,80\n"
             "6,80\n"
             "12,80\n"
             "18,80\n"
             "24,80\n"
             "30,80\n"
             "36,80\n"
             "42,80\n"
             "48,80\n"
             "54,80\n"
             "60,80\n"
             "66,80\n"
             "72,80\n"
             "78,80\n"
             "84,80\n"
             "90,80\n"
             "96,80\n"
             "102,80\n"
             "108,80\n"
             "114,80\n"
             "0,86\n"
             "6,86\n"
             "12,86\n"
             "18,86\n"
             "24,86\n"
             "30,86\n"
             "36,86\n"
             "42,86\n"
             "48,86\n"
             "54,86\n"
             "60,86\n"
             "66,86\n"
             "72,86\n"
             "78,86\n"
             "84,86\n"
             "90,86\n"
             "96,86\n"
             "102,86\n"
             "108,86\n"
             "114,86\n"
             "0,92\n"
             "6,92\n"
             "12,92\n"
             "18,92\n"
             "24,92\n"
             "30,92\n"
             "36,92\n"
             "42,92\n"
             "48,92\n"
             "54,92\n"
             "60,92\n"
             "66,92\n"
             "72,92\n"
             "78,92\n"
             "84,92\n"
             "90,92\n"
             "96,92\n"
             "102,92\n"
             "108,92\n"
             "114,92\n"
             "0,98\n"
             "6,98\n"
             "12,98\n"
             "18,98\n"
             "24,98\n"
             "30,98\n"
             "36,98\n"
             "42,98\n"
             "48,98\n"
             "54,98\n"
             "60,98\n"
             "66,98\n"
             "72,98\n"
             "78,98\n"
             "84,98\n"
             "90,98\n"
             "96,98\n"
             "102,98\n"
             "108,98\n"
             "114,98\n"
             "0,104\n"
             "6,104\n"
             "12,104\n"
             "18,104\n"
             "24,104\n"
             "30,104\n"
             "36,104\n"
             "42,104\n"
             "48,104\n"
             "54,104\n"
             "60,104\n"
             "66,104\n"
             "72,104\n"
             "78,104\n"
             "84,104\n"
             "90,104\n"
             "96,104\n"
             "102,104\n"
             "108,104\n"
             "114,104\n"
             "0,110\n"
             "6,110\n"
             "12,110\n"
             "18,110\n"
             "24,110\n"
             "30,110\n"
             "36,110\n"
             "42,110\n"
             "48,110\n"
             "54,110\n"
             "60,110\n"
             "66,110\n"
             "72,110\n"
             "78,110\n"
             "84,110\n"
             "90,110\n"
             "96,110\n"
             "102,110\n"
             "108,110\n"
             "114,110\n"
             "0,136\n"
             "6,136\n"
             "12,136\n"
             "18,136\n"
             "24,136\n"
             "30,136\n"
             "36,136\n"
             "42,136\n"
             "48,136\n"
             "54,136\n"
             "60,136\n"
             "66,136\n"
             "72,136\n"
             "78,136\n"
             "84,136\n"
             "90,136\n"
             "96,136\n"
             "102,136\n"
             "108,136\n"
             "114,136\n"
             "0,142\n"
             "6,142\n"
             "12,142\n"
             "18,142\n"
             "24,142\n"
             "30,142\n"
             "36,142\n"
             "42,142\n"
             "48,142\n"
             "54,142\n"
             "60,142\n"
             "66,142\n"
             "72,142\n"
             "78,142\n"
             "84,142\n"
             "90,142\n"
             "96,142\n"
             "102,142\n"
             "108,142\n"
             "114,142\n"
             "0,148\n"
             "6,148\n"
             "12,148\n"
             "18,148\n"
             "24,148\n"
             "30,148\n"
             "36,148\n"
             "42,148\n"
             "48,148\n"
             "54,148\n"
             "60,148\n"
             "66,148\n"
             "72,148\n"
             "78,148\n"
             "84,148\n"
             "90,148\n"
             "96,148\n"
             "102,148\n"
             "108,148\n"
             "114,148\n"
             "0,154\n"
             "6,154\n"
             "12,154\n"
             "18,154\n"
             "24,154\n"
             "30,154\n"
             "36,154\n"
             "42,154\n"
             "48,154\n"
             "54,154\n"
             "60,154\n"
             "66,154\n"
             "72,154\n"
             "78,154\n"
             "84,154\n"
             "90,154\n"
             "96,154\n"
             "102,154\n"
             "108,154\n"
             "114,154\n"
             "\n";
    gridF.close();
  }
  void TearDown() override {
    std::remove(ARCH_FN);
    std::remove(GRID_FN);
  }
};

TEST_F(TestNAArchitecture, Import) {
  na::Architecture const arch(ARCH_FN, GRID_FN);

  EXPECT_EQ(arch.getNZones(), 3);
  EXPECT_EQ(arch.getNSites(), 360);
  EXPECT_EQ(arch.getName(), "Nature");
  EXPECT_EQ(arch.getZoneLabel(arch.getZoneOfSite(0)), "entangling");
}

TEST_F(TestNAArchitecture, GateApplicability) {
  na::Architecture const arch(ARCH_FN, GRID_FN);

  EXPECT_TRUE(arch.isAllowedGlobally({qc::OpType::RY, 0}, 1));
  EXPECT_TRUE(arch.isAllowedGlobally({qc::OpType::Z, 1}, 0));
  EXPECT_TRUE(arch.isAllowedLocally({qc::OpType::RZ, 0}, 1));
}

TEST_F(TestNAArchitecture, WithConfiguration) {
  na::Architecture const  arch(ARCH_FN, GRID_FN);
  na::Configuration const config(2, 3);
  const auto              modArch = arch.withConfig(config);
  EXPECT_EQ(modArch.getNSites(), 54);
  for (std::size_t z = 0; z < modArch.getNZones(); ++z) {
    std::cout << "Zone: " << modArch.getZoneLabel(z) << std::endl;
    for (const auto s : modArch.getSitesInZone(z)) {
      std::cout << "  Site: " << s << "(" << modArch.getPositionOfSite(s) << ")" << std::endl;
    }
  }
}