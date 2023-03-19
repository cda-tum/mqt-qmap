//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"
#include "utils.hpp"

#include "gtest/gtest.h"
#include <iostream>

TEST(General, LoadCouplingMapNonexistentFile) {
  EXPECT_THROW(Architecture("path/that/does/not/exist"), QMAPException);
}

TEST(General, LoadCouplingMapEmptyFile) {
  std::ofstream ofs("test.arch");
  ofs.close();
  EXPECT_THROW(Architecture("test.arch"), QMAPException);
}

TEST(General, LoadCouplingMapNoQubitCount) {
  std::ofstream ofs("test.arch");
  ofs << "noqubits\n";
  ofs.close();
  EXPECT_THROW(Architecture("test.arch"), QMAPException);
}

TEST(General, LoadCouplingMapNoEdge) {
  std::ofstream ofs("test.arch");
  ofs << "1\n"
      << "noedge\n";
  ofs.close();
  EXPECT_THROW(Architecture("test.arch"), QMAPException);
}

TEST(General, LoadCalibrationDataNonexistentFile) {
  std::ofstream ofs("test.arch");
  ofs << "2\n"
      << "0 1\n";
  ofs.close();
  EXPECT_THROW(Architecture("test.arch", "path/that/does/not/exist"),
               QMAPException);
}

TEST(General, TestLineParsing) {
  const std::string line =
      "Entry1;Entry2;\"EscapedEntry1;EscapedEntry2\";Entry3";

  std::vector<std::string> data{};
  parseLine(line, ';', {'\"'}, {'\\'}, data);

  EXPECT_EQ(data[1], "Entry2");
  EXPECT_EQ(data[2], "EscapedEntry1;EscapedEntry2");
}

TEST(General, Dijkstra) {
  /*
            (6)
        .----------.
       \/          \/
  0 <-> 1 --> 2 <-> 3
    (1)   (2)   (3)
  */
  CouplingMap cm = {{0,1},{1,0},{1,2},{2,3},{3,2},{1,3},{3,1}};
  Matrix edgeWeights = {
    {0,1,0,0},
    {1,0,2,6},
    {0,15,0,3},
    {0,6,3,0}
  };
  Matrix distanceTable{};
  Dijkstra::buildTable(4, cm, distanceTable, edgeWeights, [](const Dijkstra::Node& n) {
    return n.cost;
  });
  EXPECT_EQ(distanceTable.size(), 4);
  EXPECT_EQ(distanceTable[0].size(), 4);
  EXPECT_EQ(distanceTable[0][0], 0);
  EXPECT_EQ(distanceTable[0][1], 1);
  EXPECT_EQ(distanceTable[0][2], 3);
  EXPECT_EQ(distanceTable[0][3], 6);
  EXPECT_EQ(distanceTable[1].size(), 4);
  EXPECT_EQ(distanceTable[1][0], 1);
  EXPECT_EQ(distanceTable[1][1], 0);
  EXPECT_EQ(distanceTable[1][2], 2);
  EXPECT_EQ(distanceTable[1][3], 5);
  EXPECT_EQ(distanceTable[2].size(), 4);
  EXPECT_EQ(distanceTable[2][0], 10);
  EXPECT_EQ(distanceTable[2][1], 9);
  EXPECT_EQ(distanceTable[2][2], 0);
  EXPECT_EQ(distanceTable[2][3], 3);
  EXPECT_EQ(distanceTable[3].size(), 4);
  EXPECT_EQ(distanceTable[3][0], 7);
  EXPECT_EQ(distanceTable[3][1], 6);
  EXPECT_EQ(distanceTable[3][2], 3);
  EXPECT_EQ(distanceTable[3][3], 0);
  
  distanceTable = {};
  Dijkstra::buildTable(4, cm, distanceTable, edgeWeights, [](const Dijkstra::Node& n) {
    return n.prevCost;
  });
  EXPECT_EQ(distanceTable.size(), 4);
  EXPECT_EQ(distanceTable[0].size(), 4);
  EXPECT_EQ(distanceTable[0][0], 0);
  EXPECT_EQ(distanceTable[0][1], 0);
  EXPECT_EQ(distanceTable[0][2], 1);
  EXPECT_EQ(distanceTable[0][3], 3);
  EXPECT_EQ(distanceTable[1].size(), 4);
  EXPECT_EQ(distanceTable[1][0], 0);
  EXPECT_EQ(distanceTable[1][1], 0);
  EXPECT_EQ(distanceTable[1][2], 0);
  EXPECT_EQ(distanceTable[1][3], 2);
  EXPECT_EQ(distanceTable[2].size(), 4);
  EXPECT_EQ(distanceTable[2][0], 9);
  EXPECT_EQ(distanceTable[2][1], 3);
  EXPECT_EQ(distanceTable[2][2], 0);
  EXPECT_EQ(distanceTable[2][3], 0);
  EXPECT_EQ(distanceTable[3].size(), 4);
  EXPECT_EQ(distanceTable[3][0], 6);
  EXPECT_EQ(distanceTable[3][1], 0);
  EXPECT_EQ(distanceTable[3][2], 0);
  EXPECT_EQ(distanceTable[3][3], 0);
}
