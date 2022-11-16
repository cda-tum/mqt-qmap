/*
 * This file is part of the MQT QMAP library which is released under the MIT
 * license. See file README.md or go to
 * https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

#include "Architecture.hpp"

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
  std::string line = "Entry1;Entry2;\"EscapedEntry1;EscapedEntry2\";Entry3";

  std::vector<std::string> data{};
  parse_line(line, ';', {'\"'}, {'\\'}, data);

  EXPECT_EQ(data[1], "Entry2");
  EXPECT_EQ(data[2], "EscapedEntry1;EscapedEntry2");
}
