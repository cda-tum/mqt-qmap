/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "Architecture.hpp"

#include "gtest/gtest.h"
#include <iostream>

TEST(General, LoadCouplingMapNonexistentFile) {
    EXPECT_THROW(Architecture("path/that/does/not/exist"), QMAPException);
}

TEST(General, LoadCouplingMapEmptyFile) {
    std::ofstream ofs("test.arch");
    EXPECT_THROW(Architecture("test.arch"), QMAPException);
}

TEST(General, LoadCouplingMapNoQubitCount) {
    std::ofstream ofs("test.arch");
    ofs << "noqubits\n";
    EXPECT_THROW(Architecture("test.arch"), QMAPException);
}

TEST(General, LoadCouplingMapNoEdge) {
    std::ofstream ofs("test.arch");
    ofs << "1\n"
        << "noedge\n";
    EXPECT_THROW(Architecture("test.arch"), QMAPException);
}

TEST(General, LoadCalibrationDataNonexistentFile) {
    std::ofstream ofs("test.arch");
    ofs << "2\n"
        << "0 1\n";
    EXPECT_THROW(Architecture("test.arch", "path/that/does/not/exist"), QMAPException);
}
