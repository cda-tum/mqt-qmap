#include "Architecture.hpp"

#include "gtest/gtest.h"

class TestArchitecture: public testing::TestWithParam<std::string> {
protected:
    std::string test_architecture_dir = "./architectures/";
    std::string test_calibration_dir  = "./calibration/";

    void SetUp() override {
        using namespace dd::literals;
    }
};

INSTANTIATE_TEST_SUITE_P(
        Architecture, TestArchitecture,
        testing::Values(
                "ibm_qx4.arch",
                "ibmq_casablanca.arch",
                "ibmq_london.arch",
                "ibmq_london.csv"));

TEST_P(TestArchitecture, QubitMap) {
    auto& arch_name = GetParam();
    Architecture arch{};
    std::stringstream ss{};
    if (arch_name.find(".arch")!=std::string::npos){
        ss << test_architecture_dir << arch_name;
        arch.loadCouplingMap(ss.str());
    } else {
        ss << test_calibration_dir << arch_name;
        arch.loadCalibrationData(ss.str());
    }

    EXPECT_EQ(Architecture::getQubitMap(arch.getCouplingMap()).size(), arch.getNqubits());
}
TEST_P(TestArchitecture, GetAllConnectedSubsets) {
    auto& arch_name = GetParam();
    Architecture arch{};
    std::stringstream ss{};
    if (arch_name.find(".arch")!=std::string::npos){
        ss << test_architecture_dir << arch_name;
        arch.loadCouplingMap(ss.str());
    } else {
        ss << test_calibration_dir << arch_name;
        arch.loadCalibrationData(ss.str());
    }

    EXPECT_EQ(arch.getAllConnectedSubsets(arch.getNqubits()).size(), 1);
    EXPECT_EQ(arch.getAllConnectedSubsets(1).size(), arch.getNqubits());
}
TEST_P(TestArchitecture, GetHighestFidelity) {
    auto& arch_name = GetParam();
    Architecture arch{};
    std::stringstream ss{};
    if (arch_name.find(".arch")!=std::string::npos){
        ss << test_architecture_dir << arch_name;
        arch.loadCouplingMap(ss.str());
    } else {
        ss << test_calibration_dir << arch_name;
        arch.loadCalibrationData(ss.str());
    }
    CouplingMap cm{};

    arch.getHighestFidelityCouplingMap(arch.getNqubits(), cm);

    EXPECT_EQ(cm, arch.getCouplingMap());

    arch.getHighestFidelityCouplingMap(1, cm);

    CouplingMap expected {};

    if (arch_name.find(".csv")==std::string::npos){
        EXPECT_EQ(cm, arch.getCouplingMap());
    } else {
        EXPECT_NE(cm, arch.getCouplingMap());
    }
}