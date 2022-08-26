#include "Architecture.hpp"
#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

class TestCliffordSynthesis: public testing::TestWithParam<std::string> {
protected:
    std::string test_architecture_dir = "./architectures/";
    std::string test_calibration_dir  = "./calibration/";

    void SetUp() override {
        using namespace dd::literals;
    }
};

INSTANTIATE_TEST_SUITE_P(
        Architecture, TestCliffordSynthesis,
        testing::Values(
                "ibm_qx4.arch",
                "ibmq_casablanca.arch",
                "ibmq_london.arch",
                "ibmq_london.csv"));

TEST_P(TestCliffordSynthesis, LoadArchitecture) {
    auto&             arch_name = GetParam();
    Architecture      arch{};
    CliffordOptimizer opt{};

    std::stringstream ss{};
    if (arch_name.find(".arch") != std::string::npos) {
        ss << test_architecture_dir << arch_name;
        arch.loadCouplingMap(ss.str());
    } else {
        ss << test_calibration_dir << arch_name;
        arch.loadProperties(ss.str());
    }
    opt.setArchitecture(arch);

    EXPECT_EQ(opt.highestFidelityMap.size(), 1);
}
