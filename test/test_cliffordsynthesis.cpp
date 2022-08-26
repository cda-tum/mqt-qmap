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

TEST(TestCliffordSynthesis, LoadTableau) {
    Tableau     tableau{};
    Tableau     tableau2{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.importString(tableau_string);

    tableau_string = "Stabilizer = ['+IZ', '+ZI']";

    tableau2.importString(tableau_string);

    EXPECT_EQ(tableau, tableau2);

    tableau_string = "Destabilizer = ['+IZ', '+XI']";

    tableau.clear();
    tableau.importString(tableau_string);

    tableau_string = "Stabilizer = ['+IX', '+ZI']";

    tableau2.clear();
    tableau2.importString(tableau_string);

    EXPECT_EQ(tableau, tableau2);
}

TEST(TestCliffordSynthesis, GetStrRepresentation) {
    Tableau     tableau{};
    Tableau     tableau2{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.importString(tableau_string);

    std::string result_string = "Tableau: 2|1|2|3|4|\n1|0|0|1|0|0|\n2|0|0|0|1|0|";

    EXPECT_EQ(tableau.getStrRepresentation(), result_string);
}