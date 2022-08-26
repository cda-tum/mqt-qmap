#include "Architecture.hpp"
#include "Tableau.hpp"
#include "QuantumComputation.hpp"

#include "gtest/gtest.h"

class TestTableau: public testing::TestWithParam<std::string> {
protected:
    std::string test_architecture_dir = "./architectures/";
    std::string test_calibration_dir  = "./calibration/";

    void SetUp() override {
        using namespace dd::literals;
    }
};

TEST(TestTableau, LoadTableau) {
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

TEST(TestTableau, GetStrRepresentation) {
    Tableau     tableau{};
    Tableau     tableau2{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.importString(tableau_string);

    std::string result_string = "2|1|2|3|4|\n1|0|0|0|1|0|\n2|0|0|1|0|0|\n";

    EXPECT_EQ(tableau.getStrRepresentation(), result_string);
}

TEST(TestTableau, AccessValues) {
    Tableau     tableau{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.importString(tableau_string);

    EXPECT_EQ(tableau[0][0], 0);
    EXPECT_EQ(tableau[0][1], 0);
    EXPECT_EQ(tableau[0][3], 1);

    auto val = tableau.at(0);

    EXPECT_EQ(val[0], 0);
    EXPECT_EQ(val[1], 0);
    EXPECT_EQ(val[3], 1);

    auto tableau2 = tableau;

    EXPECT_EQ(tableau2[0][0], 0);
    EXPECT_EQ(tableau2[0][1], 0);
    EXPECT_EQ(tableau2[0][3], 1);

    EXPECT_EQ(tableau.back(), tableau2.back());
}

TEST(TestTableau, BasicFunctions) {
    Tableau     tableau{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.importString(tableau_string);

    tableau.clear();

    EXPECT_EQ(tableau.empty(), true);

    tableau.importString(tableau_string);

    auto tableau2 = Tableau::getDiagonalTableau(2);

    EXPECT_EQ(tableau2[0][0], 0);
    EXPECT_EQ(tableau2[0][1], 0);
    EXPECT_EQ(tableau2[0][2], 1);
    EXPECT_EQ(tableau2[0][3], 0);
    EXPECT_EQ(tableau2[1][0], 0);
    EXPECT_EQ(tableau2[1][1], 0);
    EXPECT_EQ(tableau2[1][2], 0);
    EXPECT_EQ(tableau2[1][3], 1);
}

TEST(TestTableau, LoadTableauFromQC) {
    using namespace dd::literals;

    auto qc = qc::QuantumComputation(2U);
    qc.h(0);
    qc.x(1, 0_pc);

    Tableau     tableau{};
    Tableau tableau1{};
    std::string tableau_string = "Stabilizer = ['+XX', '+ZZ']";

    tableau.importString(tableau_string);
    Tableau::generateTableau(tableau1, qc);

    EXPECT_EQ(tableau, tableau1);

    tableau1.clear();

    tableau.import("cliffordexamples/basic.qasm");

    EXPECT_EQ(tableau, tableau1);
}