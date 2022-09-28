#include "Architecture.hpp"
#include "QuantumComputation.hpp"
#include "cliffordsynthesis/Tableau.hpp"

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

    tableau.fromString(tableau_string);

    tableau_string = "Stabilizer = ['+IZ', '+ZI']";

    tableau2.fromString(tableau_string);

    EXPECT_EQ(tableau, tableau2);

    tableau_string = "Destabilizer = ['+IZ', '+XI']";

    tableau.clear();
    tableau.fromString(tableau_string);

    tableau_string = "Stabilizer = ['+IX', '+ZI']";

    tableau2.clear();
    tableau2.fromString(tableau_string);

    EXPECT_EQ(tableau, tableau2);
}

TEST(TestTableau, GetStrRepresentation) {
    Tableau     tableau{};
    Tableau     tableau2{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.fromString(tableau_string);

    std::string result_string = "0;0;0;1;0;\n0;0;1;0;0;\n";

    EXPECT_EQ(tableau.toString(), result_string);
}

TEST(TestTableau, DumpTableau) {
    Tableau     tableau{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.fromString(tableau_string);

    tableau.dump("tableau_dump.txt");
}

TEST(TestTableau, AccessValues) {
    Tableau     tableau{};
    std::string tableau_string = "Destabilizer = ['+IX', '+XI']";

    tableau.fromString(tableau_string);

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
    std::string tableau_string = "Destabilizer = ['+XI', '+IX']";

    tableau.fromString(tableau_string);

    tableau.clear();

    EXPECT_EQ(tableau.empty(), true);

    tableau.fromString(tableau_string);

    auto tableau2 = Tableau::getDiagonalTableau(2);

    EXPECT_EQ(tableau2[0][0], 0);
    EXPECT_EQ(tableau2[0][1], 0);
    EXPECT_EQ(tableau2[0][2], 1);
    EXPECT_EQ(tableau2[0][3], 0);
    EXPECT_EQ(tableau2[1][0], 0);
    EXPECT_EQ(tableau2[1][1], 0);
    EXPECT_EQ(tableau2[1][2], 0);
    EXPECT_EQ(tableau2[1][3], 1);

    EXPECT_EQ(tableau2, tableau);

    EXPECT_EQ(tableau.tableauDistance(tableau2, 2), 0);
}

TEST(TestTableau, LoadTableauFrom) {
    using namespace dd::literals;

    auto qc = qc::QuantumComputation(2U);
    qc.h(0);
    qc.x(1, 0_pc);

    Tableau     tableau{};
    Tableau     tableau1{};
    std::string tableau_string = "Stabilizer = ['+XI', '+IZ']";

    tableau.fromString(tableau_string);
    Tableau::generateTableau(tableau1, qc);

    EXPECT_EQ(tableau, tableau1);

    tableau1.clear();

    tableau1.import("examples/cliffordexamples/base-tableau.tabl");

    EXPECT_EQ(tableau, tableau1);

    tableau1.clear();

    qc.s(0);
    qc.x(1);
    qc.sdag(1);
    qc.z(1);
    qc.y(1);

    Tableau::generateTableau(tableau1, qc);

    auto compOP = std::make_unique<qc::CompoundOperation>(2);
    auto h0     = std::make_unique<qc::StandardOperation>(1, 0, qc::H);
    auto x1     = std::make_unique<qc::StandardOperation>(1, 0_pc, 1, qc::X);
    compOP->emplace_back(h0);
    compOP->emplace_back(x1);

    qc::QuantumComputation qc2(2U);
    qc2.emplace_back(compOP);

    Tableau::generateTableau(tableau, qc2);
}
TEST(TestTableau, InitTableau) {
    using namespace dd::literals;

    Tableau tableau{};
    tableau.init(4);

    EXPECT_EQ(tableau[0][0], 0);
    EXPECT_EQ(tableau[0][1], 0);
    EXPECT_EQ(tableau[0][2], 0);
    EXPECT_EQ(tableau[0][3], 0);
    EXPECT_EQ(tableau[0][4], 1);
    EXPECT_EQ(tableau[0][5], 0);
    EXPECT_EQ(tableau[0][6], 0);
    EXPECT_EQ(tableau[0][7], 0);
    EXPECT_EQ(tableau[0][8], 0);
}

TEST(TestTableau, BVTableau) {
    using namespace dd::literals;

    Tableau tableau{};
    tableau.init(2);
    unsigned long bitvector1 = 0b10;
    unsigned long bitvector2 = 0b01;

    tableau.populateTableauFrom(bitvector1, 2, 0);
    tableau.populateTableauFrom(bitvector2, 2, 1);

    EXPECT_EQ(tableau[0][0], 0);
    EXPECT_EQ(tableau[1][0], 1);
    EXPECT_EQ(tableau[0][1], 1);
    EXPECT_EQ(tableau[1][1], 0);

    unsigned long bitvector3 = tableau.getBVFrom(0);
    unsigned long bitvector4 = tableau.getBVFrom(1);

    EXPECT_EQ(bitvector3, bitvector1);
    EXPECT_EQ(bitvector4, bitvector2);
}
TEST(TestTableau, EmbedTableau) {
    using namespace dd::literals;

    Tableau tableau{};
    tableau.init(2);

    Tableau embededTableau = tableau.embedTableau(3);

    EXPECT_EQ(embededTableau, Tableau::getDiagonalTableau(3));
}
