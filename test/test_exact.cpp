/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "exact/ExactMapper.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

class ExactTest: public testing::TestWithParam<std::string> {
protected:
    std::string test_example_dir      = "./examples/";
    std::string test_architecture_dir = "./architectures/";
    std::string test_calibration_dir  = "./calibration/";

    qc::QuantumComputation qc{};
    Configuration          settings{};
    Architecture           IBMQ_Yorktown{};
    Architecture           IBMQ_London{};
    Architecture           IBM_QX4{};
    ExactMapper            IBMQ_Yorktown_mapper{qc, IBMQ_Yorktown};
    ExactMapper            IBMQ_London_mapper{qc, IBMQ_London};
    ExactMapper            IBM_QX4_mapper{qc, IBM_QX4};

    void SetUp() override {
        if (::testing::UnitTest::GetInstance()->current_test_info()->value_param()) {
            qc.import(test_example_dir + GetParam() + ".qasm");
        }
        IBMQ_Yorktown.loadCouplingMap(AvailableArchitecture::IBMQ_Yorktown);
        IBMQ_London.loadCouplingMap(test_architecture_dir + "ibmq_london.arch");
        IBMQ_London.loadCalibrationData(test_calibration_dir + "ibmq_london.csv");
        IBM_QX4.loadCouplingMap(AvailableArchitecture::IBM_QX4);
        settings.verbose = true;
    }
};

INSTANTIATE_TEST_SUITE_P(Exact, ExactTest,
                         testing::Values(
                                 "3_17_13",
                                 "ex-1_166",
                                 "ham3_102",
                                 "miller_11",
                                 "4gt11_84"),
                         [](const testing::TestParamInfo<ExactTest::ParamType>& info) {
		std::string name = info.param;
		std::replace(name.begin(), name.end(), '-', '_');
		std::stringstream ss{};
		ss << name;
		return ss.str(); });

TEST_P(ExactTest, IndividualGates) {
    settings.layering = Layering::IndividualGates;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_individual.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);

    IBMQ_London_mapper.map(settings);
    IBMQ_London_mapper.dumpResult(GetParam() + "_exact_london_individual.qasm");
    IBMQ_London_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, DisjointQubits) {
    settings.layering = Layering::DisjointQubits;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_disjoint.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);

    IBMQ_London_mapper.map(settings);
    IBMQ_London_mapper.dumpResult(GetParam() + "_exact_london_disjoint.qasm");
    IBMQ_London_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, OddGates) {
    settings.layering = Layering::OddGates;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_odd.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, QubitTriangle) {
    settings.layering = Layering::QubitTriangle;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_triangle.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, CommanderEncodingfixed3) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Fixed3;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_fixed3.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingfixed2) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Fixed2;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_fixed2.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodinghalves) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Halves;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_halves.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodinglogarithm) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Logarithm;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_log.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, CommanderEncodingUnidirectionalfixed3) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Fixed3;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_fixed3.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionalfixed2) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Fixed2;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_fixed2.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionalhalves) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Halves;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_halves.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionallogarithm) {
    settings.encoding          = Encoding::Commander;
    settings.commanderGrouping = CommanderGrouping::Logarithm;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_log.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, BimanderEncodingfixed3) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Fixed3;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingfixed2) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Fixed2;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodinghalves) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Halves;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodinglogaritm) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Logarithm;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, BimanderEncodingUnidirectionalfixed3) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Fixed3;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionalfixed2) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Fixed2;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionalhalves) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Halves;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionallogarithm) {
    settings.encoding          = Encoding::Bimander;
    settings.commanderGrouping = CommanderGrouping::Logarithm;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsBidirectional) {
    settings.enableSwapLimits = true;
    settings.useSubsets       = false;
    settings.swapReduction    = SwapReduction::CouplingLimit;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_swapreduct.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalSubsetSwaps) {
    settings.enableSwapLimits = true;
    settings.useSubsets       = true;
    settings.swapReduction    = SwapReduction::CouplingLimit;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_swapreduct.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalCustomLimit) {
    settings.enableSwapLimits = true;
    settings.swapReduction    = SwapReduction::Custom;
    settings.swapLimit        = 10;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_swapreduct.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsUnidirectional) {
    settings.enableSwapLimits = true;
    settings.useSubsets       = false;
    settings.swapReduction    = SwapReduction::CouplingLimit;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_swapreduct.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsUnidirectionalSubsetSwaps) {
    settings.enableSwapLimits = true;
    settings.useSubsets       = true;
    settings.swapReduction    = SwapReduction::CouplingLimit;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_swapreduct.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsUnidirectionalCustomLimit) {
    settings.enableSwapLimits = true;
    settings.swapReduction    = SwapReduction::Custom;
    settings.swapLimit        = 10;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_swapreduct.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, IncreasingCustomLimitUnidirectional) {
    settings.enableSwapLimits = true;
    settings.swapReduction    = SwapReduction::Increasing;
    settings.swapLimit        = 3;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_swapreduct_inccustom.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, IncreasingUnidirectional) {
    settings.enableSwapLimits = true;
    settings.swapReduction    = SwapReduction::Increasing;
    settings.swapLimit        = 0;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_swapreduct_inc.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsBidirectionalBDD) {
    settings.enableSwapLimits = true;
    settings.useSubsets       = false;
    settings.useBDD           = true;
    settings.swapReduction    = SwapReduction::CouplingLimit;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_swapreduct_bdd.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalSubsetSwapsBDD) {
    settings.enableSwapLimits = true;
    settings.useSubsets       = true;
    settings.useBDD           = true;
    settings.swapReduction    = SwapReduction::CouplingLimit;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_swapreduct_bdd.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, NoSubsets) {
    settings.useSubsets       = false;
    settings.enableSwapLimits = false;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_nosubsets.qasm");
    IBM_QX4_mapper.printResult(std::cout);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, toStringMethods) {
    EXPECT_EQ(toString(InitialLayout::Identity), "identity");
    EXPECT_EQ(toString(InitialLayout::Static), "static");
    EXPECT_EQ(toString(InitialLayout::Dynamic), "dynamic");
    EXPECT_EQ(toString(InitialLayout::None), "none");

    EXPECT_EQ(toString(Layering::IndividualGates), "individual_gates");
    EXPECT_EQ(toString(Layering::DisjointQubits), "disjoint_qubits");
    EXPECT_EQ(toString(Layering::OddGates), "odd_gates");
    EXPECT_EQ(toString(Layering::QubitTriangle), "qubit_triangle");
    EXPECT_EQ(toString(Layering::None), "none");

    EXPECT_EQ(toString(Encoding::Naive), "naive");
    EXPECT_EQ(toString(Encoding::Commander), "commander");
    EXPECT_EQ(toString(Encoding::Bimander), "bimander");

    EXPECT_EQ(toString(CommanderGrouping::Fixed2), "fixed2");
    EXPECT_EQ(toString(CommanderGrouping::Fixed3), "fixed3");
    EXPECT_EQ(toString(CommanderGrouping::Logarithm), "logarithm");
    EXPECT_EQ(toString(CommanderGrouping::Halves), "halves");

    EXPECT_EQ(toString(SwapReduction::CouplingLimit), "coupling_limit");
    EXPECT_EQ(toString(SwapReduction::Custom), "custom");
    EXPECT_EQ(toString(SwapReduction::None), "none");
    EXPECT_EQ(toString(SwapReduction::Increasing), "increasing");

    SUCCEED() << "ToStringMethods working";
}

TEST_F(ExactTest, CircuitWithOnlySingleQubitGates) {
    qc.addQubitRegister(2U);
    qc.x(0);
    qc.x(1);
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(std::cout, qc::OpenQASM);
    SUCCEED() << "Mapping successful";
}

TEST_F(ExactTest, MapToSubsetNotIncludingQ0) {
    using namespace dd::literals;

    CouplingMap  cm{{0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 3}, {3, 2}, {1, 3}, {3, 1}};
    Architecture arch(4U, cm);

    qc.addQubitRegister(3U);
    qc.x(0, 1_pc);
    qc.x(1, 2_pc);
    qc.x(2, 0_pc);

    auto mapper = ExactMapper(qc, arch);
    settings.useSubsets = false;
    mapper.map(settings);

    std::ostringstream oss{};
    mapper.dumpResult(oss, qc::OpenQASM);
    std::cout << oss.str() << std::endl;
    EXPECT_STREQ(oss.str().c_str(),
                 "// i 1 2 3 0\n"
                 "// o 1 2 3\n"
                 "OPENQASM 2.0;\n"
                 "include \"qelib1.inc\";\n"
                 "qreg q[4];\n"
                 "creg c[3];\n"
                 "cx q[2], q[1];\n"
                 "cx q[3], q[2];\n"
                 "cx q[1], q[3];\n"
                 "measure q[1] -> c[0];\n"
                 "measure q[2] -> c[1];\n"
                 "measure q[3] -> c[2];\n");
}
