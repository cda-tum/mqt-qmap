/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
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
    Architecture           IBMQ_Yorktown{};
    Architecture           IBMQ_London{};
    Architecture           IBM_QX4{};
    ExactMapper            IBMQ_Yorktown_mapper{qc, IBMQ_Yorktown};
    ExactMapper            IBMQ_London_mapper{qc, IBMQ_London};
    ExactMapper            IBM_QX4_mapper{qc, IBM_QX4};

    void SetUp() override {
        qc.import(test_example_dir + GetParam() + ".qasm");
        IBMQ_Yorktown.loadCouplingMap(AvailableArchitectures::IBMQ_Yorktown);
        IBMQ_London.loadCouplingMap(test_architecture_dir + "ibmq_london.arch");
        IBMQ_London.loadCalibrationData(test_calibration_dir + "ibmq_london.csv");
        IBM_QX4.loadCouplingMap(AvailableArchitectures::IBM_QX4);
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
    MappingSettings settings{};
    settings.layeringStrategy = LayeringStrategy::IndividualGates;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_individual.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);

    IBMQ_London_mapper.map(settings);
    IBMQ_London_mapper.dumpResult(GetParam() + "_exact_london_individual.qasm");
    IBMQ_London_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, DisjointQubits) {
    MappingSettings settings{};
    settings.layeringStrategy = LayeringStrategy::DisjointQubits;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_disjoint.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);

    IBMQ_London_mapper.map(settings);
    IBMQ_London_mapper.dumpResult(GetParam() + "_exact_london_disjoint.qasm");
    IBMQ_London_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, OddGates) {
    MappingSettings settings{};
    settings.layeringStrategy = LayeringStrategy::OddGates;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_odd.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, QubitTriangle) {
    MappingSettings settings{};
    settings.layeringStrategy = LayeringStrategy::QubitTriangle;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_triangle.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, CommanderEncodingfixed3) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Fixed3;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_fixed3.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingfixed2) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Fixed2;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_fixed2.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodinghalves) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Halves;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_halves.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodinglogarithm) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Logarithm;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_commander_log.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, CommanderEncodingUnidirectionalfixed3) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Fixed3;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_fixed3.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionalfixed2) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Fixed2;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_fixed2.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionalhalves) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Halves;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_halves.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, CommanderEncodingUnidirectionallogarithm) {
    MappingSettings settings{};
    settings.encoding = Encodings::Commander;
    settings.grouping = CMDRVariableGroupings::Logarithm;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_commander_log.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, BimanderEncodingfixed3) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Fixed3;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingfixed2) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Fixed2;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodinghalves) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Halves;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodinglogaritm) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Logarithm;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bimander.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, BimanderEncodingUnidirectionalfixed3) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Fixed3;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionalfixed2) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Fixed2;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionalhalves) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Halves;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, BimanderEncodingUnidirectionallogarithm) {
    MappingSettings settings{};
    settings.encoding = Encodings::Bimander;
    settings.grouping = CMDRVariableGroupings::Logarithm;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bimander.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsBidirectional) {
    MappingSettings settings{};
    settings.enableLimits    = true;
    settings.useQubitSubsets = false;
    settings.strategy        = SwapReductionStrategy::CouplingLimit;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bdd.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalSubsetSwaps) {
    MappingSettings settings{};
    settings.enableLimits    = true;
    settings.useQubitSubsets = true;
    settings.strategy        = SwapReductionStrategy::CouplingLimit;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bdd.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsBidirectionalCustomLimit) {
    MappingSettings settings{};
    settings.enableLimits = true;
    settings.strategy     = SwapReductionStrategy::Custom;
    settings.limit        = 10;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_exact_yorktown_bdd.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, LimitsUnidirectional) {
    MappingSettings settings{};
    settings.enableLimits    = true;
    settings.useQubitSubsets = false;
    settings.strategy        = SwapReductionStrategy::CouplingLimit;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bdd.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsUnidirectionalSubsetSwaps) {
    MappingSettings settings{};
    settings.enableLimits    = true;
    settings.useQubitSubsets = true;
    settings.strategy        = SwapReductionStrategy::CouplingLimit;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bdd.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, LimitsUnidirectionalCustomLimit) {
    MappingSettings settings{};
    settings.enableLimits = true;
    settings.strategy     = SwapReductionStrategy::Custom;
    settings.limit        = 10;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bdd.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, IncreasingCustomLimitUnidirectional) {
    MappingSettings settings{};
    settings.enableLimits = true;
    settings.strategy     = SwapReductionStrategy::Increasing;
    settings.limit        = 3;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bdd_inccustom.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, IncreasingUnidirectional) {
    MappingSettings settings{};
    settings.enableLimits = true;
    settings.strategy     = SwapReductionStrategy::Increasing;
    settings.limit        = 0;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_bdd_inc.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, NoSubsets) {
    MappingSettings settings{};
    settings.useQubitSubsets = false;
    settings.enableLimits    = false;
    IBM_QX4_mapper.map(settings);
    IBM_QX4_mapper.dumpResult(GetParam() + "_exact_QX4_nosubsets.qasm");
    IBM_QX4_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
TEST_P(ExactTest, toStringMethods) {
    EXPECT_EQ(toString(InitialLayoutStrategy::Identity), "identity");
    EXPECT_EQ(toString(InitialLayoutStrategy::Static), "static");
    EXPECT_EQ(toString(InitialLayoutStrategy::Dynamic), "dynamic");
    EXPECT_EQ(toString(InitialLayoutStrategy::None), "none");

    EXPECT_EQ(toString(LayeringStrategy::IndividualGates), "individual_gates");
    EXPECT_EQ(toString(LayeringStrategy::DisjointQubits), "disjoint_qubits");
    EXPECT_EQ(toString(LayeringStrategy::OddGates), "odd_gates");
    EXPECT_EQ(toString(LayeringStrategy::QubitTriangle), "qubit_triangle");
    EXPECT_EQ(toString(LayeringStrategy::None), "none");

    EXPECT_EQ(toString(Encodings::Naive), "naive");
    EXPECT_EQ(toString(Encodings::Commander), "commander");
    EXPECT_EQ(toString(Encodings::Bimander), "bimander");

    EXPECT_EQ(toString(CMDRVariableGroupings::Fixed2), "fixed2");
    EXPECT_EQ(toString(CMDRVariableGroupings::Fixed3), "fixed3");
    EXPECT_EQ(toString(CMDRVariableGroupings::Logarithm), "logarithm");
    EXPECT_EQ(toString(CMDRVariableGroupings::Halves), "halves");

    EXPECT_EQ(toString(SwapReductionStrategy::CouplingLimit), "coupling_limit");
    EXPECT_EQ(toString(SwapReductionStrategy::Custom), "custom");
    EXPECT_EQ(toString(SwapReductionStrategy::None), "none");
    EXPECT_EQ(toString(SwapReductionStrategy::Increasing), "increasing");

    SUCCEED() << "ToStringMethods working";
}
