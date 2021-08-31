/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "heuristic/HeuristicMapper.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

class HeuristicTest5Q: public testing::TestWithParam<std::string> {
protected:
    std::string test_example_dir      = "./examples/";
    std::string test_architecture_dir = "./architectures/";
    std::string test_calibration_dir  = "./calibration/";

    qc::QuantumComputation qc{};
    Architecture           IBMQ_Yorktown{};
    Architecture           IBMQ_London{};
    HeuristicMapper        IBMQ_Yorktown_mapper{qc, IBMQ_Yorktown};
    HeuristicMapper        IBMQ_London_mapper{qc, IBMQ_London};

    void SetUp() override {
        qc.import(test_example_dir + GetParam() + ".qasm");
        IBMQ_Yorktown.loadCouplingMap(AvailableArchitectures::IBMQ_Yorktown);
        IBMQ_London.loadCouplingMap(test_architecture_dir + "ibmq_london.arch");
        IBMQ_London.loadCalibrationData(test_calibration_dir + "ibmq_london.csv");
    }
};

INSTANTIATE_TEST_SUITE_P(Heuristic, HeuristicTest5Q,
                         testing::Values(
                                 "3_17_13",
                                 "ex-1_166",
                                 "ham3_102",
                                 "miller_11",
                                 "4gt11_84",
                                 "4mod5-v0_20",
                                 "mod5d1_63"),
                         [](const testing::TestParamInfo<HeuristicTest5Q::ParamType>& info) {
	                         std::string name = info.param;
	                         std::replace( name.begin(), name.end(), '-', '_');
	                         std::stringstream ss{};
	                         ss << name;
	                         return ss.str(); });

TEST_P(HeuristicTest5Q, Identity) {
    MappingSettings settings{};
    settings.initialLayoutStrategy = InitialLayoutStrategy::Identity;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_heuristic_qx4_identity.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);

    IBMQ_London_mapper.map(settings);
    IBMQ_London_mapper.dumpResult(GetParam() + "_heuristic_london_identity.qasm");
    IBMQ_London_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Static) {
    MappingSettings settings{};
    settings.initialLayoutStrategy = InitialLayoutStrategy::Static;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_heuristic_qx4_static.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    IBMQ_London_mapper.map(settings);
    IBMQ_London_mapper.dumpResult(GetParam() + "_heuristic_london_static.qasm");
    IBMQ_London_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Dynamic) {
    MappingSettings settings{};
    settings.initialLayoutStrategy = InitialLayoutStrategy::Dynamic;
    IBMQ_Yorktown_mapper.map(settings);
    IBMQ_Yorktown_mapper.dumpResult(GetParam() + "_heuristic_qx4_dynamic.qasm");
    IBMQ_Yorktown_mapper.printResult(std::cout, true);
    IBMQ_London_mapper.map(settings);
    IBMQ_London_mapper.dumpResult(GetParam() + "_heuristic_london_dynamic.qasm");
    IBMQ_London_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

class HeuristicTest16Q: public testing::TestWithParam<std::string> {
protected:
    std::string test_example_dir      = "../../examples/";
    std::string test_architecture_dir = "../../extern/architectures/";

    qc::QuantumComputation qc{};
    Architecture           IBM_QX5{};
    HeuristicMapper        IBM_QX5_mapper{qc, IBM_QX5};

    void SetUp() override {
        qc.import(test_example_dir + GetParam() + ".qasm");
        IBM_QX5.loadCouplingMap(test_architecture_dir + "ibm_qx5.arch");
    }
};

INSTANTIATE_TEST_SUITE_P(Heuristic, HeuristicTest16Q,
                         testing::Values(
                                 "ising_model_10",
                                 "rd73_140",
                                 "cnt3-5_179",
                                 "qft_16"),
                         [](const testing::TestParamInfo<HeuristicTest16Q::ParamType>& info) {
	                         std::string name = info.param;
	                         std::replace( name.begin(), name.end(), '-', '_');
	                         std::stringstream ss{};
	                         ss << name;
	                         return ss.str(); });

TEST_P(HeuristicTest16Q, Dynamic) {
    MappingSettings settings{};
    settings.initialLayoutStrategy = InitialLayoutStrategy::Dynamic;
    IBM_QX5_mapper.map(settings);
    IBM_QX5_mapper.dumpResult(GetParam() + "_heuristic_qx5_dynamic.qasm");
    IBM_QX5_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

class HeuristicTest20Q: public testing::TestWithParam<std::string> {
protected:
    std::string test_example_dir      = "../../examples/";
    std::string test_architecture_dir = "../../extern/architectures/";

    qc::QuantumComputation qc{};
    Architecture           arch{};
    HeuristicMapper        tokyo_mapper{qc, arch};

    void SetUp() override {
        qc.import(test_example_dir + GetParam() + ".qasm");
        arch.loadCouplingMap(test_architecture_dir + "ibmq_tokyo_20qubit.arch");
    }
};

INSTANTIATE_TEST_SUITE_P(Heuristic, HeuristicTest20Q,
                         testing::Values(
                                 "ising_model_10",
                                 "rd73_140",
                                 "cnt3-5_179",
                                 "qft_16",
                                 "z4_268"),
                         [](const testing::TestParamInfo<HeuristicTest20Q::ParamType>& info) {
                             std::string name = info.param;
                             std::replace( name.begin(), name.end(), '-', '_');
                             std::stringstream ss{};
                             ss << name;
                             return ss.str(); });

TEST_P(HeuristicTest20Q, Dynamic) {
    MappingSettings settings{};
    settings.initialLayoutStrategy = InitialLayoutStrategy::Dynamic;
    tokyo_mapper.map(settings);
    tokyo_mapper.dumpResult(GetParam() + "_heuristic_tokyo_dynamic.qasm");
    tokyo_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}

class HeuristicTest20QTeleport: public testing::TestWithParam<std::tuple<unsigned long long, std::string>> {
protected:
    std::string test_example_dir      = "../../examples/";
    std::string test_architecture_dir = "../../extern/architectures/";

    qc::QuantumComputation qc{};
    Architecture           arch{};
    HeuristicMapper        tokyo_mapper{qc, arch};

    void SetUp() override {
        qc.import(test_example_dir + std::get<1>(GetParam()) + ".qasm");
        arch.loadCouplingMap(test_architecture_dir + "ibmq_tokyo_20qubit.arch");
    }
};

INSTANTIATE_TEST_SUITE_P(HeuristicTeleport, HeuristicTest20QTeleport,
                         testing::Combine(
                                 testing::Values(1, 2, 3, 1337, 1338, 3147),
                                 testing::Values("ising_model_10", "rd73_140", "cnt3-5_179", "qft_16", "z4_268")),
                         [](const testing::TestParamInfo<HeuristicTest20QTeleport::ParamType>& info) {
                             std::string name = std::get<1>(info.param);
                             std::replace( name.begin(), name.end(), '-', '_');
                             std::stringstream ss{};
                             ss << name << "_seed" << std::get<0>(info.param);
                             return ss.str(); });

TEST_P(HeuristicTest20QTeleport, Teleportation) {
    MappingSettings settings{};
    settings.initialLayoutStrategy = InitialLayoutStrategy::Dynamic;
    settings.teleportationQubits   = std::min((arch.getNqubits() - qc.getNqubits()) & ~1u, 8u);
    settings.teleportationSeed     = std::get<0>(GetParam());
    tokyo_mapper.map(settings);
    tokyo_mapper.dumpResult(std::get<1>(GetParam()) + "_heuristic_tokyo_teleport.qasm");
    tokyo_mapper.printResult(std::cout, true);
    SUCCEED() << "Mapping successful";
}
