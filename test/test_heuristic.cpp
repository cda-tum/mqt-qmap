/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "heuristic/HeuristicMapper.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

class HeuristicTest5Q: public testing::TestWithParam<std::string> {
protected:
	std::string test_example_dir = "../../examples/";
	std::string test_architecture_dir = "../../extern/architectures/";
	std::string test_calibration_dir = "../../extern/calibration/";
};

INSTANTIATE_TEST_SUITE_P(Heuristic, HeuristicTest5Q,
                         testing::Values(
		                         "3_17_13",
		                         "ex-1_166",
		                         "ham3_102",
		                         "miller_11",
		                         "4gt11_84",
		                         "4mod5-v0_20",
		                         "mod5d1_63"
                                        ),
                         [](const testing::TestParamInfo<HeuristicTest5Q::ParamType>& info) {
	                         std::string name = info.param;
	                         std::replace( name.begin(), name.end(), '-', '_');
	                         std::stringstream ss{};
	                         ss << name;
	                         return ss.str();});

TEST_P(HeuristicTest5Q, Identity) {
	auto IBM_QX4_mapper = HeuristicMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx4.arch");
	auto IBMQ_London_mapper = HeuristicMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibmq_london.arch", test_calibration_dir + "ibmq_london.csv");
	MappingSettings settings{};
	settings.initialLayoutStrategy = InitialLayoutStrategy::Identity;
	IBM_QX4_mapper.map(settings);
	IBM_QX4_mapper.dumpResult(GetParam() + "_heuristic_qx4_identity.qasm");
	IBM_QX4_mapper.printResult(std::cout, true);

	IBMQ_London_mapper.map(settings);
	IBMQ_London_mapper.dumpResult(GetParam() + "_heuristic_london_identity.qasm");
	IBMQ_London_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Static) {
	auto IBM_QX4_mapper = HeuristicMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx4.arch");
	auto IBMQ_London_mapper = HeuristicMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibmq_london.arch", test_calibration_dir + "ibmq_london.csv");
	MappingSettings settings{};
	settings.initialLayoutStrategy = InitialLayoutStrategy::Static;
	IBM_QX4_mapper.map(settings);
	IBM_QX4_mapper.dumpResult(GetParam() + "_heuristic_qx4_static.qasm");
	IBM_QX4_mapper.printResult(std::cout, true);
	IBMQ_London_mapper.map(settings);
	IBMQ_London_mapper.dumpResult(GetParam() + "_heuristic_london_static.qasm");
	IBMQ_London_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}

TEST_P(HeuristicTest5Q, Dynamic) {
	auto IBM_QX4_mapper = HeuristicMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx4.arch");
	auto IBMQ_London_mapper = HeuristicMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibmq_london.arch", test_calibration_dir + "ibmq_london.csv");
	MappingSettings settings{};
	settings.initialLayoutStrategy = InitialLayoutStrategy::Dynamic;
	IBM_QX4_mapper.map(settings);
	IBM_QX4_mapper.dumpResult(GetParam() + "_heuristic_qx4_dynamic.qasm");
	IBM_QX4_mapper.printResult(std::cout, true);
	IBMQ_London_mapper.map(settings);
	IBMQ_London_mapper.dumpResult(GetParam() + "_heuristic_london_dynamic.qasm");
	IBMQ_London_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}

class HeuristicTest16Q: public testing::TestWithParam<std::string> {
protected:

	std::string test_example_dir = "../../examples/";
	std::string test_architecture_dir = "../../extern/architectures/";

};

INSTANTIATE_TEST_SUITE_P(Heuristic, HeuristicTest16Q,
                         testing::Values(
		                         "ising_model_10",
		                         "rd73_140",
		                         "cnt3-5_179",
		                         "qft_16"
                                        ),
                         [](const testing::TestParamInfo<HeuristicTest16Q::ParamType>& info) {
	                         std::string name = info.param;
	                         std::replace( name.begin(), name.end(), '-', '_');
	                         std::stringstream ss{};
	                         ss << name;
	                         return ss.str();});

TEST_P(HeuristicTest16Q, Dynamic) {
	auto IBM_QX5_mapper = HeuristicMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx5.arch");
	MappingSettings settings{};
	settings.initialLayoutStrategy = InitialLayoutStrategy::Dynamic;
	IBM_QX5_mapper.map(settings);
	IBM_QX5_mapper.dumpResult(GetParam() + "_heuristic_qx5_dynamic.qasm");
	IBM_QX5_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}
