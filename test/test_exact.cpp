/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "exact/ExactMapper.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

class ExactTest: public testing::TestWithParam<std::string> {
protected:

	std::string test_example_dir = "../../examples/";
	std::string test_architecture_dir = "../../extern/architectures/";
	std::string test_calibration_dir = "../../extern/calibration/";

	void SetUp() override {

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
	                         std::replace( name.begin(), name.end(), '-', '_');
	                         std::stringstream ss{};
	                         ss << name;
	                         return ss.str();});

TEST_P(ExactTest, IndividualGates) {
	auto IBM_QX4_mapper = ExactMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx4.arch");
	auto IBMQ_London_mapper = ExactMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibmq_london.arch", test_calibration_dir + "ibmq_london.csv");

	MappingSettings settings{};
	settings.layeringStrategy = LayeringStrategy::IndividualGates;
	IBM_QX4_mapper.map(settings);
	IBM_QX4_mapper.dumpResult(GetParam() + "_exact_qx4_individual.qasm");
	IBM_QX4_mapper.printResult(std::cout, true);

	IBMQ_London_mapper.map(settings);
	IBMQ_London_mapper.dumpResult(GetParam() + "_exact_london_individual.qasm");
	IBMQ_London_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, DisjointQubits) {
	auto IBM_QX4_mapper = ExactMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx4.arch");
	auto IBMQ_London_mapper = ExactMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibmq_london.arch", test_calibration_dir + "ibmq_london.csv");

	MappingSettings settings{};
	settings.layeringStrategy = LayeringStrategy::DisjointQubits;
	IBM_QX4_mapper.map(settings);
	IBM_QX4_mapper.dumpResult(GetParam() + "_exact_qx4_disjoint.qasm");
	IBM_QX4_mapper.printResult(std::cout, true);

	IBMQ_London_mapper.map(settings);
	IBMQ_London_mapper.dumpResult(GetParam() + "_exact_london_disjoint.qasm");
	IBMQ_London_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, OddGates) {
	auto IBM_QX4_mapper = ExactMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx4.arch");
	MappingSettings settings{};
	settings.layeringStrategy = LayeringStrategy::OddGates;
	IBM_QX4_mapper.map(settings);
	IBM_QX4_mapper.dumpResult(GetParam() + "_exact_qx4_odd.qasm");
	IBM_QX4_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}

TEST_P(ExactTest, QubitTriangle) {
	auto IBM_QX4_mapper = ExactMapper(test_example_dir + GetParam() + ".qasm", test_architecture_dir + "ibm_qx4.arch");
	MappingSettings settings{};
	settings.layeringStrategy = LayeringStrategy::QubitTriangle;
	IBM_QX4_mapper.map(settings);
	IBM_QX4_mapper.dumpResult(GetParam() + "_exact_qx4_triangle.qasm");
	IBM_QX4_mapper.printResult(std::cout, true);
	SUCCEED() << "Mapping successful";
}

