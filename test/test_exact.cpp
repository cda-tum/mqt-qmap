/*
 * This file is part of the JKQ QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "exact/ExactMapper.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

class ExactTest: public testing::TestWithParam<std::string> {
protected:

	std::string test_example_dir = "./examples/";
	std::string test_architecture_dir = "./architectures/";
	std::string test_calibration_dir = "./calibration/";

	qc::QuantumComputation qc{};
	Architecture IBMQ_Yorktown{};
	Architecture IBMQ_London{};
	ExactMapper IBMQ_Yorktown_mapper{qc, IBMQ_Yorktown};
	ExactMapper IBMQ_London_mapper{qc, IBMQ_London};

	void SetUp() override {
		qc.import(test_example_dir + GetParam() + ".qasm");
		IBMQ_Yorktown.loadCouplingMap(AvailableArchitectures::IBMQ_Yorktown);
		IBMQ_London.loadCouplingMap(test_architecture_dir + "ibmq_london.arch");
		IBMQ_London.loadCalibrationData(test_calibration_dir + "ibmq_london.csv");
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

