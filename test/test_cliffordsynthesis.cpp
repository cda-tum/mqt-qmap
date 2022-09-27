#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

using namespace cs;

class TestCliffordSynthesis: public testing::TestWithParam<std::string> {
protected:
    std::string testArchitectureDir = "./architectures/";
    std::string testCalibrationDir  = "./calibration/";
    std::string testExampleDir      = "./examples/cliffordexamples/";

    Architecture ibmqYorktown{};
    Architecture ibmqLondon{};
    Architecture ibmQX4{};

    std::unique_ptr<CliffordSynthesizer> qx4Optimizer{};
    std::unique_ptr<CliffordSynthesizer> yorktownOptimizer{};
    std::unique_ptr<CliffordSynthesizer> londonOptimizer{};

    void SetUp() override {
        using namespace dd::literals;
        util::init();
        ibmqYorktown.loadCouplingMap(AvailableArchitecture::IBMQ_Yorktown);
        ibmqLondon.loadCouplingMap(testArchitectureDir + "ibmq_london.arch");
        ibmqLondon.loadProperties(testCalibrationDir + "ibmq_london.csv");
        ibmQX4.loadCouplingMap(AvailableArchitecture::IBM_QX4);

        yorktownOptimizer = std::make_unique<CliffordSynthesizer>();

        londonOptimizer = std::make_unique<CliffordSynthesizer>();

        qx4Optimizer = std::make_unique<CliffordSynthesizer>();
    }
};

INSTANTIATE_TEST_SUITE_P(
        CliffordSynthesizer, TestCliffordSynthesis,
        testing::Values(
                "destabilizer.txt"));

TEST_P(TestCliffordSynthesis, SimpleSynthesis) {
    const auto& inputFile = GetParam();
    Tableau     tableau{};
    if (inputFile.find(".txt") != std::string::npos) {
        auto is = std::ifstream(testExampleDir + inputFile);
        if (!is.good()) {
            FATAL() << "Error opening file " << inputFile;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.fromString(line);
            Configuration configuration{};

            configuration.nqubits          = 2;
            configuration.initialTimesteps = 10;
            Tableau::initTableau(configuration.initialTableau, 2);
            configuration.targetTableau = tableau;
            qx4Optimizer->synthesize(configuration);
            tableau.clear();

            qx4Optimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(qx4Optimizer->optimalResults.result, logicbase::Result::SAT);
        }
    }
}

TEST(TestCliffordSynthesis, SanityCheck) {
    util::init();
    qc::QuantumComputation qc{};
    CliffordSynthesizer    cs{};
    Configuration          configuration{false, false, 2, 10, OptimizationStrategy::UseMinimizer, TargetMetric::DEPTH};
    configuration.verbosity = 5;
    qc.addQubitRegister(2U);
    qc.h(0);
    qc.h(0);
    qc.h(0);
    qc.h(0);
    qc.h(0);

    configuration.targetCircuit = qc.clone();

    cs.optimize(configuration);

    EXPECT_EQ(cs.optimalResults.depth, 1);
}

TEST_P(TestCliffordSynthesis, TestDepthOpt) {
    const auto& inputFile = GetParam();
    Tableau     tableau{};
    if (inputFile.find(".txt") != std::string::npos) {
        auto is = std::ifstream(testExampleDir + inputFile);
        if (!is.good()) {
            FATAL() << "Error opening file " << inputFile;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.fromString(line);
            Configuration configuration{};

            configuration.nqubits          = 2;
            configuration.initialTimesteps = 10;
            configuration.target           = TargetMetric::DEPTH;
            Tableau::initTableau(configuration.initialTableau, 2);
            configuration.targetTableau = tableau;
            qx4Optimizer->synthesize(configuration);
            tableau.clear();

            qx4Optimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(qx4Optimizer->optimalResults.result, logicbase::Result::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestFidelityOpt) {
    const auto& inputFile = GetParam();
    Tableau     tableau{};
    if (inputFile.find(".txt") != std::string::npos) {
        auto is = std::ifstream(testExampleDir + inputFile);
        if (!is.good()) {
            FATAL() << "Error opening file " << inputFile;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.fromString(line);
            Configuration configuration{};

            configuration.nqubits          = 2;
            configuration.initialTimesteps = 10;
            configuration.target           = TargetMetric::FIDELITY;
            configuration.architecture     = ibmqLondon;
            Tableau::initTableau(configuration.initialTableau, 2);
            configuration.targetTableau = tableau;
            londonOptimizer->synthesize(configuration);
            tableau.clear();

            londonOptimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(londonOptimizer->optimalResults.result, logicbase::Result::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestCNOTONLYOpt) {
    const auto& inputFile = GetParam();
    Tableau     tableau{};
    if (inputFile.find(".txt") != std::string::npos) {
        auto is = std::ifstream(testExampleDir + inputFile);
        if (!is.good()) {
            FATAL() << "Error opening file " << inputFile;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.fromString(line);
            Configuration configuration{};

            configuration.nqubits          = 2;
            configuration.initialTimesteps = 10;
            configuration.target           = TargetMetric::GATES_ONLY_CNOT;
            Tableau::initTableau(configuration.initialTableau, 2);
            configuration.targetTableau = tableau;
            qx4Optimizer->synthesize(configuration);
            tableau.clear();

            qx4Optimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(qx4Optimizer->optimalResults.result, logicbase::Result::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestStartLow) {
    const auto& inputFile = GetParam();
    Tableau     tableau{};
    if (inputFile.find(".txt") != std::string::npos) {
        auto is = std::ifstream(testExampleDir + inputFile);
        if (!is.good()) {
            FATAL() << "Error opening file " << inputFile;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.fromString(line);
            Configuration configuration{};

            configuration.nqubits          = 2;
            configuration.initialTimesteps = 10;
            configuration.target           = TargetMetric::GATES;
            configuration.strategy         = OptimizationStrategy::StartLow;
            Tableau::initTableau(configuration.initialTableau, 2);
            configuration.targetTableau = tableau;
            qx4Optimizer->synthesize(configuration);
            tableau.clear();

            qx4Optimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(qx4Optimizer->optimalResults.result, logicbase::Result::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestStartHigh) {
    const auto& inputFile = GetParam();
    Tableau     tableau{};
    if (inputFile.find(".txt") != std::string::npos) {
        auto is = std::ifstream(testExampleDir + inputFile);
        if (!is.good()) {
            FATAL() << "Error opening file " << inputFile;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.fromString(line);
            Configuration configuration{};

            configuration.nqubits          = 2;
            configuration.initialTimesteps = 50;
            configuration.target           = TargetMetric::GATES;
            configuration.strategy         = OptimizationStrategy::StartHigh;
            Tableau::initTableau(configuration.initialTableau, 2);
            configuration.targetTableau = tableau;
            qx4Optimizer->synthesize(configuration);
            tableau.clear();

            qx4Optimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(qx4Optimizer->optimalResults.result, logicbase::Result::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestMinMax) {
    const auto& inputFile = GetParam();
    Tableau     tableau{};
    if (inputFile.find(".txt") != std::string::npos) {
        auto is = std::ifstream(testExampleDir + inputFile);
        if (!is.good()) {
            FATAL() << "Error opening file " << inputFile;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.fromString(line);
            Configuration configuration{};

            configuration.nqubits          = 2;
            configuration.initialTimesteps = 10;
            configuration.target           = TargetMetric::GATES;
            configuration.strategy         = OptimizationStrategy::MinMax;
            Tableau::initTableau(configuration.initialTableau, 2);
            configuration.targetTableau = tableau;
            qx4Optimizer->synthesize(configuration);
            tableau.clear();

            qx4Optimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(qx4Optimizer->optimalResults.result, logicbase::Result::SAT);
        }
    }
}

TEST(TestCliffordSynthesis, TestSplitIter) {
    util::init();
    using namespace dd::literals;
    qc::QuantumComputation qc{};
    CliffordSynthesizer    optimizer{};
    qc.addQubitRegister(2U);

    qc.x(1);
    qc.z(0);
    qc.y(0);
    qc.x(1, 0_pc);
    qc.h(1);
    qc.s(0);
    qc.x(1);
    qc.sdag(1);
    qc.z(0);
    qc.y(1);
    qc.x(1);
    qc.z(0);
    qc.y(0);
    qc.x(1, 0_pc);
    qc.h(1);
    qc.s(0);
    qc.x(1);
    qc.sdag(1);
    qc.z(0);
    qc.y(1);
    qc.x(1);
    qc.z(0);
    qc.y(0);
    qc.x(1, 0_pc);
    qc.h(1);
    qc.s(0);
    qc.x(1);
    qc.sdag(1);
    qc.z(0);
    qc.y(1);
    qc.x(1);
    qc.z(0);
    qc.y(0);

    Configuration configuration{};
    configuration.nqubits          = 2;
    configuration.initialTimesteps = 20;
    configuration.nThreads         = 1;
    configuration.targetCircuit    = qc.clone();
    configuration.target           = TargetMetric::GATES;
    configuration.strategy         = OptimizationStrategy::SplitIter;

    optimizer.synthesize(configuration);
    EXPECT_EQ(optimizer.optimalResults.result, logicbase::Result::SAT);
}
