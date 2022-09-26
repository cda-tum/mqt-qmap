#include "depthsynthesis/DepthSynthesizer.hpp"
#include "fidelitysynthesis/FidelitySynthesizer.hpp"
#include "gatesynthesis/GateSynthesizer.hpp"

#include "gtest/gtest.h"

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

        yorktownOptimizer = std::make_unique<GateSynthesizer>();
        yorktownOptimizer->setArchitecture(ibmqYorktown);

        londonOptimizer = std::make_unique<CliffordSynthesizer>();
        londonOptimizer->setArchitecture(ibmqLondon);

        qx4Optimizer = std::make_unique<CliffordSynthesizer>();
        qx4Optimizer->setArchitecture(ibmQX4);
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
            qx4Optimizer->nqubits          = 2;
            qx4Optimizer->initialTimesteps = 10;
            Tableau::initTableau(qx4Optimizer->initialTableau, 2);
            qx4Optimizer->targetTableau = tableau;
            qx4Optimizer->optimize();
            tableau.clear();

            londonOptimizer->optimalResults.dump(std::cout);

            qx4Optimizer->optimalResults.dump(std::cout);
            EXPECT_EQ(qx4Optimizer->optimalResults.result, SynthesisResult::SAT);
        }
    }
}

TEST(TestCliffordSynthesis, SanityCheck) {
    util::init();
    qc::QuantumComputation qc{};
    CliffordSynthesizer    optimizer{false, false, 2, 10, SynthesisStrategy::UseMinimizer, SynthesisTarget::DEPTH};
    qc.addQubitRegister(2U);
    qc.h(0);
    qc.h(0);
    qc.h(0);
    qc.h(0);
    qc.h(0);

    optimizer.setCircuit(qc);

    optimizer.optimize();

    EXPECT_EQ(optimizer.optimalResults.depth, 1);
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
            qx4Optimizer->nqubits          = 2;
            qx4Optimizer->initialTimesteps = 10;
            qx4Optimizer->target           = SynthesisTarget::DEPTH;
            Tableau::initTableau(qx4Optimizer->initialTableau, 2);
            qx4Optimizer->targetTableau = tableau;
            qx4Optimizer->optimize();
            tableau.clear();

            londonOptimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(qx4Optimizer->optimalResults.result, SynthesisResult::SAT);
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
            londonOptimizer->nqubits          = 2;
            londonOptimizer->initialTimesteps = 5;
            londonOptimizer->target           = SynthesisTarget::FIDELITY;
            Tableau::initTableau(londonOptimizer->initialTableau, 2);
            londonOptimizer->targetTableau = tableau;
            londonOptimizer->optimize();
            tableau.clear();

            londonOptimizer->optimalResults.dump(std::cout);

            EXPECT_EQ(londonOptimizer->optimalResults.result, SynthesisResult::SAT);
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
            qx4Optimizer->nqubits          = 2;
            qx4Optimizer->initialTimesteps = 10;
            qx4Optimizer->target           = SynthesisTarget::GATES_ONLY_CNOT;
            Tableau::initTableau(qx4Optimizer->initialTableau, 2);
            qx4Optimizer->targetTableau = tableau;
            qx4Optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4Optimizer->optimalResults.result, SynthesisResult::SAT);
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
            qx4Optimizer->nqubits          = 2;
            qx4Optimizer->initialTimesteps = 10;
            qx4Optimizer->target           = SynthesisTarget::GATES;
            qx4Optimizer->strategy         = SynthesisStrategy::StartLow;
            Tableau::initTableau(qx4Optimizer->initialTableau, 2);
            qx4Optimizer->targetTableau = tableau;
            qx4Optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4Optimizer->optimalResults.result, SynthesisResult::SAT);
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
            qx4Optimizer->nqubits          = 2;
            qx4Optimizer->initialTimesteps = 25;
            qx4Optimizer->target           = SynthesisTarget::GATES;
            qx4Optimizer->strategy         = SynthesisStrategy::StartHigh;
            Tableau::initTableau(qx4Optimizer->initialTableau, 2);
            qx4Optimizer->targetTableau = tableau;
            qx4Optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4Optimizer->optimalResults.result, SynthesisResult::SAT);
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
            qx4Optimizer->nqubits          = 2;
            qx4Optimizer->initialTimesteps = 10;
            qx4Optimizer->target           = SynthesisTarget::GATES;
            qx4Optimizer->strategy         = SynthesisStrategy::MinMax;
            Tableau::initTableau(qx4Optimizer->initialTableau, 2);
            qx4Optimizer->targetTableau = tableau;
            qx4Optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4Optimizer->optimalResults.result, SynthesisResult::SAT);
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

    optimizer.nqubits          = 2;
    optimizer.initialTimesteps = 20;
    optimizer.nthreads         = 1;
    optimizer.circuit          = qc.clone();
    optimizer.target           = SynthesisTarget::GATES;
    optimizer.strategy         = SynthesisStrategy::SplitIter;

    optimizer.optimize();
    EXPECT_EQ(optimizer.optimalResults.result, SynthesisResult::SAT);
}
