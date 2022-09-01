#include "Architecture.hpp"
#include "Tableau.hpp"
#include "cliffordsynthesis/CliffordSynthesizer.hpp"

#include "gtest/gtest.h"

class TestCliffordSynthesis: public testing::TestWithParam<std::string> {
protected:
    std::string test_architecture_dir = "./architectures/";
    std::string test_calibration_dir  = "./calibration/";
    std::string test_example_dir      = "./examples/cliffordexamples/";

    Architecture IBMQ_Yorktown{};
    Architecture IBMQ_London{};
    Architecture IBM_QX4{};

    std::unique_ptr<CliffordOptimizer> qx4_optimizer{};
    std::unique_ptr<CliffordOptimizer> yorktown_optimizer{};
    std::unique_ptr<CliffordOptimizer> london_optimizer{};

    void SetUp() override {
        using namespace dd::literals;
        IBMQ_Yorktown.loadCouplingMap(AvailableArchitecture::IBMQ_Yorktown);
        IBMQ_London.loadCouplingMap(test_architecture_dir + "ibmq_london.arch");
        IBMQ_London.loadProperties(test_calibration_dir + "ibmq_london.csv");
        IBM_QX4.loadCouplingMap(AvailableArchitecture::IBM_QX4);

        yorktown_optimizer = std::make_unique<CliffordOptimizer>();
        yorktown_optimizer->setArchitecture(IBMQ_Yorktown);

        london_optimizer = std::make_unique<CliffordOptimizer>();
        london_optimizer->setArchitecture(IBMQ_London);

        qx4_optimizer = std::make_unique<CliffordOptimizer>();
        qx4_optimizer->setArchitecture(IBM_QX4);
    }
};

INSTANTIATE_TEST_SUITE_P(
        CliffordOptimizer, TestCliffordSynthesis,
        testing::Values(
                "destabilizer.txt"));

TEST_P(TestCliffordSynthesis, SimpleOptimization) {
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            qx4_optimizer->nqubits           = 2;
            qx4_optimizer->initial_timesteps = 10;
            Tableau::initTableau(qx4_optimizer->initialTableau, 2);
            qx4_optimizer->targetTableau = tableau;
            qx4_optimizer->optimize();
            tableau.clear();


            qx4_optimizer->optimal_results.dump(std::cout);
            EXPECT_EQ(qx4_optimizer->optimal_results.result, OptResult::SAT);
        }
    }
}

TEST(TestCliffordSynthesis, SanityCheck) {
    qc::QuantumComputation qc{};
    CliffordOptimizer      optimizer{};
    qc.addQubitRegister(2U);
    qc.h(0);
    qc.h(0);
    qc.h(0);
    qc.h(0);
    qc.h(0);

    optimizer.nqubits           = 2;
    optimizer.initial_timesteps = 2;
    Tableau::initTableau(optimizer.initialTableau, optimizer.nqubits);
    Tableau::generateTableau(optimizer.targetTableau, qc);

    optimizer.optimize();

    EXPECT_EQ(optimizer.optimal_results.depth, 1);
}

TEST_P(TestCliffordSynthesis, TestDepthOpt) {
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            qx4_optimizer->nqubits           = 2;
            qx4_optimizer->initial_timesteps = 10;
            qx4_optimizer->target            = OptTarget::DEPTH;
            Tableau::initTableau(qx4_optimizer->initialTableau, 2);
            qx4_optimizer->targetTableau = tableau;
            qx4_optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4_optimizer->optimal_results.result, OptResult::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestFidelityOpt) {
    //util::init();
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            london_optimizer->nqubits           = 2;
            london_optimizer->initial_timesteps = 5;
            london_optimizer->target            = OptTarget::FIDELITY;
            Tableau::initTableau(london_optimizer->initialTableau, 2);
            london_optimizer->targetTableau = tableau;
            london_optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(london_optimizer->optimal_results.result, OptResult::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestCNOTONLYOpt) {
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            qx4_optimizer->nqubits           = 2;
            qx4_optimizer->initial_timesteps = 10;
            qx4_optimizer->target            = OptTarget::GATES_ONLY_CNOT;
            Tableau::initTableau(qx4_optimizer->initialTableau, 2);
            qx4_optimizer->targetTableau = tableau;
            qx4_optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4_optimizer->optimal_results.result, OptResult::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestStartLow) {
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            qx4_optimizer->nqubits           = 2;
            qx4_optimizer->initial_timesteps = 10;
            qx4_optimizer->target            = OptTarget::GATES;
            qx4_optimizer->strategy          = OptimizingStrategy::StartLow;
            Tableau::initTableau(qx4_optimizer->initialTableau, 2);
            qx4_optimizer->targetTableau = tableau;
            qx4_optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4_optimizer->optimal_results.result, OptResult::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestStartHigh) {
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            qx4_optimizer->nqubits           = 2;
            qx4_optimizer->initial_timesteps = 50;
            qx4_optimizer->target            = OptTarget::GATES;
            qx4_optimizer->strategy          = OptimizingStrategy::StartHigh;
            Tableau::initTableau(qx4_optimizer->initialTableau, 2);
            qx4_optimizer->targetTableau = tableau;
            qx4_optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4_optimizer->optimal_results.result, OptResult::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestMinMax) {
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            qx4_optimizer->nqubits           = 2;
            qx4_optimizer->initial_timesteps = 10;
            qx4_optimizer->target            = OptTarget::GATES;
            qx4_optimizer->strategy          = OptimizingStrategy::MinMax;
            Tableau::initTableau(qx4_optimizer->initialTableau, 2);
            qx4_optimizer->targetTableau = tableau;
            qx4_optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4_optimizer->optimal_results.result, OptResult::SAT);
        }
    }
}

TEST_P(TestCliffordSynthesis, TestSplitIter) {
    auto&   input_file = GetParam();
    Tableau tableau{};
    if (input_file.find(".txt") != std::string::npos) {
        auto is = std::ifstream(test_example_dir + input_file);
        if (!is.good()) {
            FATAL() << "Error opening file " << input_file;
        }

        std::string line;
        while (std::getline(is, line)) {
            tableau.importString(line);
            qx4_optimizer->nqubits           = 2;
            qx4_optimizer->initial_timesteps = 10;
            qx4_optimizer->target            = OptTarget::GATES;
            qx4_optimizer->strategy          = OptimizingStrategy::SplitIter;
            Tableau::initTableau(qx4_optimizer->initialTableau, 2);
            qx4_optimizer->targetTableau = tableau;
            qx4_optimizer->optimize();
            tableau.clear();

            EXPECT_EQ(qx4_optimizer->optimal_results.result, OptResult::UNDEF);
        }
    }
}
