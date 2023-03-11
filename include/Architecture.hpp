//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#pragma once

#include "configuration/AvailableArchitecture.hpp"
#include "nlohmann/json.hpp"
#include "utils.hpp"

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <regex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

constexpr std::uint8_t GATES_OF_BIDIRECTIONAL_SWAP  = 3U;
constexpr std::uint8_t GATES_OF_UNIDIRECTIONAL_SWAP = 7U;
constexpr std::uint8_t GATES_OF_DIRECTION_REVERSE   = 4U;
constexpr std::uint8_t GATES_OF_TELEPORTATION       = 7U;

constexpr std::uint32_t COST_SINGLE_QUBIT_GATE = 1;
constexpr std::uint32_t COST_CNOT_GATE         = 10;
constexpr std::uint32_t COST_MEASUREMENT       = 10;
constexpr std::uint32_t COST_UNIDIRECTIONAL_SWAP =
    3 * COST_CNOT_GATE + 4 * COST_SINGLE_QUBIT_GATE;
constexpr std::uint32_t COST_BIDIRECTIONAL_SWAP = 3 * COST_CNOT_GATE;
constexpr std::uint32_t COST_TELEPORTATION =
    2 * COST_CNOT_GATE + COST_MEASUREMENT + 4 * COST_SINGLE_QUBIT_GATE;
constexpr std::uint32_t COST_DIRECTION_REVERSE = 4 * COST_SINGLE_QUBIT_GATE;

enum class DirectionReversalStrategy { Identity, Hadamard, NotApplicable };

class Architecture {
public:
  class Properties {
  protected:
    template <class KeyType, class ValueType> class Property {
    public:
      Property() = default;

      [[nodiscard]] auto&       get(const KeyType& key) { return props[key]; }
      [[nodiscard]] const auto& get(const KeyType& key) const {
        return props.at(key);
      }
      [[nodiscard]] const auto& get() const { return props; }
      void set(const KeyType& key, const ValueType& value) {
        props[key] = value;
      }
      [[nodiscard, gnu::pure]] bool available(const KeyType& key) const {
        return props.find(key) != props.end();
      }
      void               clear() { props.clear(); }
      [[nodiscard]] bool empty() const { return props.empty(); }

    protected:
      std::map<KeyType, ValueType> props{};
    };

  public:
    Properties() = default;

    [[nodiscard]] std::string getName() const { return name; }
    void setName(const std::string& propertiesName) { name = propertiesName; }

    [[nodiscard]] std::uint16_t getNqubits() const { return nq; }
    void                        setNqubits(std::uint16_t nqs) { nq = nqs; }

    Property<std::uint16_t, Property<qc::OpType, double>>
        singleQubitErrorRate{};
    Property<std::uint16_t,
             Property<std::uint16_t, Property<qc::OpType, double>>>
                                         twoQubitErrorRate{};
    Property<std::uint16_t, double>      readoutErrorRate{};
    Property<std::uint16_t, double>      t1Time{};
    Property<std::uint16_t, double>      t2Time{};
    Property<std::uint16_t, double>      qubitFrequency{};
    Property<std::uint16_t, std::string> calibrationDate{};

    // convenience functions
    void setSingleQubitErrorRate(std::uint16_t      qubit,
                                 const std::string& operation,
                                 double             errorRate) {
      singleQubitErrorRate.get(qubit).set(qc::opTypeFromString(operation),
                                          errorRate);
    }
    [[nodiscard]] double
    getSingleQubitErrorRate(std::uint16_t      qubit,
                            const std::string& operation) const {
      return singleQubitErrorRate.get(qubit).get(
          qc::opTypeFromString(operation));
    }
    [[nodiscard]] double
    getAverageSingleQubitErrorRate(const std::uint16_t qubit) const {
      double avgErrorRate = 0.0;
      for (const auto& [opType, error] :
           singleQubitErrorRate.get(qubit).get()) {
        avgErrorRate += error;
      }
      return avgErrorRate /
             static_cast<double>(singleQubitErrorRate.get(qubit).get().size());
    }

    void setTwoQubitErrorRate(std::uint16_t qubit1, std::uint16_t qubit2,
                              double             errorRate,
                              const std::string& operation = "cx") {
      twoQubitErrorRate.get(qubit1).get(qubit2).set(
          qc::opTypeFromString(operation), errorRate);
    }
    [[nodiscard]] double
    getTwoQubitErrorRate(std::uint16_t qubit1, std::uint16_t qubit2,
                         const std::string& operation = "cx") const {
      return twoQubitErrorRate.get(qubit1).get(qubit2).get(
          qc::opTypeFromString(operation));
    }
    [[nodiscard]] bool
    twoQubitErrorRateAvailable(std::uint16_t qubit1, std::uint16_t qubit2,
                               const std::string& operation = "cx") const {
      return twoQubitErrorRate.available(qubit1) &&
             twoQubitErrorRate.get(qubit1).available(qubit2) &&
             twoQubitErrorRate.get(qubit1).get(qubit2).available(
                 qc::opTypeFromString(operation));
    }

    void clear() {
      singleQubitErrorRate.clear();
      twoQubitErrorRate.clear();
      readoutErrorRate.clear();
      t1Time.clear();
      t2Time.clear();
      qubitFrequency.clear();
      calibrationDate.clear();
    }

    [[nodiscard]] bool empty() const {
      return singleQubitErrorRate.empty() && twoQubitErrorRate.empty() &&
             readoutErrorRate.empty() && t1Time.empty() && t2Time.empty() &&
             qubitFrequency.empty() && calibrationDate.empty();
    }

    [[nodiscard]] nlohmann::json json() const {
      nlohmann::json json;
      if (empty()) {
        return json;
      }

      json["name"]   = name;
      json["qubits"] = {};
      for (std::uint16_t i = 0U; i < nq; ++i) {
        auto& qubitProperties = json["qubits"][std::to_string(i)];

        if (singleQubitErrorRate.available(i)) {
          auto& singleQubitErrorRates =
              qubitProperties["single_qubit_error_rate"];
          for (const auto& [operation, error] :
               singleQubitErrorRate.get(i).get()) {
            singleQubitErrorRates[qc::toString(operation)] = error;
          }
        }

        if (t1Time.available(i)) {
          qubitProperties["t1_time"] = t1Time.get(i);
        }
        if (t2Time.available(i)) {
          qubitProperties["t2_time"] = t2Time.get(i);
        }
        if (qubitFrequency.available(i)) {
          qubitProperties["frequency"] = qubitFrequency.get(i);
        }
        if (calibrationDate.available(i)) {
          qubitProperties["calibration_date"] = calibrationDate.get(i);
        }
        if (readoutErrorRate.available(i)) {
          qubitProperties["readout_error_rate"] = readoutErrorRate.get(i);
        }

        if (twoQubitErrorRate.available(i)) {
          auto& twoQubitErrorRates = qubitProperties["two_qubit_error_rate"];
          for (const auto& [qubit2, errorRates] :
               twoQubitErrorRate.get(i).get()) {
            const std::string pair =
                '(' + std::to_string(i) + ',' + std::to_string(qubit2) + ')';
            auto& qubits = twoQubitErrorRates[pair];
            for (const auto& [operation, error] : errorRates.get()) {
              qubits[qc::toString(operation)] = error;
            }
          }
        }
      }

      return json;
    }
    [[nodiscard]] std::string toString() const { return json().dump(2); }

  protected:
    std::string   name{};
    std::uint16_t nq{};
  };

  void loadCouplingMap(std::istream& is);
  void loadCouplingMap(std::istream&& is);
  void loadCouplingMap(const std::string& filename);
  void loadCouplingMap(std::uint16_t nQ, const CouplingMap& cm);
  void loadCouplingMap(AvailableArchitecture architecture);
  void loadProperties(std::istream& is);
  void loadProperties(std::istream&& is);
  void loadProperties(const std::string& filename);
  void loadProperties(const Properties& props);

  Architecture() = default;
  explicit Architecture(const std::string& cmFilename) {
    loadCouplingMap(cmFilename);
  }
  Architecture(const std::string& cmFilename, const std::string& propsFilename)
      : Architecture(cmFilename) {
    loadProperties(propsFilename);
  }

  Architecture(std::uint16_t nQ, const CouplingMap& cm);
  Architecture(std::uint16_t nQ, const CouplingMap& cm,
               const Properties& props);

  [[nodiscard]] std::uint16_t getNqubits() const { return nqubits; }
  void                        setNqubits(std::uint16_t nQ) { nqubits = nQ; }

  [[nodiscard]] const std::string& getName() const { return name; }
  void setName(const std::string& architectureName) { name = architectureName; }

  [[nodiscard]] const CouplingMap& getCouplingMap() const {
    return couplingMap;
  }
  [[nodiscard]] CouplingMap& getCouplingMap() { return couplingMap; }

  void setCouplingMap(const CouplingMap& cm) {
    couplingMap = cm;
    createDistanceTable();
  }

  CouplingMap& getCurrentTeleportations() { return currentTeleportations; }
  std::vector<std::pair<std::int16_t, std::int16_t>>& getTeleportationQubits() {
    return teleportationQubits;
  }

  [[nodiscard]] const Matrix& getDistanceTable() const { return distanceTable; }

  [[nodiscard]] const Properties& getProperties() const { return properties; }

  [[nodiscard]] Properties& getProperties() { return properties; }

  void setProperties(const Properties& props) {
    properties = props;
    createFidelityTable();
  }

  [[nodiscard]] const Matrix& getFidelityTable() const { return fidelityTable; }

  [[nodiscard]] const std::vector<double>& getSingleQubitFidelities() const {
    return singleQubitFidelities;
  }

  [[nodiscard]] bool bidirectional() const { return isBidirectional; }

  [[nodiscard]] bool isArchitectureAvailable() const {
    return !(name.empty()) && nqubits != 0;
  }
  [[nodiscard]] bool isCalibrationDataAvailable() const {
    return !(name.empty()) && !properties.empty();
  }

  void reset() {
    name    = "";
    nqubits = 0;
    couplingMap.clear();
    distanceTable.clear();
    isBidirectional = true;
    properties.clear();
    fidelityTable.clear();
    singleQubitFidelities.clear();
  }

  [[nodiscard]] double distance(std::uint16_t control,
                                std::uint16_t target) const {
    if (currentTeleportations.empty()) {
      return distanceTable.at(control).at(target);
    }
    return static_cast<double>(bfs(control, target, currentTeleportations));
  }

  [[nodiscard]] std::set<std::uint16_t> getQubitSet() const {
    std::set<std::uint16_t> result{};
    for (std::uint16_t i = 0; i < nqubits; ++i) {
      result.insert(result.end(),
                    i); // should be constant with gcc, or at most O(nqubits)
    }
    return result;
  }

  std::uint64_t minimumNumberOfSwaps(std::vector<std::uint16_t>& permutation,
                                     std::int64_t                limit = -1);
  void          minimumNumberOfSwaps(std::vector<std::uint16_t>& permutation,
                                     std::vector<Edge>&          swaps);

  struct Node {
    std::uint64_t                                    nswaps = 0U;
    std::vector<Edge>                                swaps{};
    std::unordered_map<std::uint16_t, std::uint16_t> permutation{};

    void print(std::ostream& out) {
      out << swaps.size() << ": ";
      for (const auto& p : permutation) {
        out << p.first << "->" << p.second << " ";
      }
      out << " | ";
      for (const auto& swap : swaps) {
        out << swap.first << "<->" << swap.second << " ";
      }
      out << std::endl;
    }
  };

  [[nodiscard]] std::size_t getCouplingLimit() const;
  [[nodiscard]] std::size_t
  getCouplingLimit(const std::set<std::uint16_t>& qubitChoice) const;

  void getHighestFidelityCouplingMap(std::uint16_t subsetSize,
                                     CouplingMap&  reducedMap) const;
  [[nodiscard]] std::vector<QubitSubset>
       getAllConnectedSubsets(std::uint16_t subsetSize) const;
  void getReducedCouplingMaps(std::uint16_t             subsetSize,
                              std::vector<CouplingMap>& couplingMaps) const;
  void getReducedCouplingMap(const QubitSubset& qubitChoice,
                             CouplingMap&       reducedMap) const;
  [[nodiscard]] static double
  getAverageArchitectureFidelity(const CouplingMap& cm,
                                 const QubitSubset& qubitChoice,
                                 const Properties&  props);

  [[nodiscard]] static QubitSubset getQubitSet(const CouplingMap& cm);
  [[nodiscard]] static std::vector<std::uint16_t>
  getQubitList(const CouplingMap& cm) {
    const auto qubitSet = getQubitSet(cm);
    return {qubitSet.begin(), qubitSet.end()};
  }

  static bool isConnected(const QubitSubset& qubitChoice,
                          const CouplingMap& reducedCouplingMap);

  static void printCouplingMap(const CouplingMap& cm, std::ostream& os);

  static DirectionReversalStrategy
                       getDirectionReversalStrategy(qc::OpType opType);
  static std::uint32_t computeCostDirectionReverse(qc::OpType opType);
  static std::uint32_t computeGatesDirectionReverse(qc::OpType opType);
  static bool          supportsDirectionReversal(qc::OpType opType);

protected:
  std::string                                        name;
  std::uint16_t                                      nqubits               = 0;
  CouplingMap                                        couplingMap           = {};
  CouplingMap                                        currentTeleportations = {};
  bool                                               isBidirectional = true;
  Matrix                                             distanceTable   = {};
  std::vector<std::pair<std::int16_t, std::int16_t>> teleportationQubits{};
  Properties                                         properties            = {};
  Matrix                                             fidelityTable         = {};
  std::vector<double>                                singleQubitFidelities = {};

  void createDistanceTable();
  void createFidelityTable();

  static double costHeuristicBidirectional(const Dijkstra::Node& node) {
    auto length = node.cost - 1;
    if (node.containsCorrectEdge) {
      return length * COST_BIDIRECTIONAL_SWAP;
    }
    throw QMAPException("In a bidrectional architecture it should not happen "
                        "that a node does not contain the right edge.");
  }

  static double costHeuristicUnidirectional(const Dijkstra::Node& node) {
    auto length = node.cost - 1;
    if (node.containsCorrectEdge) {
      return length * COST_UNIDIRECTIONAL_SWAP;
    }
    // TODO: distance has no gate context, but gates have different costs
    return length * COST_UNIDIRECTIONAL_SWAP +
           Architecture::computeCostDirectionReverse(qc::OpType::X);
  }

  // added for teleportation
  static bool contains(const std::vector<int>& v, const int e) {
    return std::find(v.begin(), v.end(), e) != v.end();
  }
  [[nodiscard]] std::uint64_t bfs(std::uint16_t start, std::uint16_t goal,
                                  const std::set<Edge>& teleportations) const;

  static std::size_t findCouplingLimit(const CouplingMap& cm,
                                       std::uint16_t      nQubits);
  static std::size_t
              findCouplingLimit(const CouplingMap& cm, std::uint16_t nQubits,
                                const std::set<std::uint16_t>& qubitChoice);
  static void findCouplingLimit(
      std::uint16_t node, std::uint16_t curSum,
      const std::vector<std::unordered_set<std::uint16_t>>& connections,
      std::vector<std::uint16_t>& d, std::vector<bool>& visited);
};
