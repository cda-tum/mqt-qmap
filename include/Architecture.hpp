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
#include <vector>

constexpr unsigned short GATES_OF_BIDIRECTIONAL_SWAP  = 3U;
constexpr unsigned short GATES_OF_UNIDIRECTIONAL_SWAP = 7U;
constexpr unsigned short GATES_OF_DIRECTION_REVERSE   = 4U;
constexpr unsigned short GATES_OF_TELEPORTATION       = 7U;

constexpr int COST_SINGLE_QUBIT_GATE = 1;
constexpr int COST_CNOT_GATE         = 10;
constexpr int COST_MEASUREMENT       = 10;
constexpr int COST_UNIDIRECTIONAL_SWAP =
    3 * COST_CNOT_GATE + 4 * COST_SINGLE_QUBIT_GATE;
constexpr int COST_BIDIRECTIONAL_SWAP = 3 * COST_CNOT_GATE;
constexpr int COST_TELEPORTATION =
    2 * COST_CNOT_GATE + COST_MEASUREMENT + 4 * COST_SINGLE_QUBIT_GATE;
constexpr int COST_DIRECTION_REVERSE = 4 * COST_SINGLE_QUBIT_GATE;

class Architecture {
  static constexpr bool VERBOSE = false;

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
      [[nodiscard]] bool available(const KeyType& key) const {
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

    [[nodiscard]] unsigned short getNqubits() const { return nq; }
    void                         setNqubits(unsigned short nqs) { nq = nqs; }

    Property<unsigned short, Property<qc::OpType, double>>
        singleQubitErrorRate{};
    Property<unsigned short,
             Property<unsigned short, Property<qc::OpType, double>>>
                                          twoQubitErrorRate{};
    Property<unsigned short, double>      readoutErrorRate{};
    Property<unsigned short, double>      t1Time{};
    Property<unsigned short, double>      t2Time{};
    Property<unsigned short, double>      qubitFrequency{};
    Property<unsigned short, std::string> calibrationDate{};

    // convenience functions
    void setSingleQubitErrorRate(unsigned short     qubit,
                                 const std::string& operation,
                                 double             errorRate) {
      singleQubitErrorRate.get(qubit).set(qc::opTypeFromString(operation),
                                          errorRate);
    }
    [[nodiscard]] double
    getSingleQubitErrorRate(unsigned short     qubit,
                            const std::string& operation) const {
      return singleQubitErrorRate.get(qubit).get(
          qc::opTypeFromString(operation));
    }
    [[nodiscard]] double
    getAverageSingleQubitErrorRate(const unsigned short qubit) const {
      double avgErrorRate = 0.0;
      for (const auto& [opType, error] :
           singleQubitErrorRate.get(qubit).get()) {
        avgErrorRate += error;
      }
      return avgErrorRate /
             static_cast<double>(singleQubitErrorRate.get(qubit).get().size());
    }

    void setTwoQubitErrorRate(unsigned short qubit1, unsigned short qubit2,
                              double             errorRate,
                              const std::string& operation = "cx") {
      twoQubitErrorRate.get(qubit1).get(qubit2).set(
          qc::opTypeFromString(operation), errorRate);
    }
    [[nodiscard]] double
    getTwoQubitErrorRate(unsigned short qubit1, unsigned short qubit2,
                         const std::string& operation = "cx") const {
      return twoQubitErrorRate.get(qubit1).get(qubit2).get(
          qc::opTypeFromString(operation));
    }
    [[nodiscard]] bool
    twoQubitErrorRateAvailable(unsigned short qubit1, unsigned short qubit2,
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
      for (unsigned short i = 0; i < nq; ++i) {
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
    std::string    name{};
    unsigned short nq{};
  };

  void loadCouplingMap(std::istream& is);
  void loadCouplingMap(std::istream&& is);
  void loadCouplingMap(const std::string& filename);
  void loadCouplingMap(unsigned short nQ, const CouplingMap& cm);
  void loadCouplingMap(AvailableArchitecture architecture);
  void loadProperties(std::istream& is);
  void loadProperties(std::istream&& is);
  void loadProperties(const std::string& filename);
  void loadProperties(const Properties& properties);

  Architecture() = default;
  explicit Architecture(const std::string& cm_filename) {
    loadCouplingMap(cm_filename);
  }
  Architecture(const std::string& cm_filename,
               const std::string& props_filename)
      : Architecture(cm_filename) {
    loadProperties(props_filename);
  }

  Architecture(unsigned short nQ, const CouplingMap& couplingMap);
  Architecture(unsigned short nQ, const CouplingMap& couplingMap,
               const Properties& properties);

  [[nodiscard]] unsigned short getNqubits() const { return nqubits; }
  void                         setNqubits(unsigned short nQ) { nqubits = nQ; }

  [[nodiscard]] const std::string& getName() const { return name; }
  void setName(const std::string& architectureName) { name = architectureName; }

  [[nodiscard]] const CouplingMap& getCouplingMap() const {
    return couplingMap;
  }
  [[nodiscard]] CouplingMap& getCouplingMap() { return couplingMap; }
  void                       setCouplingMap(const CouplingMap& cm) {
                          couplingMap = cm;
                          createDistanceTable();
  }

  CouplingMap& getCurrentTeleportations() { return current_teleportations; }
  std::vector<std::pair<short, short>>& getTeleportationQubits() {
    return teleportationQubits;
  }

  [[nodiscard]] const Matrix& getDistanceTable() const { return distanceTable; }

  [[nodiscard]] const Properties& getProperties() const { return properties; }
  [[nodiscard]] Properties&       getProperties() { return properties; }
  void                            setProperties(const Properties& props) {
                               properties = props;
                               createFidelityTable();
  }

  [[nodiscard]] const Matrix& getFidelityTable() const { return fidelityTable; }

  [[nodiscard]] const std::vector<double>& getSingleQubitFidelities() const {
    return singleQubitFidelities;
  }

  [[nodiscard]] bool bidirectional() const { return isBidirectional; }

  [[nodiscard]] bool isArchitectureAvailable() {
    return !(name.empty()) && nqubits != 0;
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

  [[nodiscard]] double distance(unsigned short control,
                                unsigned short target) const {
    if (current_teleportations.empty()) {
      return distanceTable.at(control).at(target);
    } else {
      return static_cast<double>(bfs(control, target, current_teleportations));
    }
  }

  [[nodiscard]] std::set<unsigned short> getQubitSet() const {
    std::set<unsigned short> result{};
    for (int i = 0; i < nqubits; ++i) {
      result.insert(result.end(),
                    i); // should be constant with gcc, or at most O(nqubits)
    }
    return result;
  }

  unsigned long minimumNumberOfSwaps(std::vector<unsigned short>& permutation,
                                     long                         limit = -1);
  void          minimumNumberOfSwaps(
               std::vector<unsigned short>&                            permutation,
               std::vector<std::pair<unsigned short, unsigned short>>& swaps);

  struct Node {
    unsigned long                                          nswaps = 0;
    std::vector<std::pair<unsigned short, unsigned short>> swaps{};
    std::unordered_map<unsigned short, unsigned short>     permutation{};

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
  getCouplingLimit(const std::set<unsigned short>& qubitChoice) const;

  void getHighestFidelityCouplingMap(unsigned short subsetSize,
                                     CouplingMap&   couplingMap);
  [[nodiscard]] std::vector<std::set<unsigned short>>
       getAllConnectedSubsets(unsigned short subsetSize);
  void getReducedCouplingMaps(unsigned short            subsetSize,
                              std::vector<CouplingMap>& couplingMaps);
  void getReducedCouplingMap(const std::set<unsigned short>& qubitChoice,
                             CouplingMap&                    couplingMap);
  [[nodiscard]] static double
  getAverageArchitectureFidelity(const CouplingMap&              couplingMap,
                                 const std::set<unsigned short>& qubitChoice,
                                 const Properties&               props);

  [[nodiscard]] static std::vector<unsigned short>
  getQubitList(const CouplingMap& couplingMap);

  static bool isConnected(const std::set<unsigned short>& qubitChoice,
                          const CouplingMap&              reducedCouplingMap);

protected:
  std::string                          name;
  unsigned short                       nqubits                = 0;
  CouplingMap                          couplingMap            = {};
  CouplingMap                          current_teleportations = {};
  bool                                 isBidirectional        = true;
  Matrix                               distanceTable          = {};
  std::vector<std::pair<short, short>> teleportationQubits{};

  Properties          properties            = {};
  Matrix              fidelityTable         = {};
  std::vector<double> singleQubitFidelities = {};

  void createDistanceTable();
  void createFidelityTable();

  static double cost_heuristic_bidirectional(const Dijkstra::Node& node) {
    auto length = node.cost - 1;
    if (node.contains_correct_edge) {
      return length * COST_BIDIRECTIONAL_SWAP;
    } else {
      throw QMAPException("In a bidrectional architecture it should not happen "
                          "that a node does not contain the right edge.");
    }
  }

  static double cost_heuristic_unidirectional(const Dijkstra::Node& node) {
    auto length = node.cost - 1;
    if (node.contains_correct_edge) {
      return length * COST_UNIDIRECTIONAL_SWAP;
    } else {
      return length * COST_UNIDIRECTIONAL_SWAP + COST_DIRECTION_REVERSE;
    }
  }

  // added for teleportation
  static bool contains(const std::vector<int>& v, const int e) {
    return std::find(v.begin(), v.end(), e) != v.end();
  }
  [[nodiscard]] unsigned long bfs(unsigned short start, unsigned short goal,
                                  const std::set<Edge>& teleportations) const;

  static std::size_t findCouplingLimit(const CouplingMap& cm, int nQubits);
  static std::size_t
  findCouplingLimit(const CouplingMap& cm, int nQubits,
                    const std::set<unsigned short>& qubitChoice);
  static void
  findCouplingLimit(unsigned short node, int curSum,
                    const std::vector<std::vector<unsigned short>>& connections,
                    std::vector<int>& d, std::vector<bool>& visited);
};
