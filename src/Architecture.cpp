//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include "Architecture.hpp"

#include "utils.hpp"

#include <utility>

void Architecture::loadCouplingMap(AvailableArchitecture architecture) {
  std::stringstream ss{getCouplingMapSpecification(architecture)};
  name = toString(architecture);
  loadCouplingMap(ss);
}

void Architecture::loadCouplingMap(std::istream& is) {
  loadCouplingMap(std::move(is));
}

void Architecture::loadCouplingMap(const std::string& filename) {
  const std::size_t slash = filename.find_last_of('/');
  const std::size_t dot   = filename.find_last_of('.');
  name                    = filename.substr(slash + 1, dot - slash - 1);
  auto ifs                = std::ifstream(filename);
  if (ifs.good()) {
    this->loadCouplingMap(ifs);
  } else {
    throw QMAPException("Error opening coupling map file.");
  }
}

void Architecture::loadCouplingMap(std::istream&& is) {
  couplingMap.clear();
  properties.clear();
  std::string line;

  const auto  rNqubits = std::regex("([0-9]+)");
  const auto  rEdge    = std::regex("([0-9]+) ([0-9]+)");
  std::smatch m;

  // get number of qubits
  if (std::getline(is, line)) {
    if (std::regex_search(line, m, rNqubits)) {
      nqubits = static_cast<std::uint16_t>(std::stoul(m.str(1)));
    } else {
      throw QMAPException("No qubit count found in coupling map file: " + line);
    }
  } else {
    throw QMAPException("Error reading coupling map file.");
  }
  // load edges
  while (std::getline(is, line)) {
    if (std::regex_search(line, m, rEdge)) {
      auto v1 = static_cast<std::uint16_t>(std::stoul(m.str(1)));
      auto v2 = static_cast<std::uint16_t>(std::stoul(m.str(2)));
      couplingMap.emplace(v1, v2);
    } else {
      throw QMAPException("Could not identify edge in coupling map file: " +
                          line);
    }
  }
  createDistanceTable();
  createFidelityDistanceTable();
}

void Architecture::loadCouplingMap(std::uint16_t nQ, const CouplingMap& cm) {
  nqubits     = nQ;
  couplingMap = cm;
  properties.clear();
  name = "generic_" + std::to_string(nQ);
  createDistanceTable();
  createFidelityDistanceTable();
}

void Architecture::loadProperties(std::istream& is) {
  loadProperties(std::move(is));
}

void Architecture::loadProperties(const std::string& filename) {
  const std::size_t slash = filename.find_last_of('/');
  const std::size_t dot   = filename.find_last_of('.');
  properties.setName(filename.substr(slash + 1, dot - slash - 1));
  if (!isArchitectureAvailable()) {
    name = properties.getName();
  }
  auto ifs = std::ifstream(filename);
  if (ifs.good()) {
    this->loadProperties(ifs);
  } else {
    throw QMAPException("Error opening properties file.");
  }
}

void Architecture::loadProperties(std::istream&& is) {
  static const auto SINGLE_QUBIT_GATES = {"id", "u1", "u2", "u3",
                                          "rz", "sx", "x"};

  properties.clear();

  double averageCNOTFidelity = 0.0;
  int    numCNOTFidelities   = 0;

  std::string      line;
  const std::regex regexDoubleFidelity = std::regex(
      R"(((\d+).?(\d+):\W*?(-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?)))");
  std::smatch sMatch;
  std::getline(is, line); // skip first line
  // load edges
  std::uint16_t qubitNumber = 0U;
  while (std::getline(is, line)) {
    std::vector<std::string> data{};
    parseLine(line, ',', {'\"'}, {'\\'}, data);
    properties.t1Time.set(qubitNumber, std::stod(data.at(1U)));
    properties.t2Time.set(qubitNumber, std::stod(data.at(2U)));
    properties.qubitFrequency.set(qubitNumber, std::stod(data.at(3U)));
    properties.readoutErrorRate.set(qubitNumber, std::stod(data.at(4U)));
    // csv file reports average single qubit fidelities
    for (const auto& operation : SINGLE_QUBIT_GATES) {
      properties.setSingleQubitErrorRate(qubitNumber, operation,
                                         std::stod(data.at(5)));
    }
    // only try to parse CNOT fidelities if there are any
    if (data.size() >= 7U) {
      std::string s = data[6U];
      while (std::regex_search(s, sMatch, regexDoubleFidelity)) {
        auto a = static_cast<std::uint16_t>(std::stoul(sMatch.str(2U)));
        auto b = static_cast<std::uint16_t>(std::stoul(sMatch.str(3U)));
        if (!isArchitectureAvailable()) {
          couplingMap.emplace(a, b);
        }
        // calc moving average
        averageCNOTFidelity = averageCNOTFidelity + (std::stod(sMatch.str(4U)) -
                                                     averageCNOTFidelity) /
                                                        ++numCNOTFidelities;
        properties.setTwoQubitErrorRate(a, b, std::stod(sMatch.str(4U)));
        s = sMatch.suffix().str();
      }
    }
    if (data.size() == 8U) {
      properties.calibrationDate.set(qubitNumber, data[7U]);
    }
    ++qubitNumber;
  }

  if (isArchitectureAvailable()) {
    for (const auto& edge : couplingMap) {
      if (!properties.twoQubitErrorRateAvailable(edge.first, edge.second)) {
        properties.setTwoQubitErrorRate(edge.first, edge.second,
                                        averageCNOTFidelity);
      }
    }
  }
  properties.setNqubits(qubitNumber);
  if (!isArchitectureAvailable()) {
    nqubits = qubitNumber;
    createDistanceTable();
  }

  createFidelityTable();
  createFidelityDistanceTable();
}

void Architecture::loadProperties(const Properties& props) {
  if (!isArchitectureAvailable()) {
    for (const auto& [control, targetProps] : props.twoQubitErrorRate.get()) {
      for (const auto& [target, errorRate] : targetProps.get()) {
        couplingMap.emplace(control, target);
      }
    }
    nqubits = props.getNqubits();
    name    = "generic_" + std::to_string(nqubits);
    createDistanceTable();
  }
  properties = props;
  createFidelityTable();
  createFidelityDistanceTable();
}

Architecture::Architecture(const std::uint16_t nQ, const CouplingMap& cm) {
  loadCouplingMap(nQ, cm);
}

Architecture::Architecture(const std::uint16_t nQ, const CouplingMap& cm,
                           const Properties& props)
    : Architecture(nQ, cm) {
  loadProperties(props);
}

void Architecture::createDistanceTable() {
  Matrix<double> edgeWeights(
      nqubits,
      std::vector<double>(nqubits, std::numeric_limits<double>::max()));
  for (const auto& edge : couplingMap) {
    if (couplingMap.find({edge.second, edge.first}) == couplingMap.end()) {
      edgeWeights.at(edge.second).at(edge.first) = COST_UNIDIRECTIONAL_SWAP;
      edgeWeights.at(edge.first).at(edge.second) = COST_UNIDIRECTIONAL_SWAP;
    } else {
      edgeWeights.at(edge.first).at(edge.second) = COST_BIDIRECTIONAL_SWAP;
    }
  }

  Dijkstra::buildTable(nqubits, couplingMap, distanceTable, edgeWeights,
                       Architecture::dijkstraNodeToCostNonFidelity);
}

void Architecture::createFidelityDistanceTable() {
  for (const auto& edge : couplingMap) {
    if (!properties.twoQubitErrorRateAvailable(edge.first, edge.second)) {
      fidelityDistanceTable = std::vector<std::vector<double>>(
          nqubits, std::vector<double>(nqubits, 0));
      return;
    }
  }

  Dijkstra::buildTable(nqubits, couplingMap, fidelityDistanceTable,
                       swapFidelityCost,
                       Architecture::dijkstraNodeToCostFidelity);
}

void Architecture::createFidelityTable() {
  fidelityTable.clear();
  fidelityTable.resize(nqubits, std::vector<double>(nqubits, 0.0));
  twoQubitFidelityCost.clear();
  twoQubitFidelityCost.resize(
      nqubits,
      std::vector<double>(nqubits, std::numeric_limits<double>::max()));
  swapFidelityCost.clear();
  swapFidelityCost.resize(
      nqubits,
      std::vector<double>(nqubits, std::numeric_limits<double>::max()));

  singleQubitFidelities.resize(nqubits, 1.0);
  singleQubitFidelityCost.resize(nqubits, 0.0);

  for (const auto& [qubit, operationProps] :
       properties.singleQubitErrorRate.get()) {
    singleQubitFidelities[qubit] =
        1.0 - properties.getAverageSingleQubitErrorRate(qubit);
    singleQubitFidelityCost[qubit] = -std::log2(singleQubitFidelities[qubit]);
  }

  for (const auto& [first, second] : couplingMap) {
    if (properties.twoQubitErrorRateAvailable(first, second)) {
      fidelityTable[first][second] =
          1.0 - properties.getTwoQubitErrorRate(first, second);
      twoQubitFidelityCost[first][second] =
          -std::log2(fidelityTable[first][second]);
      if (couplingMap.find({second, first}) == couplingMap.end()) {
        // CNOT reversal (unidirectional edge q1 -> q2):
        // CX(q2,q1) = H(q1) H(q2) CX(q1,q2) H(q1) H(q2)
        twoQubitFidelityCost[second][first] =
            twoQubitFidelityCost[first][second] +
            2 * singleQubitFidelityCost[first] +
            2 * singleQubitFidelityCost[second];
        // SWAP decomposition (unidirectional edge q1 -> q2):
        // SWAP(q1,q2) = CX(q1,q2) H(q1) H(q2) CX(q1,q2) H(q1) H(q2) CX(q1,q2)
        swapFidelityCost[first][second] =
            3 * twoQubitFidelityCost[first][second] +
            2 * singleQubitFidelityCost[first] +
            2 * singleQubitFidelityCost[second];
        swapFidelityCost[second][first] = swapFidelityCost[first][second];
      } else {
        // SWAP decomposition (bidirectional edge q1 <-> q2):
        // SWAP(q1,q2) = CX(q1,q2) CX(q2,q1) CX(q1,q2)
        swapFidelityCost[first][second] =
            3 * twoQubitFidelityCost[first][second];
      }
    }
  }
}

std::uint64_t
Architecture::minimumNumberOfSwaps(std::vector<std::uint16_t>& permutation,
                                   std::int64_t                limit) {
  const bool tryToAbortEarly = (limit != -1);

  // consolidate used qubits
  QubitSubset qubits{};
  for (const auto& q : permutation) {
    qubits.insert(q);
  }

  // create map for goal permutation
  std::unordered_map<std::uint16_t, std::uint16_t> goalPermutation{};
  std::uint16_t                                    count    = 0U;
  bool                                             identity = true;
  for (const auto q : qubits) {
    goalPermutation.emplace(q, permutation.at(count));
    if (q != permutation.at(count)) {
      identity = false;
    }
    ++count;
  }

  if (identity) {
    return 0U;
  }

  // create selection of swap possibilities
  std::set<Edge> possibleSwaps{};
  for (const auto& edge : couplingMap) {
    // only use SWAPs between qubits that are currently being considered
    if (qubits.count(edge.first) == 0 || qubits.count(edge.second) == 0) {
      continue;
    }

    if (!bidirectional() ||
        (possibleSwaps.count(edge) == 0 &&
         possibleSwaps.count({edge.second, edge.first}) == 0)) {
      possibleSwaps.emplace(edge);
    }
  }

  Node start{};
  // start with identity permutation
  for (std::uint16_t i = 0U; i < nqubits; ++i) {
    start.permutation.emplace(i, i);
  }

  auto priority = [](const Node& x, const Node& y) {
    return x.nswaps > y.nswaps;
  };
  std::priority_queue<Node, std::vector<Node>, decltype(priority)> queue(
      priority);
  queue.push(start);

  while (!queue.empty()) {
    const Node current = queue.top();
    queue.pop();

    // in case no solution has been found using less than `limit` swaps, search
    // can be aborted
    if (tryToAbortEarly &&
        current.nswaps >= static_cast<std::uint64_t>(limit)) {
      return static_cast<std::uint64_t>(limit + 1U);
    }

    for (const auto& swap : possibleSwaps) {
      Node next = current;

      // apply and insert swap
      std::swap(next.permutation.at(swap.first),
                next.permutation.at(swap.second));
      next.nswaps++;
      bool done = true;
      for (const auto& assignment : goalPermutation) {
        if (next.permutation.at(assignment.first) != assignment.second) {
          done = false;
          break;
        }
      }

      if (done) {
        return next.nswaps;
      }
      queue.push(next);
    }
  }

  return start.nswaps;
}

void Architecture::minimumNumberOfSwaps(std::vector<std::uint16_t>& permutation,
                                        std::vector<Edge>&          swaps) {
  // consolidate used qubits
  QubitSubset qubits{};
  for (const auto& q : permutation) {
    qubits.insert(q);
  }

  if (qubits.size() != permutation.size()) {
    throw std::runtime_error(
        "Architecture::minimumNumberOfSwaps: permutation contains duplicates");
  }

  // create map for goal permutation
  std::unordered_map<std::uint16_t, std::uint16_t> goalPermutation{};
  std::uint16_t                                    count    = 0;
  bool                                             identity = true;
  for (const auto q : qubits) {
    goalPermutation.emplace(q, permutation.at(count));
    if (q != permutation.at(count)) {
      identity = false;
    }
    ++count;
  }

  if (identity) {
    return;
  }

  // create selection of swap possibilities
  std::set<Edge> possibleSwaps{};
  for (const auto& edge : couplingMap) {
    // only use SWAPs between qubits that are currently being considered
    if (qubits.count(edge.first) == 0 || qubits.count(edge.second) == 0) {
      continue;
    }

    if (!bidirectional() ||
        (possibleSwaps.count(edge) == 0 &&
         possibleSwaps.count({edge.second, edge.first}) == 0)) {
      possibleSwaps.emplace(edge);
    }
  }

  swaps.clear();
  Node start{};

  // start with identity permutation
  for (std::uint16_t i = 0U; i < nqubits; ++i) {
    start.permutation.emplace(i, i);
  }

  auto priority = [](const Node& x, const Node& y) {
    return x.swaps.size() > y.swaps.size();
  };
  std::priority_queue<Node, std::vector<Node>, decltype(priority)> queue(
      priority);
  queue.push(start);

  while (!queue.empty()) {
    const Node current = queue.top();
    queue.pop();

    for (const auto& swap : possibleSwaps) {
      Node next = current;
      // continue if the same swap was applied earlier
      if (!next.swaps.empty() && next.swaps.back() == swap) {
        continue;
      }

      // apply and insert swap
      std::swap(next.permutation.at(swap.first),
                next.permutation.at(swap.second));
      next.swaps.emplace_back(swap);
      next.nswaps++;

      bool done = true;
      for (const auto& assignment : goalPermutation) {
        if (next.permutation.at(assignment.first) != assignment.second) {
          done = false;
          break;
        }
      }

      if (done) {
        swaps = next.swaps;
        return;
      }
      queue.push(next);
    }
  }
}

std::size_t Architecture::getCouplingLimit() const {
  return findCouplingLimit(getCouplingMap(), getNqubits());
}

std::size_t
Architecture::getCouplingLimit(const QubitSubset& qubitChoice) const {
  return findCouplingLimit(getCouplingMap(), getNqubits(), qubitChoice);
}

std::uint64_t Architecture::bfs(const std::uint16_t   start,
                                const std::uint16_t   goal,
                                const std::set<Edge>& teleportations) const {
  std::queue<std::vector<int>> queue;
  std::vector<int>             v;
  v.push_back(start);
  queue.push(v);
  std::vector<std::vector<int>> solutions;

  std::uint64_t length = 0U;
  std::set<int> successors;
  while (!queue.empty()) {
    v = queue.front();
    queue.pop();
    const int current = v[v.size() - 1];
    if (current == goal) {
      length = v.size();
      solutions.push_back(v);
      break;
    }
    successors.clear();
    for (const auto& edge : getCouplingMap()) {
      if (edge.first == current && !contains(v, edge.second)) {
        successors.insert(edge.second);
      }
      if (edge.second == current && !contains(v, edge.first)) {
        successors.insert(edge.first);
      }
    }
    for (const auto& edge : teleportations) {
      if (edge.first == current && !contains(v, edge.second)) {
        successors.insert(edge.second);
      }
      if (edge.second == current && !contains(v, edge.first)) {
        successors.insert(edge.first); // was v2 but this is probably wrong
      }
    }

    for (const auto successor : successors) {
      auto v2 = v;
      v2.push_back(successor);
      queue.push(v2);
    }
  }
  while (!queue.empty() && queue.front().size() == length) {
    if (queue.front()[queue.front().size() - 1] == goal) {
      solutions.push_back(queue.front());
    }
    queue.pop();
  }

  // TODO: different weight if this contains a teleportation
  for (const auto& s : solutions) {
    for (std::size_t j = 0; j < s.size() - 1; j++) {
      const Edge e{s[j], s[j + 1]};
      if (getCouplingMap().find(e) != getCouplingMap().end()) {
        return (length - 2) * 7;
      }
    }
  }

  if (length == 2 &&
      getCouplingMap().find(Edge{start, goal}) == getCouplingMap().end() &&
      getCouplingMap().find(Edge{goal, start}) == getCouplingMap().end()) {
    return 7;
  }

  return (length - 2) * 7 + 4;
}

std::size_t Architecture::findCouplingLimit(const CouplingMap&  cm,
                                            const std::uint16_t nQubits) {
  std::vector<std::vector<std::uint16_t>> connections;
  std::vector<std::uint16_t>              d;
  std::vector<bool>                       visited;
  connections.resize(nQubits);
  std::uint16_t maxSum = 0;
  for (const auto& edge : cm) {
    connections.at(edge.first).emplace_back(edge.second);
  }
  for (std::uint16_t q = 0; q < nQubits; ++q) {
    d.clear();
    d.resize(nQubits);
    std::fill(d.begin(), d.end(), 0);
    visited.clear();
    visited.resize(nQubits);
    std::fill(visited.begin(), visited.end(), false);
    findCouplingLimit(q, 0, connections, d, visited);
    auto it = std::max_element(d.begin(), d.end());
    if ((*it) > maxSum) {
      maxSum = (*it);
    }
  }
  return maxSum;
}

std::size_t Architecture::findCouplingLimit(const CouplingMap&  cm,
                                            const std::uint16_t nQubits,
                                            const QubitSubset&  qubitChoice) {
  std::vector<std::vector<std::uint16_t>> connections;
  std::vector<std::uint16_t>              d;
  std::vector<bool>                       visited;
  connections.resize(nQubits);
  std::uint16_t maxSum = 0;
  for (const auto& edge : cm) {
    if ((qubitChoice.count(edge.first) != 0U) &&
        (qubitChoice.count(edge.second) != 0U)) {
      connections.at(edge.first).emplace_back(edge.second);
    }
  }
  for (std::uint16_t q = 0; q < nQubits; ++q) {
    if (connections.at(q).empty()) {
      continue;
    }
    d.clear();
    d.resize(nQubits);
    std::fill(d.begin(), d.end(), 0);
    visited.clear();
    visited.resize(nQubits);
    std::fill(visited.begin(), visited.end(), false);
    findCouplingLimit(q, 0, connections, d, visited);
    auto it = std::max_element(d.begin(), d.end());
    if ((*it) > maxSum) {
      maxSum = (*it);
    }
  }
  return maxSum;
}

void Architecture::findCouplingLimit(
    const std::uint16_t node, const std::uint16_t curSum,
    const std::vector<std::vector<std::uint16_t>>& connections,
    std::vector<std::uint16_t>& d, std::vector<bool>& visited) {
  if (visited.at(node)) {
    return;
  }
  visited[node] = true;

  auto& elem = d[node];
  if (elem == 0 || elem > curSum) {
    elem = curSum;
  }
  if (connections.at(node).empty()) {
    visited[node] = false;
    return;
  }

  for (auto child : connections.at(node)) {
    findCouplingLimit(child, curSum + 1, connections, d, visited);
  }

  visited[node] = false;
}

void Architecture::getHighestFidelityCouplingMap(
    std::uint16_t subsetSize, CouplingMap& reducedMap) const {
  if (!isArchitectureAvailable()) {
    reducedMap = getFullyConnectedMap(subsetSize);
    return;
  }

  if (nqubits == subsetSize) {
    reducedMap = couplingMap;
    return;
  }

  double bestFidelity        = std::numeric_limits<double>::lowest();
  auto   allConnectedSubsets = getAllConnectedSubsets(subsetSize);

  for (const auto& qubitChoice : allConnectedSubsets) {
    CouplingMap map{};
    getReducedCouplingMap(qubitChoice, map);
    const auto currentFidelity =
        getAverageArchitectureFidelity(map, qubitChoice, properties);
    if (currentFidelity > bestFidelity) {
      reducedMap   = map;
      bestFidelity = currentFidelity;
    }
  }
}
std::vector<QubitSubset>
Architecture::getAllConnectedSubsets(std::uint16_t subsetSize) const {
  std::vector<QubitSubset> result{};
  if (!isArchitectureAvailable() || nqubits == subsetSize) {
    result.emplace_back(getQubitSet());
  } else if (nqubits < subsetSize) {
    throw QMAPException("Architecture too small!");
  } else {
    auto filter = [&](const QubitSubset& subset) {
      CouplingMap cm = {};
      Architecture::getReducedCouplingMap(subset, cm);
      return isConnected(subset, cm);
    };
    for (const auto& subset : subsets(getQubitSet(), subsetSize, filter)) {
      result.emplace_back(subset);
    }
  }
  return result;
}
void Architecture::getReducedCouplingMaps(
    std::uint16_t subsetSize, std::vector<CouplingMap>& couplingMaps) const {
  couplingMaps.clear();
  if (!isArchitectureAvailable()) {
    couplingMaps.emplace_back(getFullyConnectedMap(subsetSize));
  } else {
    for (const auto& qubitChoice : getAllConnectedSubsets(subsetSize)) {
      couplingMaps.emplace_back();
      getReducedCouplingMap(qubitChoice, couplingMaps.back());
    }
  }
}
void Architecture::getReducedCouplingMap(const QubitSubset& qubitChoice,
                                         CouplingMap&       reducedMap) const {
  reducedMap.clear();
  if (!isArchitectureAvailable()) {
    reducedMap =
        getFullyConnectedMap(static_cast<std::uint16_t>(qubitChoice.size()));
  } else {
    for (const auto& [q0, q1] : couplingMap) {
      if (qubitChoice.find(q0) != qubitChoice.end() &&
          qubitChoice.find(q1) != qubitChoice.end()) {
        reducedMap.emplace(q0, q1);
      }
    }
  }
}

double
Architecture::getAverageArchitectureFidelity(const CouplingMap& cm,
                                             const QubitSubset& qubitChoice,
                                             const Properties&  props) {
  if (props.empty()) {
    return 0.0;
  }
  double result = 1.0;
  for (const auto& [control, target] : cm) {
    if (props.twoQubitErrorRateAvailable(control, target)) {
      result *= (1.0 - props.getTwoQubitErrorRate(control, target));
    }
  }

  for (const auto& qubit : qubitChoice) {
    if (props.singleQubitErrorRate.available(qubit)) {
      result *= (1.0 - props.getAverageSingleQubitErrorRate(qubit));
    }
  }
  return result;
}

QubitSubset Architecture::getQubitSet(const CouplingMap& cm) {
  QubitSubset result{};
  for (const auto& [q0, q1] : cm) {
    result.emplace(q0);
    result.emplace(q1);
  }
  return result;
}

bool Architecture::isConnected(const QubitSubset& qubitChoice,
                               const CouplingMap& reducedCouplingMap) {
  QubitSubset reachedQubits{};
  reachedQubits.emplace(*(qubitChoice.begin()));
  dfs(*(qubitChoice.begin()), reachedQubits, reducedCouplingMap);
  return (reachedQubits == qubitChoice);
}

void Architecture::printCouplingMap(const CouplingMap& cm, std::ostream& os) {
  os << "{ ";
  for (const auto& edge : cm) {
    os << "(" << edge.first << " " << edge.second << ") ";
  }
  os << "}" << std::endl;
}
