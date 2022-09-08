/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://www.cda.cit.tum.de/research/ibm_qx_mapping/ for more information.
 */

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
    size_t slash = filename.find_last_of('/');
    size_t dot   = filename.find_last_of('.');
    name         = filename.substr(slash + 1, dot - slash - 1);
    auto ifs     = std::ifstream(filename);
    if (ifs.good())
        this->loadCouplingMap(ifs);
    else
        throw QMAPException("Error opening coupling map file.");
}

void Architecture::loadCouplingMap(std::istream&& is) {
    couplingMap.clear();
    properties.clear();
    std::string line;

    std::regex  r_nqubits = std::regex("([0-9]+)");
    std::regex  r_edge    = std::regex("([0-9]+) ([0-9]+)");
    std::smatch m;

    // get number of qubits
    if (std::getline(is, line)) {
        if (std::regex_match(line, m, r_nqubits)) {
            nqubits = static_cast<unsigned short>(std::stoul(m.str(1)));
        } else {
            throw QMAPException("No qubit count found in coupling map file: " + line);
        }
    } else {
        throw QMAPException("Error reading coupling map file.");
    }
    // load edges
    while (std::getline(is, line)) {
        if (std::regex_match(line, m, r_edge)) {
            auto v1 = static_cast<unsigned short>(std::stoul(m.str(1)));
            auto v2 = static_cast<unsigned short>(std::stoul(m.str(2)));
            couplingMap.emplace(v1, v2);
        } else {
            throw QMAPException("Could not identify edge in coupling map file: " + line);
        }
    }

    if (VERBOSE) {
        std::cout << "Coupling Map (" << nqubits << " qubits): ";
        for (const auto& edge: couplingMap) {
            std::cout << "(" << edge.first << "-" << edge.second << ") ";
        }
        std::cout << std::endl;
    }

    createDistanceTable();
}

void Architecture::loadCouplingMap(unsigned short nQ, const CouplingMap& cm) {
    nqubits     = nQ;
    couplingMap = cm;
    properties.clear();
    name = "generic_" + std::to_string(nQ);
    createDistanceTable();
}

void Architecture::loadProperties(std::istream& is) {
    loadProperties(std::move(is));
}

void Architecture::loadProperties(const std::string& filename) {
    size_t slash = filename.find_last_of('/');
    size_t dot   = filename.find_last_of('.');
    properties.setName(filename.substr(slash + 1, dot - slash - 1));
    if (!isArchitectureAvailable())
        name = properties.getName();
    auto ifs = std::ifstream(filename);
    if (ifs.good())
        this->loadProperties(ifs);
    else
        throw QMAPException("Error opening properties file.");
}

void Architecture::loadProperties(std::istream&& is) {
    static const auto SingleQubitGates = {"id", "u1", "u2", "u3", "rz", "sx", "x"};

    properties.clear();

    double averageCNOTFidelity = 0.0;
    int    numCNOTFidelities   = 0;

    std::string line;
    std::string word;
    std::regex  regexDoubleFidelity =
            std::regex(R"(((\d+).?(\d+):\W*?(-?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+-]?\d+)?)))");
    std::smatch sMatch;
    std::getline(is, line); //skip first line
    // load edges
    int qubitNumber = 0;
    while (std::getline(is, line)) {
        std::stringstream        ss(line);
        std::vector<std::string> data{};
        parse_line(line, ',', {'\"'}, {'\\'}, data);
        properties.t1Time.set(qubitNumber, std::stod(data.at(1)));
        properties.t2Time.set(qubitNumber, std::stod(data.at(2)));
        properties.qubitFrequency.set(qubitNumber, std::stod(data.at(3)));
        properties.readoutErrorRate.set(qubitNumber, std::stod(data.at(4)));
        // csv file reports average single qubit fidelities
        for (const auto& operation: SingleQubitGates) {
            properties.setSingleQubitErrorRate(qubitNumber, operation, std::stod(data.at(5)));
        }
        // only try to parse CNOT fidelities if there are any
        if (data.size() >= 7U) {
            std::string s = data[6];
            while (std::regex_search(s, sMatch, regexDoubleFidelity)) {
                auto a = static_cast<unsigned short>(std::stoul(sMatch.str(2U)));
                auto b = static_cast<unsigned short>(std::stoul(sMatch.str(3U)));
                if (!isArchitectureAvailable()) {
                    couplingMap.emplace(a, b);
                }
                // calc moving average
                averageCNOTFidelity = averageCNOTFidelity + (std::stod(sMatch.str(4U)) - averageCNOTFidelity) / ++numCNOTFidelities;
                properties.setTwoQubitErrorRate(a, b, std::stod(sMatch.str(4U)));
                s = sMatch.suffix().str();
            }
        }
        if (data.size() == 8U) {
            properties.calibrationDate.set(qubitNumber, data[7]);
        }
        ++qubitNumber;
    }

    if (isArchitectureAvailable()) {
        for (const auto& edge: couplingMap) {
            if (!properties.twoQubitErrorRateAvailable(edge.first, edge.second)) {
                properties.setTwoQubitErrorRate(edge.first, edge.second, averageCNOTFidelity);
            }
        }
    }
    properties.setNqubits(qubitNumber);
    if (!isArchitectureAvailable()) {
        nqubits = static_cast<unsigned short>(qubitNumber);
        createDistanceTable();
    }

    createFidelityTable();
}

void Architecture::loadProperties(const Properties& props) {
    if (!isArchitectureAvailable()) {
        for (const auto& [control, targetProps]: props.twoQubitErrorRate.get()) {
            for (const auto& [target, errorRate]: targetProps.get()) {
                couplingMap.emplace(control, target);
            }
        }
        nqubits = props.getNqubits();
        name    = "generic_" + std::to_string(nqubits);
        createDistanceTable();
    }
    properties = props;
    createFidelityTable();
}

Architecture::Architecture(unsigned short nQ, const CouplingMap& couplingMap) {
    loadCouplingMap(nQ, couplingMap);
}

Architecture::Architecture(unsigned short nQ, const CouplingMap& couplingMap, const Properties& props):
    Architecture(nQ, couplingMap) {
    loadProperties(props);
}

void Architecture::createDistanceTable() {
    for (const auto& edge: couplingMap) {
        if (couplingMap.find({edge.second, edge.first}) == couplingMap.end()) {
            isBidirectional = false;
            break;
        }
    }

    if (VERBOSE) {
        std::cout << "Architecture is bidirectional: " << (bidirectional() ? "yes" : "no") << std::endl;
    }

    if (isBidirectional) {
        Dijkstra::build_table(nqubits, couplingMap, distanceTable, Architecture::cost_heuristic_bidirectional);
    } else {
        Dijkstra::build_table(nqubits, couplingMap, distanceTable, Architecture::cost_heuristic_unidirectional);
    }
}

void Architecture::createFidelityTable() {
    fidelityTable.clear();
    fidelityTable.resize(nqubits, std::vector<double>(nqubits, 0.0));

    singleQubitFidelities.resize(nqubits, 0.0);

    for (const auto& [first, second]: couplingMap) {
        if (properties.twoQubitErrorRateAvailable(first, second)) {
            fidelityTable[first][second] = 1.0 - properties.getTwoQubitErrorRate(first, second);
        }
    }

    for (const auto& [qubit, operationProps]: properties.singleQubitErrorRate.get())
        singleQubitFidelities[qubit] = 1.0 - properties.getAverageSingleQubitErrorRate(qubit);
}

unsigned long Architecture::minimumNumberOfSwaps(std::vector<unsigned short>& permutation, long limit) {
    bool tryToAbortEarly = (limit != -1);

    // consolidate used qubits
    std::set<unsigned short> qubits{};
    for (const auto& q: permutation) {
        qubits.insert(q);
    }

    // create map for goal permutation
    std::unordered_map<unsigned short, unsigned short> goalPermutation{};
    unsigned short                                     count    = 0;
    bool                                               identity = true;
    for (const auto q: qubits) {
        goalPermutation.emplace(q, permutation.at(count));
        if (q != permutation.at(count)) {
            identity = false;
        }
        ++count;
    }

    if (identity) {
        return 0;
    }

    // create selection of swap possibilities
    std::set<std::pair<unsigned short, unsigned short>> possibleSwaps{};
    for (const auto& edge: couplingMap) {
        // only use SWAPs between qubits that are currently being considered
        if (qubits.count(edge.first) == 0 || qubits.count(edge.second) == 0) {
            continue;
        }

        if (!bidirectional() || (possibleSwaps.count(edge) == 0 && possibleSwaps.count({edge.second, edge.first}) == 0)) {
            possibleSwaps.emplace(edge);
        }
    }

    Node start{};
    // start with identity permutation
    for (unsigned short i = 0; i < nqubits; ++i) {
        start.permutation.emplace(i, i);
    }

    auto                                                             priority = [](const Node& x, const Node& y) { return x.nswaps > y.nswaps; };
    std::priority_queue<Node, std::vector<Node>, decltype(priority)> queue(priority);
    queue.push(start);

    while (!queue.empty()) {
        Node current = queue.top();
        queue.pop();

        // in case no solution has been found using less than `limit` swaps, search can be aborted
        if (tryToAbortEarly && current.nswaps >= static_cast<unsigned long>(limit)) {
            return limit + 1U;
        }

        for (const auto& swap: possibleSwaps) {
            Node next = current;

            // apply and insert swap
            std::swap(next.permutation.at(swap.first), next.permutation.at(swap.second));
            next.nswaps++;
            bool done = true;
            for (const auto& assignment: goalPermutation) {
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

void Architecture::minimumNumberOfSwaps(std::vector<unsigned short>& permutation, std::vector<std::pair<unsigned short, unsigned short>>& swaps) {
    // consolidate used qubits
    std::set<unsigned short> qubits{};
    for (const auto& q: permutation) {
        qubits.insert(q);
    }

    if (qubits.size() != permutation.size()) {
        throw std::runtime_error("Architecture::minimumNumberOfSwaps: permutation contains duplicates");
    }

    // create map for goal permutation
    std::unordered_map<unsigned short, unsigned short> goalPermutation{};
    unsigned short                                     count    = 0;
    bool                                               identity = true;
    for (const auto q: qubits) {
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
    std::set<std::pair<unsigned short, unsigned short>> possibleSwaps{};
    for (const auto& edge: couplingMap) {
        // only use SWAPs between qubits that are currently being considered
        if (qubits.count(edge.first) == 0 || qubits.count(edge.second) == 0) {
            continue;
        }

        if (!bidirectional() || (possibleSwaps.count(edge) == 0 && possibleSwaps.count({edge.second, edge.first}) == 0)) {
            possibleSwaps.emplace(edge);
        }
    }

    swaps.clear();
    Node start{};

    // start with identity permutation
    for (unsigned short i = 0; i < nqubits; ++i) {
        start.permutation.emplace(i, i);
    }

    auto                                                             priority = [](const Node& x, const Node& y) { return x.swaps.size() > y.swaps.size(); };
    std::priority_queue<Node, std::vector<Node>, decltype(priority)> queue(priority);
    queue.push(start);

    while (!queue.empty()) {
        Node current = queue.top();
        queue.pop();

        for (const auto& swap: possibleSwaps) {
            Node next = current;
            // continue if the same swap was applied earlier
            if (!next.swaps.empty() && next.swaps.back() == swap)
                continue;

            // apply and insert swap
            std::swap(next.permutation.at(swap.first), next.permutation.at(swap.second));
            next.swaps.emplace_back(swap);
            next.nswaps++;

            bool done = true;
            for (const auto& assignment: goalPermutation) {
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

std::size_t Architecture::getCouplingLimit(const std::set<unsigned short>& qubitChoice) const {
    return findCouplingLimit(getCouplingMap(), getNqubits(), qubitChoice);
}

unsigned long Architecture::bfs(unsigned short start, unsigned short goal, const std::set<Edge>& teleportations) const {
    std::queue<std::vector<int>> queue;
    std::vector<int>             v;
    v.push_back(start);
    queue.push(v);
    std::vector<std::vector<int>> solutions;

    unsigned long length = 0;
    std::set<int> successors;
    while (!queue.empty()) {
        v = queue.front();
        queue.pop();
        int current = v[v.size() - 1];
        if (current == goal) {
            length = v.size();
            solutions.push_back(v);
            break;
        } else {
            successors.clear();
            for (const auto& edge: getCouplingMap()) {
                if (edge.first == current && !contains(v, edge.second)) {
                    successors.insert(edge.second);
                }
                if (edge.second == current && !contains(v, edge.first)) {
                    successors.insert(edge.first);
                }
            }
            for (const auto& edge: teleportations) {
                if (edge.first == current && !contains(v, edge.second)) {
                    successors.insert(edge.second);
                }
                if (edge.second == current && !contains(v, edge.first)) {
                    successors.insert(edge.first); // was v2 but this is probably wrong
                }
            }

            for (int successor: successors) {
                std::vector<int> v2 = v;
                v2.push_back(successor);
                queue.push(v2);
            }
        }
    }
    while (!queue.empty() && queue.front().size() == length) {
        if (queue.front()[queue.front().size() - 1] == goal) {
            solutions.push_back(queue.front());
        }
        queue.pop();
    }

    //TODO: different weight if this contains a teleportation
    for (const auto& s: solutions) {
        for (std::size_t j = 0; j < s.size() - 1; j++) {
            Edge e{s[j], s[j + 1]};
            if (getCouplingMap().find(e) != getCouplingMap().end()) {
                return (length - 2) * 7;
            }
        }
    }

    if (length == 2 && getCouplingMap().find(Edge{start, goal}) == getCouplingMap().end() && getCouplingMap().find(Edge{goal, start}) == getCouplingMap().end()) {
        return 7;
    }

    return (length - 2) * 7 + 4;
}

std::size_t Architecture::findCouplingLimit(const CouplingMap& cm, int nQubits) {
    std::vector<std::vector<unsigned short>> connections;
    std::vector<int>                         d;
    std::vector<bool>                        visited;
    connections.resize(nQubits);
    int maxSum = -1;
    for (auto edge: cm) {
        connections.at(edge.first).emplace_back(edge.second);
    }
    for (int q = 0; q < nQubits; ++q) {
        d.clear();
        d.resize(nQubits);
        std::fill(d.begin(), d.end(), 0);
        visited.clear();
        visited.resize(nQubits);
        std::fill(visited.begin(), visited.end(), false);
        findCouplingLimit(q, 0, connections, d, visited);
        auto it = std::max_element(d.begin(), d.end());
        if ((*it) > maxSum)
            maxSum = (*it);
    }
    return maxSum;
}

std::size_t Architecture::findCouplingLimit(const CouplingMap& cm, int nQubits, const std::set<unsigned short>& qubitChoice) {
    std::vector<std::vector<unsigned short>> connections;
    std::vector<int>                         d;
    std::vector<bool>                        visited;
    connections.resize(nQubits);
    int maxSum = -1;
    for (auto edge: cm) {
        if (qubitChoice.count(edge.first) && qubitChoice.count(edge.second))
            connections.at(edge.first).emplace_back(edge.second);
    }
    for (int q = 0; q < nQubits; ++q) {
        if (connections.at(q).empty())
            continue;
        d.clear();
        d.resize(nQubits);
        std::fill(d.begin(), d.end(), 0);
        visited.clear();
        visited.resize(nQubits);
        std::fill(visited.begin(), visited.end(), false);
        findCouplingLimit(q, 0, connections, d, visited);
        auto it = std::max_element(d.begin(), d.end());
        if ((*it) > maxSum)
            maxSum = (*it);
    }
    return maxSum;
}

void Architecture::findCouplingLimit(unsigned short node, int curSum, const std::vector<std::vector<unsigned short>>& connections, std::vector<int>& d, std::vector<bool>& visited) {
    if (visited.at(node))
        return;
    visited[node] = true;

    if (d.at(node) < curSum)
        d[node] = curSum;
    if (connections.at(node).empty()) {
        visited[node] = false;
        return;
    }

    for (auto child: connections.at(node)) {
        findCouplingLimit(child, curSum + 1, connections, d, visited);
    }

    visited[node] = false;
}

void Architecture::getHighestFidelityCouplingMap(unsigned short subsetSize, CouplingMap& reducedMap) {
    if (!isArchitectureAvailable() || nqubits == subsetSize || properties.empty()) {
        reducedMap = couplingMap;
    } else {
        double bestFidelity        = 0.0;
        auto   allConnectedSubsets = getAllConnectedSubsets(subsetSize);

        for (const auto& qubitChoice: allConnectedSubsets) {
            double      currentFidelity{};
            CouplingMap map{};
            getReducedCouplingMap(qubitChoice, map);
            currentFidelity = getAverageArchitectureFidelity(map, qubitChoice, properties);
            if (currentFidelity > bestFidelity) {
                reducedMap   = map;
                bestFidelity = currentFidelity;
            }
        }
    }
}
std::vector<std::set<unsigned short>> Architecture::getAllConnectedSubsets(unsigned short subsetSize) {
    std::vector<std::set<unsigned short>> result{};
    if (!isArchitectureAvailable() || nqubits == subsetSize) {
        result.emplace_back(getQubitSet());
    } else if (nqubits < subsetSize) {
        throw QMAPException("Architecture too small!");
    } else {
        auto filter = [&](const std::set<unsigned short>& subset) {
            CouplingMap cm = {};
            Architecture::getReducedCouplingMap(subset, cm);
            return isConnected(subset, cm);
        };
        for (const auto& subset: subsets(getQubitSet(), subsetSize, filter)) {
            result.emplace_back(subset);
        }
    }
    return result;
}
void Architecture::getReducedCouplingMaps(unsigned short subsetSize, std::vector<CouplingMap>& couplingMaps) {
    couplingMaps.clear();
    if (!isArchitectureAvailable()) {
        couplingMaps.emplace_back(getFullyConnectedMap(subsetSize));
    } else {
        for (const auto& qubitChoice: getAllConnectedSubsets(subsetSize)) {
            couplingMaps.emplace_back();
            getReducedCouplingMap(qubitChoice, couplingMaps.back());
        }
    }
}
void Architecture::getReducedCouplingMap(const std::set<unsigned short>& qubitChoice, CouplingMap& reducedMap) {
    reducedMap.clear();
    if (!isArchitectureAvailable()) {
        reducedMap = getFullyConnectedMap(qubitChoice.size());
    } else {
        for (const auto& [q0, q1]: couplingMap) {
            if (qubitChoice.find(q0) != qubitChoice.end() && qubitChoice.find(q1) != qubitChoice.end()) {
                reducedMap.emplace(q0, q1);
            }
        }
    }
}

double Architecture::getAverageArchitectureFidelity(const CouplingMap& cm, const std::set<unsigned short>& qubitChoice, const Properties& props) {
    if (props.empty()) {
        return 0.0;
    }
    double result = 1.0;
    for (const auto& [control, target]: cm) {
        if (props.twoQubitErrorRateAvailable(control, target)) {
            result *= (1.0 - props.getTwoQubitErrorRate(control, target));
        }
    }

    for (const auto& qubit: qubitChoice) {
        if (props.singleQubitErrorRate.available(qubit)) {
            result *= (1.0 - props.getAverageSingleQubitErrorRate(qubit));
        }
    }
    return result;
}

std::vector<unsigned short> Architecture::getQubitList(const CouplingMap& couplingMap) {
    std::set<unsigned short> result{};
    for (const auto& edge: couplingMap) {
        result.emplace(edge.first);
        result.emplace(edge.second);
    }
    return {result.begin(), result.end()};
}

bool Architecture::isConnected(const std::set<unsigned short>& qubitChoice, const CouplingMap& reducedCouplingMap) {
    std::set<unsigned short> reachedQubits{};
    reachedQubits.insert(*(qubitChoice.begin()));
    dfs(*(qubitChoice.begin()), reachedQubits, reducedCouplingMap);
    return (reachedQubits == qubitChoice);
}
