/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include "Architecture.hpp"

void Architecture::loadCouplingMap(AvailableArchitecture architecture) {
    std::stringstream ss{getCouplingMapSpecification(architecture)};
    architectureName = toString(architecture);
    loadCouplingMap(ss);
}

void Architecture::loadCouplingMap(std::istream& is) {
    loadCouplingMap(std::move(is));
}

void Architecture::loadCouplingMap(const std::string& filename) {
    size_t slash     = filename.find_last_of('/');
    size_t dot       = filename.find_last_of('.');
    architectureName = filename.substr(slash + 1, dot - slash - 1);
    auto ifs         = std::ifstream(filename);
    if (ifs.good())
        this->loadCouplingMap(ifs);
    else
        throw QMAPException("Error opening coupling map file.");
}

void Architecture::loadCouplingMap(std::istream&& is) {
    couplingMap.clear();
    calibrationData.clear();
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
    calibrationData.clear();
    architectureName = "generic_" + std::to_string(nQ);
    createDistanceTable();
}

void Architecture::loadCalibrationData(std::istream& is) {
    loadCalibrationData(std::move(is));
}

void Architecture::loadCalibrationData(const std::string& filename) {
    size_t slash    = filename.find_last_of('/');
    size_t dot      = filename.find_last_of('.');
    calibrationName = filename.substr(slash + 1, dot - slash - 1);
    auto ifs        = std::ifstream(filename);
    if (ifs.good())
        this->loadCalibrationData(ifs);
    else
        throw QMAPException("Error opening calibration data file.");
}

void Architecture::loadCalibrationData([[maybe_unused]] std::istream&& is) {
    calibrationData.clear();
    calibrationData.reserve(nqubits);
    // TODO: Read in IBM csv Format into calibrationData
    // Qubit,T1 (µs),T2 (µs),Frequency (GHz),Readout error,Single-qubit U2 error rate,CNOT error rate,Date
    // Q0,43.113390027720456,80.43170304899623,5.2540548421267825,3.0000000000000027e-2,6.471748738850555e-4,"cx0_1: 1.286e-2",Fri Feb 28 2020 09:27:24 GMT+0100 (CET)
    // Q1,57.036278779745814,44.25327200058967,5.04877069952246,5.499999999999994e-2,6.227987234977298e-4,"cx1_0: 1.286e-2, cx1_2: 1.117e-2, cx1_3: 1.006e-2",Fri Feb 28 2020 09:27:24 GMT+0100 (CET)
    // Q2,81.96027845904459,87.47080521889637,5.230569177436213,2.1666666666666723e-2,5.16080288896023e-4,"cx2_1: 1.117e-2",Fri Feb 28 2020 09:27:24 GMT+0100 (CET)
    // Q3,67.95013744985455,86.03784184638282,5.200948312700143,1.8333333333333313e-2,3.530411094454146e-4,"cx3_1: 1.006e-2, cx3_4: 1.670e-2",Fri Feb 28 2020 09:27:24 GMT+0100 (CET)
    // Q4,62.156360190670306,3.731765939150933,5.065777885569538,1.5000000000000013e-2,1.106360359568795e-3,"cx4_3: 1.670e-2",Fri Feb 28 2020 09:27:24 GMT+0100 (CET)

    createFidelityTable();
}

void Architecture::loadCalibrationData(const std::vector<CalibrationData>& calData) {
    calibrationData  = calData;
    architectureName = "generic_" + std::to_string(nqubits);
    createFidelityTable();
}

Architecture::Architecture(unsigned short nQ, const CouplingMap& couplingMap) {
    loadCouplingMap(nQ, couplingMap);
}

Architecture::Architecture(unsigned short nQ, const CouplingMap& couplingMap, const std::vector<CalibrationData>& calibrationData):
    Architecture(nQ, couplingMap) {
    loadCalibrationData(calibrationData);
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
    fidelityTable.resize(nqubits, std::vector<double>(nqubits, 1.0));

    singleQubitFidelities.resize(nqubits, 1.0);

    for (const auto& qubit: calibrationData) {
        for (const auto& entry: qubit.cnotErrorRate) {
            fidelityTable.at(entry.first.first).at(entry.first.second) -= entry.second;
        }
        singleQubitFidelities.at(qubit.qubit) -= qubit.singleQubitErrorRate;
    }
}

unsigned long Architecture::minimumNumberOfSwaps(std::vector<unsigned short>& permutation) {
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
