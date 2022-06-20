/*
 * This file is part of the MQT QMAP library which is released under the MIT license.
 * See file README.md or go to https://iic.jku.at/eda/research/ibm_qx_mapping/ for more information.
 */

#include <utils.hpp>

void Dijkstra::build_table(unsigned short n, const std::set<Edge>& couplingMap, Matrix& distanceTable, const std::function<double(const Node&)>& cost) {
    distanceTable.clear();
    distanceTable.resize(n, std::vector<double>(n, -1.));

    for (unsigned short i = 0; i < n; ++i) {
        std::vector<Dijkstra::Node> nodes(n);
        for (unsigned short j = 0; j < n; ++j) {
            nodes.at(j).contains_correct_edge = false;
            nodes.at(j).visited               = false;
            nodes.at(j).pos                   = j;
            nodes.at(j).cost                  = -1.;
        }

        nodes.at(i).cost = 0.;

        dijkstra(couplingMap, nodes, i);

        if (VERBOSE) {
            for (const auto& node: nodes) {
                std::cout << node.cost << " ";
            }
            std::cout << std::endl;
        }

        for (int j = 0; j < n; ++j) {
            if (i == j) {
                distanceTable.at(i).at(j) = 0;
            } else {
                distanceTable.at(i).at(j) = cost(nodes.at(j));
            }
        }
    }
}

void Dijkstra::dijkstra(const CouplingMap& couplingMap, std::vector<Node>& nodes, unsigned short start) {
    std::priority_queue<Node*> queue{};
    queue.push(&nodes.at(start));
    while (!queue.empty()) {
        auto current     = queue.top();
        current->visited = true;
        queue.pop();
        auto pos = current->pos;

        for (const auto& edge: couplingMap) {
            short to          = -1;
            bool  correctEdge = false;
            if (pos == edge.first) {
                to          = edge.second;
                correctEdge = true;
            } else if (pos == edge.second) {
                to = edge.first;
            }
            if (to != -1) {
                if (nodes.at(to).visited)
                    continue;

                Node new_node;
                new_node.cost                  = current->cost + 1.0;
                new_node.pos                   = to;
                new_node.contains_correct_edge = correctEdge;
                if (nodes.at(to).cost < 0 || new_node < nodes.at(to)) {
                    nodes.at(to) = new_node;
                    queue.push(&nodes.at(to));
                }
            }
        }
    }
}

/// Create a string representation of a given permutation
/// \param pi permutation
/// \return string representation of pi
std::string printPi(std::vector<unsigned short>& pi) {
    if (std::is_sorted(pi.begin(), pi.end())) {
        return "( )";
    }

    std::stringstream perm{};
    perm << '(';
    for (unsigned long i = 0; i < pi.size() - 1; i++) {
        perm << pi[i] << ',';
    }
    perm << pi[pi.size() - 1] << ')';

    return perm.str();
}

/// Simple depth-first-search implementation used to check whether a given subset of qubits is
/// connected on the given architecture
/// \param current index of current qubit
/// \param visited visited qubits
/// \param cm coupling map of architecture
void dfs(unsigned short current, std::set<unsigned short>& visited, const CouplingMap& rcm) {
    for (auto edge: rcm) {
        if (edge.first == current) {
            if (!visited.count(edge.second)) {
                visited.insert(edge.second);
                dfs(edge.second, visited, rcm);
            }
        } else if (edge.second == current) {
            if (!visited.count(edge.first)) {
                visited.insert(edge.first);
                dfs(edge.first, visited, rcm);
            }
        }
    }
}

/// Helper function returning correct 1D array index for 3D array
/// \param k first index
/// \param i second index
/// \param j third index
/// \return index in 1D array
unsigned long idx(unsigned int k, unsigned short i, unsigned short j, const std::set<unsigned short>& iValues, const std::set<unsigned short>& jValues) {
    unsigned short counti = 0;
    for (unsigned short iVal: iValues) {
        if (iVal == i) break;
        counti++;
    }
    unsigned short countj = 0;
    for (unsigned short jVal: jValues) {
        if (jVal == j) break;
        countj++;
    }

    return k * jValues.size() * iValues.size() + counti * jValues.size() + countj;
}

unsigned long idx(unsigned int k, unsigned short i, unsigned short j, const std::set<unsigned short>& iValues, unsigned short nj) {
    unsigned short counti = 0;
    for (unsigned short iVal: iValues) {
        if (iVal == i) break;
        counti++;
    }

    return k * static_cast<std::size_t>(nj) * iValues.size() + counti * static_cast<std::size_t>(nj) + j;
}

std::vector<std::set<unsigned short>>
subsets(const std::set<unsigned short>& input, int length, filter_function filter) {
    std::size_t                           n = input.size();
    std::vector<std::set<unsigned short>> result;

    if (length == 1) {
        for (const auto& item: input) {
            result.emplace_back();
            result.back().emplace(item);
        }
    } else {
        std::size_t i = (1U << length) - 1U;

        while (!(i >> n)) {
            std::set<unsigned short> v{};
            auto                     it = input.begin();

            for (std::size_t j = 0U; j < n; j++, ++it) {
                if (i & (1U << j)) {
                    v.emplace(*it);
                }
            }
            if (filter == nullptr || filter(v)) {
                result.emplace_back(v);
            }

            //this computes the lexographical next bitset from a set.
            //the unsigned int t = v | (v - 1); // t gets v's least significant 0 bits set to 1
            //// Next set to 1 the most significant bit to change,
            //// set to 0 the least significant ones, and add the necessary 1 bits.
            //w = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1))
            // is the original, which involves counting the leading zeros via __builtin_ctz, the version below
            // uses division to counteract that problem and might be slower on architectures that have a fast
            // variant of ctz, but more convenient on others
            i = (i + (i & (-i))) | (((i ^ (i + (i & (-i)))) >> 2) / (i & (-i)));
        }
    }

    return result;
}

void parse_line(const std::string& line, char separator, const std::set<char>& escape_chars,
                const std::set<char>& ignored_chars, std::vector<std::string>& result) {
    std::string word;
    bool        in_escape = false;
    for (char c: line) {
        if (ignored_chars.find(c) != ignored_chars.end()) {
            continue;
        }
        if (in_escape) {
            if (escape_chars.find(c) != escape_chars.end()) {
                in_escape = false;
            } else {
                word += c;
            }
        } else {
            if (escape_chars.find(c) != escape_chars.end()) {
                in_escape = true;
            } else if (c == separator) {
                result.push_back(word);
                word = "";
            } else {
                word += c;
            }
        }
    }
    result.push_back(word);
}

std::set<std::pair<unsigned short, unsigned short>>
getFullyConnectedMap(unsigned short nQubits) {
    std::set<std::pair<unsigned short, unsigned short>> result{};
    for (int q = 0; q < nQubits; ++q) {
        for (int p = q + 1; p < nQubits; ++p) {
            result.emplace(q, p);
            result.emplace(p, q);
        }
    }
    return result;
}

std::string escapeChars(const std::string& s, const std::string& chars) {
    std::stringstream ss;
    for (auto c: s) {
        if (chars.find(c) != std::string::npos) {
            ss << "\\" << c;
        } else if (c == '\n') {
            ss << "\\n";
        } else if (c == '\t') {
            ss << "\\t";
        } else {
            ss << c;
        }
    }
    return ss.str();
}