//
// This file is part of the MQT QMAP library released under the MIT license.
// See README.md or go to https://github.com/cda-tum/qmap for more information.
//

#include <utils.hpp>

void Dijkstra::buildTable(const std::uint16_t n, const CouplingMap& couplingMap,
                          Matrix& distanceTable,
                          const std::function<double(const Node&)>& cost) {
  distanceTable.clear();
  distanceTable.resize(n, std::vector<double>(n, -1.));

  for (std::uint16_t i = 0; i < n; ++i) {
    std::vector<Dijkstra::Node> nodes(n);
    for (std::uint16_t j = 0; j < n; ++j) {
      nodes.at(j).containsCorrectEdge = false;
      nodes.at(j).visited             = false;
      nodes.at(j).pos                 = j;
      nodes.at(j).cost                = -1.;
    }

    nodes.at(i).cost = 0.;

    dijkstra(couplingMap, nodes, i);

    for (std::uint16_t j = 0; j < n; ++j) {
      if (i == j) {
        distanceTable.at(i).at(j) = 0;
      } else {
        distanceTable.at(i).at(j) = cost(nodes.at(j));
      }
    }
  }
}

void Dijkstra::dijkstra(const CouplingMap& couplingMap,
                        std::vector<Node>& nodes, std::uint16_t start) {
  std::priority_queue<Node*> queue{};
  queue.push(&nodes.at(start));
  while (!queue.empty()) {
    auto* current    = queue.top();
    current->visited = true;
    queue.pop();
    auto pos = current->pos;

    for (const auto& edge : couplingMap) {
      std::optional<std::uint16_t> to          = std::nullopt;
      bool                         correctEdge = false;
      if (pos == edge.first) {
        to          = edge.second;
        correctEdge = true;
      } else if (pos == edge.second) {
        to = edge.first;
      }
      if (to.has_value()) {
        if (nodes.at(*to).visited) {
          continue;
        }

        Node newNode;
        newNode.cost                = current->cost + 1.0;
        newNode.pos                 = to;
        newNode.containsCorrectEdge = correctEdge;
        if (nodes.at(*to).cost < 0 || newNode < nodes.at(*to)) {
          nodes.at(*to) = newNode;
          queue.push(&nodes.at(*to));
        }
      }
    }
  }
}

/// Create a string representation of a given permutation
/// \param pi permutation
/// \return string representation of pi
std::string printPi(std::vector<std::uint16_t>& pi) {
  if (std::is_sorted(pi.begin(), pi.end())) {
    return "( )";
  }

  std::stringstream perm{};
  perm << '(';
  for (std::uint64_t i = 0; i < pi.size() - 1; i++) {
    perm << pi[i] << ',';
  }
  perm << pi[pi.size() - 1] << ')';

  return perm.str();
}

/// Simple depth-first-search implementation used to check whether a given
/// subset of qubits is connected on the given architecture \param current index
/// of current qubit \param visited visited qubits \param cm coupling map of
/// architecture
void dfs(std::uint16_t current, std::set<std::uint16_t>& visited,
         const CouplingMap& rcm) {
  for (auto edge : rcm) {
    if (edge.first == current) {
      if (visited.count(edge.second) == 0U) {
        visited.insert(edge.second);
        dfs(edge.second, visited, rcm);
      }
    } else if (edge.second == current) {
      if (visited.count(edge.first) == 0U) {
        visited.insert(edge.first);
        dfs(edge.first, visited, rcm);
      }
    }
  }
}

std::vector<QubitSubset> subsets(const QubitSubset&     input,
                                 const std::size_t      size,
                                 const filter_function& filter) {
  const std::size_t        n = input.size();
  std::vector<QubitSubset> result{};

  if (size == 0) {
    throw std::invalid_argument("Length of subset must be greater than 0");
  }

  if (size > n) {
    throw std::invalid_argument("Length of subset must be less than or equal "
                                "to the size of the input set");
  }

  if (size == 1) {
    for (const auto& item : input) {
      result.emplace_back();
      result.back().emplace(item);
    }
  } else {
    std::uint64_t i = (1U << size) - 1U;

    while ((i >> n) == 0U) {
      std::set<std::uint16_t> v{};
      auto                    it = input.begin();

      for (std::size_t j = 0U; j < n; j++, ++it) {
        if ((i & (1U << j)) != 0U) {
          v.emplace(*it);
        }
      }
      if (filter == nullptr || filter(v)) {
        result.emplace_back(v);
      }

      // this computes the lexicographical next bitset from a set, see
      // https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
      const std::uint64_t t = (i | (i - 1)) + 1;
      i                     = t | ((((t & -t) / (i & -i)) >> 1) - 1);
    }
  }

  return result;
}

void parseLine(const std::string& line, char separator,
               const std::set<char>&     escapeChars,
               const std::set<char>&     ignoredChars,
               std::vector<std::string>& result) {
  result.clear();
  std::string word;
  bool        inEscape = false;
  for (const char c : line) {
    if (ignoredChars.find(c) != ignoredChars.end()) {
      continue;
    }
    if (inEscape) {
      if (escapeChars.find(c) != escapeChars.end()) {
        inEscape = false;
      } else {
        word += c;
      }
    } else {
      if (escapeChars.find(c) != escapeChars.end()) {
        inEscape = true;
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

CouplingMap getFullyConnectedMap(const std::uint16_t nQubits) {
  CouplingMap result{};
  for (std::uint16_t q = 0; q < nQubits; ++q) {
    for (std::uint16_t p = q + 1; p < nQubits; ++p) {
      result.emplace(q, p);
      result.emplace(p, q);
    }
  }
  return result;
}
