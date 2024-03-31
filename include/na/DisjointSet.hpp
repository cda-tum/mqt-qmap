#pragma once

#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <valarray>
#include <vector>

template <class T> struct DisjointSet {
  std::unordered_map<T, T>           parent;
  std::unordered_map<T, std::size_t> rank;

  template <class Iterator>
  explicit DisjointSet(const Iterator& begin, const Iterator& end) {
    std::for_each(begin, end, [&](const auto& element) {
      parent[element] = element;
      rank[element]   = 0;
    });
  }

  T findSet(T v) {
    if (parent[v] != v) {
      parent[v] = findSet(parent[v]);
    }
    return parent[v];
  }

  void unionSet(T x, T y) {
    x = findSet(x);
    y = findSet(y);
    if (rank[x] > rank[y]) {
      parent[y] = x;
    } else {
      parent[x] = y;
      if (rank[x] == rank[y]) {
        rank[y]++;
      }
    }
  }
};