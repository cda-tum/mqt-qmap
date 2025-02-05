#pragma once

#include "na/NAComputation.hpp"
#include "na/azac/CompilerBase.hpp"
#include "na/operations/NALocalOperation.hpp"
#include "na/operations/NAOperation.hpp"
#include "na/operations/NAShuttlingOperation.hpp"

#include <algorithm>
#include <chrono>
#include <type_traits>
#include <utility>

namespace na {

template <typename T> class Router {
private:
  /// constant, the distance of AOD row and col to some trap. We use 1um here.
  static constexpr std::size_t PARKING_DIST = 1;

  std::priority_queue<std::pair<std::size_t, std::size_t>, std::vector<std::pair<std::size_t, std::size_t>>, std::greater<std::pair<std::size_t, std::size_t>>>
      aod_end_time;
  std::vector<std::size_t> aod_dependency;
  std::vector<std::size_t> rydberg_dependency;
  std::vector<std::size_t> qubit_dependency;
  std::unordered_map<std::tuple<const SLM*, std::size_t, std::size_t>, std::size_t> site_dependency;

protected:
  /// generate rearrangement layers between two Rydberg layers
  auto route_qubit() -> void {
    for (std::size_t i = 0;
         i < static_cast<T*>(this)->architecture.dict_AOD.size(); ++i) {
      aod_end_time.emplace(0, i);
    }
    aod_dependency(
        static_cast<T*>(this)->architecture.dict_AOD.size(), 0);
    rydberg_dependency(
        static_cast<T*>(this)->architecture.entanglement_zone.size(), 0);
    std::chrono::microseconds time_mis{};
    qubit_dependency(static_cast<T*>(this)->n_q, 0);
    write_initial_instruction();

    for (std::size_t layer = 0;
         layer < static_cast<T*>(this)->gate_scheduling.size(); ++layer) {
      // extract sets of movement that can be performed simultaneously
      const auto t_s = std::chrono::system_clock::now();
      route_qubit_mis(layer);
      time_mis += (std::chrono::system_clock::now() - t_s);
      std::cout << "[INFO] ZAC: Solve for Rydberg stage " << (layer + 1) << "/"
                << static_cast<T*>(this)->gate_scheduling.size()
                << ". mis time=" << time_mis << "\n";
    }
    static_cast<T*>(this)->runtime_analysis["routing"] = time_mis;
  }

  /// process layers of movement from storage zone to Rydberg and back to
  /// storage zone
  auto route_qubit_mis(const std::size_t layer) -> void {
    const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
        initial_mapping = static_cast<T*>(this)->qubit_mapping[2 * layer];
    const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
        gate_mapping = static_cast<T*>(this)->qubit_mapping[2 * layer + 1];
    const std::optional<
        std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>&
        final_mapping =
            layer + 2 < static_cast<T*>(this)->qubit_mapping.size()
                ? static_cast<T*>(this)->qubit_mapping[2 * layer + 2]
                : std::nullopt;

    // sort remain_graph based on qubit distance if using maximal is
    std::vector<std::size_t> remain_graph; // consist qubits to be moved
    for (const std::pair<qc::Qubit, qc::Qubit>* gate :
         static_cast<T*>(this)->gate_scheduling[layer]) {
      if (initial_mapping[gate->first] != gate_mapping[gate->first]) {
        remain_graph.emplace_back(gate->first);
      }
      if (initial_mapping[gate->second] != gate_mapping[gate->second]) {
        remain_graph.emplace_back(gate->second);
      }
    }

    if (static_cast<T*>(this)->routing_strategy !=
            CompilerBase<T>::RoutingStrategy::MAXIMAL_IS) {
      std::sort(remain_graph.begin(), remain_graph.end(),
                [&](const std::size_t a, const std::size_t b) {
                  const auto [a_init_x, a_init_y] = std::apply(
                      static_cast<T*>(this)->architecture.exact_SLM_location,
                      initial_mapping[a]);
                  const auto [b_init_x, b_init_y] = std::apply(
                      static_cast<T*>(this)->architecture.exact_SLM_location,
                      initial_mapping[b]);
                  const auto [a_gate_x, a_gate_y] = std::apply(
                      static_cast<T*>(this)->architecture.exact_SLM_location,
                      gate_mapping[a]);
                  const auto [b_gate_x, b_gate_y] = std::apply(
                      static_cast<T*>(this)->architecture.exact_SLM_location,
                      gate_mapping[b]);
                  const auto a_dx = static_cast<double>(a_gate_x) -
                                    static_cast<double>(a_init_x);
                  const auto a_dy = static_cast<double>(a_gate_y) -
                                    static_cast<double>(a_init_y);
                  const auto b_dx = static_cast<double>(b_gate_x) -
                                    static_cast<double>(b_init_x);
                  const auto b_dy = static_cast<double>(b_gate_y) -
                                    static_cast<double>(b_init_y);
                  return a_dx * a_dx + a_dy * a_dy > b_dx * b_dx + b_dy * b_dy;
                });
    }
    const auto id_layer_start =
        static_cast<T*>(this)->result.instructions.size();
    std::size_t batch = 0;
    while (!remain_graph.empty()) {
      // graph construction
      const auto& vectors =
          graph_construction(remain_graph, initial_mapping, gate_mapping);
      // collect violation
      const auto& violations = collect_violation(vectors);
      // solve MIS
      const auto& moved_qubits = maximalis_solve(vectors.size(), violations);

      std::unordered_set<std::size_t>
          set_aod{}; // use to record aods per movement layer
      set_aod.reserve(moved_qubits.size());
      for (const auto i : moved_qubits) {
        set_aod.emplace(remain_graph[i]);
      }
      process_movement_layer(set_aod, initial_mapping, gate_mapping);
      std::vector<std::size_t> tmp;
      for (const auto q : remain_graph) {
        if (set_aod.find(q) == set_aod.end()) {
          tmp.emplace_back(q);
        }
      }
      remain_graph = tmp;
      batch = batch + 1;
    }

    // append a layer for gate execution
    process_gate_layer(layer, gate_mapping);
    // move qubit back to the final location
    if (!final_mapping) {
      if (static_cast<T*>(this)->dynamic_placement ||
          static_cast<T*>(this)->reuse) {
        remain_graph.clear(); // consist qubits to be moved
        for (const std::pair<qc::Qubit, qc::Qubit>* gate :
             static_cast<T*>(this)->gate_scheduling[layer]) {
          if ((*final_mapping)[gate->first] != gate_mapping[gate->first]) {
            remain_graph.emplace_back(gate->first);
          }
          if ((*final_mapping)[gate->second] != gate_mapping[gate->second]) {
            remain_graph.emplace_back(gate->second);
          }
        }

        if (static_cast<T*>(this)->routing_strategy !=
                CompilerBase<T>::RoutingStrategy::MAXIMAL_IS) {
        std::sort(remain_graph.begin(), remain_graph.end(),
               [&](const std::size_t a, const std::size_t b) {
                 const auto [a_final_x, a_final_y] = std::apply(
                     static_cast<T*>(this)->architecture.exact_SLM_location,
                     (*final_mapping)[a]);
                 const auto [b_final_x, b_final_y] = std::apply(
                     static_cast<T*>(this)->architecture.exact_SLM_location,
                     (*final_mapping)[b]);
                 const auto [a_gate_x, a_gate_y] = std::apply(
                     static_cast<T*>(this)->architecture.exact_SLM_location,
                     gate_mapping[a]);
                 const auto [b_gate_x, b_gate_y] = std::apply(
                     static_cast<T*>(this)->architecture.exact_SLM_location,
                     gate_mapping[b]);
                 const auto a_dx = static_cast<double>(a_gate_x) -
                                   static_cast<double>(a_final_x);
                 const auto a_dy = static_cast<double>(a_gate_y) -
                                   static_cast<double>(a_final_y);
                 const auto b_dx = static_cast<double>(b_gate_x) -
                                   static_cast<double>(b_final_x);
                 const auto b_dy = static_cast<double>(b_gate_y) -
                                   static_cast<double>(b_final_y);
                 return a_dx * a_dx + a_dy * a_dy > b_dx * b_dx + b_dy * b_dy;
               });
        }
        while (!remain_graph.empty()) {
          // graph construction
          const auto& vectors =
              graph_construction(remain_graph, *final_mapping, gate_mapping);
          // collect violation
          const auto& violations = collect_violation(vectors);

          const auto& moved_qubits = maximalis_solve(vectors.size(), violations);
          // todo: add layer
          std::unordered_set<std::size_t>
              set_aod; // use to record aods per movement layer
          set_aod.reserve(moved_qubits.size());
          for (const auto i : moved_qubits) {
            set_aod.emplace(remain_graph[i]);
          }

          process_movement_layer(set_aod, gate_mapping, final_mapping);

          std::vector<std::size_t> tmp;
          for (const auto q : remain_graph) {
            if (set_aod.find(q) == set_aod.end()) {
              tmp.emplace_back(q);
            }
          }
          remain_graph = tmp;
          ++batch;
        }
      } else {
        // construct reverse layer
        construct_reverse_layer(id_layer_start, gate_mapping, final_mapping);
      }
      aod_assignment(id_layer_start);
    }
  }

  auto graph_construction(
      const std::vector<std::size_t>& remainGraph,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& initialMapping,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& finalMapping)
      -> std::vector<
          std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> {
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
        vectors{};
    const std::size_t vector_length =
        static_cast<T*>(this)->use_window
            ? std::min(static_cast<T*>(this)->window_size, remainGraph.size())
            : remainGraph.size();
    for (std::size_t i = 0; i < vector_length; ++i) {
      const auto q = remainGraph[i];
      const auto& [q_x, q_y] =
          std::apply(static_cast<T*>(this)->architecture.exact_SLM_location,
                     initialMapping[q]);
      const auto& [site_x, site_y] =
          std::apply(static_cast<T*>(this)->architecture.exact_SLM_location,
                     finalMapping[q]);
      vectors.emplace_back(q_x, site_x, q_y, site_y);
    }
    return vectors;
  }

  auto collect_violation(
      std::vector<std::vector<
          std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>>
          vectors) -> std::vector<std::pair<std::size_t, std::size_t>> {
    std::vector<std::pair<std::size_t, std::size_t>> violations{};
    for (std::size_t i = 0; i < vectors.size(); ++i) {
      for (std::size_t j = i + 1; j < vectors.size(); ++j) {
        if (!compatible_2D(vectors[i], vectors[j])) {
          violations.emplace_back(i, j);
        }
      }
    }
    return violations;
  }

  /// solve maximal independent set
  auto maximalis_solve(const std::size_t n,
                       const std::vector<std::pair<std::size_t, std::size_t>>& edges)
      -> std::vector<std::size_t> {
    // assume the vertices are sorted based on qubit distance
    std::vector is_node_conflict(n, false);
    std::vector node_neighbors(n, std::vector<std::size_t>{});
    for (const auto& edge : edges) {
      node_neighbors[edge.first].emplace_back(edge.second);
      node_neighbors[edge.second].emplace_back(edge.first);
    }
    std::vector<std::size_t> result{};
    for (std::size_t i = 0; i < is_node_conflict.size(); ++i) {
      if (!is_node_conflict[i]) {
        result.emplace_back(i);
        for (const auto j : node_neighbors[i]) {
          is_node_conflict[j] = true;
        }
      }
    }
    return result;
  }

  /// check if move a and b can be performed simultaneously
  auto compatible_2D(
      std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> a,
      std::tuple<std::size_t, std::size_t, std::size_t, std::size_t> b)
      -> bool {
    if (std::get<0>(a) == std::get<0>(b) && std::get<1>(a) != std::get<1>(b)) {
      return true;
    }
    if (std::get<1>(a) == std::get<1>(b) && std::get<0>(a) != std::get<0>(b)) {
      return true;
    }
    if (std::get<0>(a) < std::get<0>(b) && std::get<1>(a) >= std::get<1>(b)) {
      return true;
    }
    if (std::get<0>(a) > std::get<0>(b) && std::get<1>(a) <= std::get<1>(b)) {
      return true;
    }
    if (std::get<2>(a) == std::get<2>(b) && std::get<3>(a) != std::get<3>(b)) {
      return true;
    }
    if (std::get<3>(a) == std::get<3>(b) && std::get<2>(a) != std::get<2>(b)) {
      return true;
    }
    if (std::get<2>(a) < std::get<2>(b) && std::get<3>(a) >= std::get<3>(b)) {
      return true;
    }
    if (std::get<2>(a) > std::get<2>(b) && std::get<3>(a) <= std::get<3>(b)) {
      return true;
    }
    return true;
  }

  auto write_initial_instruction() -> void {
    static_cast<T*>(this)->result.instructions.clear();
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
        initial_locs{};
    initial_locs.reserve(static_cast<T*>(this)->n_q);
    for (std::size_t i = 0; i < static_cast<T*>(this)->n_q; ++i) {
      initial_locs.emplace_back(
          i, std::get<0>(static_cast<T*>(this)->qubit_mapping[0][i])->id,
          std::get<1>(static_cast<T*>(this)->qubit_mapping[0][i]),
          std::get<2>(static_cast<T*>(this)->qubit_mapping[0][i]));
    }
    static_cast<T*>(this)->result.instructions.emplace_back(
        nlohmann::json{{"type", "init"},
         {"id", 0},
         {"begin_time", 0},
         {"end_time", 0},
         {"init_locs", initial_locs}});

    // process single-qubit gates
    std::unordered_set<std::size_t> set_qubit_dependency{};
    const std::size_t inst_idx =
        static_cast<T*>(this)->result.instructions.size();
    const std::vector<qc::StandardOperation>& list_1q_gate =
        static_cast<T*>(this)->dict_g_1q_parent[nullptr];
    std::vector<nlohmann::json> result_gate{};
    for (const auto& gate_info : list_1q_gate) {
      // collect qubit dependency
      const auto qubit = gate_info.getTargets().front();
      set_qubit_dependency.emplace(qubit_dependency[qubit]);
      qubit_dependency[qubit] = inst_idx;
      result_gate.emplace_back(
          nlohmann::json{{"name", gate_info.getName()}, {"q", qubit}});
    }
    nlohmann::json dependency = {"qubit", std::vector<std::size_t>()};
    dependency["qubit"] = set_qubit_dependency;
    if (!result_gate.empty()) {
      write_1q_gate_instruction(inst_idx, result_gate, dependency, static_cast<T*>(this)->qubit_mapping[0]);
      static_cast<T*>(this)->result.instructions.back()["begin_time"] = 0;
      static_cast<T*>(this)->result.instructions.back()["end_time"] = static_cast<T*>(this)->architecture.time_1qGate * result_gate.size(); // due to sequential execution
    }
  }

        /// generate layers for row-by-row based atom transfer
    auto process_movement_layer(const std::unordered_set<std::size_t>& set_aod_qubit,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& initial_mapping,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& final_mapping) -> void {
      // seperate qubits in list_aod_qubit into multiple lists where qubits in one list can pick up simultaneously
      // we use row-based pick up
      std::unordered_map<std::size_t, std::vector<std::size_t>> pickup_dict; // key: row, value: a list of qubit in the same row
      for (const auto q : set_aod_qubit) {
        const auto& [x, y] = std::apply(static_cast<T*>(this)->architecture.exact_SLM_location, initial_mapping[q]);
        if (pickup_dict.find(y) == pickup_dict.end()) {
          pickup_dict.emplace(y, std::vector<std::size_t>{});
        }
        pickup_dict[y].emplace_back(q);
      }
      std::vector<std::vector<std::size_t>> list_aod_qubits{};
      std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>> list_end_location{};
      std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>> list_begin_location{};
      nlohmann::json dependency = {{"qubit", std::vector<std::size_t>()}, {"site", std::vector<std::size_t>()}};
      // process aod dependency
      const std::size_t inst_idx = static_cast<T*>(this)->result.instructions.size();

      std::unordered_set<std::size_t> set_qubit_dependency;
      std::unordered_set<std::size_t> set_site_dependency;
      for (const auto& [dict_key, dict_value] : pickup_dict) {
        // collect set of aod qubits to pick up
        list_aod_qubits.emplace_back(dict_value);
        std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>> row_begin_location{};
        std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>> row_end_location{};
        for (const auto q : dict_value) {
          // collect qubit begin location
          row_begin_location.emplace_back(q, std::get<0>(initial_mapping[q]), std::get<1>(initial_mapping[q]), std::get<2>(initial_mapping[q]));
          // collect qubit end location
          row_end_location.emplace_back(q, std::get<0>(final_mapping[q]), std::get<1>(final_mapping[q]), std::get<2>(final_mapping[q]));
          const auto& site_key = final_mapping[q];
          if (site_dependency.find(site_key) != site_dependency.end()) {
            set_site_dependency.emplace(site_dependency[site_key]);
          }
          site_dependency[initial_mapping[q]] = inst_idx;

          // collect qubit dependency
          set_qubit_dependency.emplace(qubit_dependency[q]);
          qubit_dependency[q] = inst_idx;
        }
        list_begin_location.emplace_back(std::move(row_begin_location));
        list_end_location.emplace_back(std::move(row_end_location));
      }
      dependency["qubit"] = set_qubit_dependency;
      dependency["site"] = set_site_dependency;
      write_rearrangement_instruction(inst_idx, list_aod_qubits, list_begin_location, list_end_location, dependency);
    }

    auto write_rearrangement_instruction(const std::size_t inst_idx, const std::vector<std::vector<size_t>>& aod_qubits,
      const std::vector<std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>& begin_location,
      const std::vector<std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>& end_location,
      nlohmann::json dependency) -> void {
      std::vector<std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>> begin_location_id{};
      begin_location_id.reserve(begin_location.size());
        for (const auto& row : begin_location) {
          std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> row_tuples{};
          row_tuples.reserve(row.size());
          for (const auto& qubit : row) {
            row_tuples.emplace_back(std::get<0>(qubit), std::get<1>(qubit)->id, std::get<2>(qubit), std::get<3>(qubit));
          }
          begin_location_id.emplace_back(std::move(row_tuples));
        }
      std::vector<std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>> end_location_id{};
      end_location_id.reserve(end_location.size());
      for (const auto& row : end_location) {
        std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>> row_tuples{};
        row_tuples.reserve(row.size());
        for (const auto& qubit : row) {
          row_tuples.emplace_back(std::get<0>(qubit), std::get<1>(qubit)->id, std::get<2>(qubit), std::get<3>(qubit));
        }
        end_location_id.emplace_back(std::move(row_tuples));
      }
      nlohmann::json inst = {
                {"type", "rearrangeJob"},
                {"id", inst_idx},
                {"aod_id", -1},
                {"aod_qubits", aod_qubits},
                {"begin_locs", begin_location_id},
                {"end_locs", end_location_id},
                {"dependency", dependency}
            };
      inst["insts"] = expand_arrangement(inst, begin_location, end_location);
      static_cast<T*>(this)->result.instructions.emplace_back(inst);
    }

  auto flatten_rearrangment_instruction() -> void {
    for (nlohmann::json& inst : static_cast<T*>(this)->result.instructions) {
      if (inst["type"] == "rearrangeJob") {
        nlohmann::json flattened{};
        for (const auto& row : inst["aod_qubits"]) {
          flattened.insert(flattened.end(), row.cbegin(), row.cend());
        }
        inst["aod_qubits"] = flattened;
        flattened.clear();
        for (const auto& row : inst["begin_locs"]) {
          flattened.insert(flattened.end(), row.cbegin(), row.cend());
        }
        inst["begin_locs"] = flattened;
        flattened.clear();
        for (const auto& row : inst["end_locs"]) {
          flattened.insert(flattened.end(), row.cbegin(), row.cend());
        }
        inst["end_locs"] = flattened;
      }
    }
  }

  /// generate a layer for gate execution
  auto process_gate_layer(
      const std::size_t layer,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& gate_mapping)
      -> void {
    const std::vector<std::size_t>& list_gate_idx =
        static_cast<T*>(this)->gate_scheduling_idx[layer];
    const std::vector<const std::pair<qc::Qubit, qc::Qubit>*>& list_gate =
        static_cast<T*>(this)->gate_scheduling[layer];
    const std::vector<const qc::StandardOperation*>& list_1q_gate =
        static_cast<T*>(this)->gate_1q_scheduling[layer];
    std::unordered_map<std::size_t, std::vector<std::size_t>> dict_gate_zone{};
    for (std::size_t i = 0; i < list_gate.size(); ++i) {
      const auto* slm_idx = std::get<0>(gate_mapping[list_gate[i]->first]);
      const auto zone_idx = slm_idx->entanglement_id->front()->id;
      if (dict_gate_zone.find(zone_idx) == dict_gate_zone.end()) {
        dict_gate_zone.emplace(zone_idx, std::vector<std::size_t>{});
      }
      dict_gate_zone[zone_idx].emplace_back(i);
    }
    for (const auto& [rydberg_idx, gate_idxs] : dict_gate_zone) {
      std::vector<nlohmann::json> result_gate{};
      for (const auto i : gate_idxs) {
        result_gate.emplace_back(nlohmann::json{{"id", list_gate_idx[i]},
                                                {"q0", list_gate[i]->first},
                                                {"q1", list_gate[i]->second}});
      }
      std::unordered_set<std::size_t> set_qubit_dependency{};
      const std::size_t inst_idx =
          static_cast<T*>(this)->result.instructions.size();
      for (const auto gate_idx : dict_gate_zone[rydberg_idx]) {
        const auto gate = list_gate[gate_idx];
        // collect qubit dependency
        set_qubit_dependency.emplace(qubit_dependency[gate->first]);
        qubit_dependency[gate->first] = inst_idx;
        set_qubit_dependency.emplace(qubit_dependency[gate->second]);
        qubit_dependency[gate->second] = inst_idx;
      }
      nlohmann::json dependency = {
          {"qubit", std::vector<std::size_t>{}},
          {"rydberg", rydberg_dependency[rydberg_idx]}};
      rydberg_dependency[rydberg_idx] = inst_idx;
      dependency["qubit"] = set_qubit_dependency;
      write_gate_instruction(inst_idx, rydberg_idx, result_gate, dependency);
    }

    // process single-qubit gates
    const std::size_t inst_idx =
        static_cast<T*>(this)->result.instructions.size();
    std::vector<nlohmann::json> result_gate;
    std::unordered_set<std::size_t> set_qubit_dependency;
    for (const auto* gate_info : list_1q_gate) {
      // collect qubit dependency
      const auto qubit = gate_info->getTargets().front();
      set_qubit_dependency.emplace(qubit_dependency[qubit]);
      qubit_dependency[qubit] = inst_idx;
      result_gate.emplace_back(
          nlohmann::json{{"name", gate_info->getName()}, {"q", qubit}});
    }
    nlohmann::json dependency = {"qubit", set_qubit_dependency};
    if (!result_gate.empty()) {
      write_1q_gate_instruction(inst_idx, result_gate, dependency,
                                gate_mapping);
    }
  }

  auto write_gate_instruction(const std::size_t inst_idx,
                              const std::size_t rydberg_idx,
                              const std::vector<nlohmann::json>& result_gate,
                              const nlohmann::json& dependency) -> void {
    static_cast<T*>(this)->result.instructions.emplace_back(
        nlohmann::json{{"type", "rydberg"},
                       {"id", inst_idx},
                       {"zone_id", rydberg_idx},
                       {"gates", result_gate},
                       {"dependency", dependency}});
  }

  auto write_1q_gate_instruction(
      std::size_t inst_idx, const std::vector<nlohmann::json>& result_gate,
      const nlohmann::json& dependency,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& gate_mapping)
      -> void {
    std::vector<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>>
        locs{};
    for (const auto& gate : result_gate) {
      locs.emplace_back(gate["q"], std::get<0>(gate_mapping[gate["q"]])->id,
                        std::get<1>(gate_mapping[gate["q"]]),
                        std::get<2>(gate_mapping[gate["q"]]));
    }

    static_cast<T*>(this)->result.instructions.emplace(
        nlohmann::json{{"type", "1qGate"},
                       {"unitary", "u3"},
                       {"id", inst_idx},
                       {"locs", locs},
                       {"gates", result_gate},
                       {"dependency", dependency}});
  }

  /// construct reverse movement layer by processing the forward movement
  auto construct_reverse_layer(
      const std::size_t id_layer_start,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>&
          initial_mapping,
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& final_mapping)
      -> void {
    const std::size_t id_layer_end =
        static_cast<T*>(this)->result.instructions.size();
    for (std::size_t layer = id_layer_start; layer < id_layer_end; ++layer) {
      if (static_cast<T*>(this)->result.instructions[layer]["type"] ==
          "rydberg") {
        break;
      } else {
        // process a rearrangement layer
        const std::size_t inst_idx =
            static_cast<T*>(this)->result.instructions.size();
        nlohmann::json dependency = {{"qubit", std::vector<std::size_t>{}},
                                     {"site", std::vector<std::size_t>{}}};
        // process aod dependency
        std::unordered_set<std::size_t> set_qubit_dependency{};
        std::unordered_set<std::size_t> set_site_dependency{};
        nlohmann::json list_aod_qubits =
            static_cast<T*>(this)->result.instructions[layer]["aod_qubits"];
        std::vector<std::vector<
            std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>
            list_end_location{};
        std::vector<std::vector<
            std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>>
            list_begin_location{};
        for (const auto& sub_list_qubits : list_aod_qubits) {
          std::vector<
              std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>
              row_begin_location{};
          std::vector<
              std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>>
              row_end_location{};
          for (const auto& q : sub_list_qubits) {
            row_begin_location.emplace_back(q, std::get<0>(initial_mapping[q]),
                                             std::get<1>(initial_mapping[q]),
                                             std::get<2>(initial_mapping[q]));
            row_end_location.emplace_back(q, std::get<0>(final_mapping[q]),
                                           std::get<1>(final_mapping[q]),
                                           std::get<2>(final_mapping[q]));
            // process site dependency
            std::tuple site_key{std::get<0>(final_mapping[q]),
                                std::get<1>(final_mapping[q]),
                                std::get<2>(final_mapping[q])};
            if (const auto& it = site_dependency.find(site_key);
                it != site_dependency.end()) {
              set_site_dependency.emplace(it->second);
            }
            site_key = {std::get<0>(initial_mapping[q]),
                        std::get<1>(initial_mapping[q]),
                        std::get<2>(initial_mapping[q])};
            site_dependency[site_key] = inst_idx;
            // collect qubit dependency
            set_qubit_dependency.emplace(qubit_dependency[q]);
            qubit_dependency[q] = inst_idx;
          }

          list_begin_location.emplace_back(row_begin_location);
          list_end_location.emplace_back(row_end_location);
        }
        dependency["qubit"] = set_qubit_dependency;
        dependency["site"] = set_site_dependency;
        write_rearrangement_instruction(inst_idx, list_aod_qubits,
                                        list_begin_location, list_end_location,
                                        dependency);
      }
    }
  }

          /// processs the aod assignment between two Rydberg stages
        auto aod_assignment(const std::size_t id_layer_start) -> void {
          std::vector list_instruction_duration(2, std::vector<std::pair<double, std::size_t>>{});
          const std::size_t id_layer_end = static_cast<T*>(this)->result.instructions.size();
          std::size_t duration_idx = 0;
          std::vector<std::size_t> list_gate_layer_idx{};
          for (std::size_t idx = id_layer_start; idx < id_layer_end; ++idx) {
            if (static_cast<T*>(this)->result.instructions[idx]["type"] != "rearrangeJob") {
              duration_idx = 1;
              list_gate_layer_idx.emplace_back(idx);
              continue;
            }
            const double duration = get_duration(static_cast<T*>(this)->result.instructions[idx]);
            list_instruction_duration[duration_idx].emplace_back(duration, idx);
          }
          std::sort(list_instruction_duration[0].begin(), list_instruction_duration[0].end(), std::greater<std::pair<double, std::size_t>>());
          std::sort(list_instruction_duration[1].begin(), list_instruction_duration[1].end(), std::greater<std::pair<double, std::size_t>>());
          // assign instruction according to the duration in descending order
          for (const std::size_t i : {0UL, 1UL}) {
            for (const auto& item : list_instruction_duration[i]) {
              const double duration = item.first;
              const auto inst = static_cast<T*>(this)->result.instructions[item.second];
              const auto& [begin_time, aod_id] = aod_end_time.top();
              aod_end_time.pop();
              begin_time = std::max(begin_time, get_begin_time(item.second, inst["dependency"]));
              const double end_time = begin_time + duration;
              inst["dependency"]["aod"] = aod_dependency[aod_id];
              aod_dependency[aod_id] = item.second;
              inst["begin_time"] = begin_time;
              inst["end_time"] = end_time;
              inst["aod_id"] = aod_id;
              aod_end_time.emplace(end_time, aod_id);
              for (auto& detail_inst : inst["insts"]) {
                detail_inst["begin_time"] += begin_time;
                detail_inst["end_time"] += begin_time;
              }
              if (static_cast<T*>(this)->result.runtime < end_time) {
                static_cast<T*>(this)->result.runtime = end_time;
              }
            }
            if (i == 0) {
              for (const auto gate_layer_idx : list_gate_layer_idx) {
                // laser scheduling
                nlohmann::json& inst = static_cast<T*>(this)->result.instructions[gate_layer_idx];
                const auto begin_time = get_begin_time(gate_layer_idx, inst["dependency"]);
                const double end_time = inst["type"] == "rydberg" ? begin_time + static_cast<T*>(this)->architecture.time_rydberg : begin_time + (static_cast<T*>(this)->architecture.time_1qGate * inst["gates"].size()); // for sequential gate execution
                if (static_cast<T*>(this)->result.runtime < end_time) {
                  static_cast<T*>(this)->result.runtime = end_time;
                }
                inst["begin_time"] = begin_time;
                inst["end_time"] = end_time;
              }
            }
          }
  }

  auto get_begin_time(const std::size_t cur_inst_idx, nlohmann::json dependency)
      -> double {
    double begin_time = 0;
    for (const auto& dependency_type : dependency.items()) {
      if (dependency_type.value().is_number_integer()) {
        const auto inst_idx = dependency[dependency_type];
        if (begin_time <
            static_cast<T*>(this)->result.instructions[inst_idx]["end_time"]) {
          begin_time =
              static_cast<T*>(this)->result.instructions[inst_idx]["end_time"];
        }
      } else {
        if (dependency_type.key() == "site") {
          for (const std::size_t inst_idx : dependency_type.value()) {
            if (static_cast<T*>(this)->result.instructions[inst_idx]["type"] ==
                "rearrangeJob") {
              // find the time that the instruction finish atom transfer
              double atom_transfer_finish_time = 0;
              for (const auto& detail_inst :
                         static_cast<T*>(this)
                             ->result.instructions[inst_idx]["insts"]) {
                std::string inst_type = detail_inst["type"];
                inst_type = inst_type.substr(0, inst_type.find(":"));
                if (inst_type == "activate") {
                  atom_transfer_finish_time = std::max(
                      detail_inst["end_time"], atom_transfer_finish_time);
                }
              }
              // find the time until dropping of the qubits
              double atom_transfer_begin_time = 0;
              for (const auto& detail_inst :
                   static_cast<T*>(this)
                       ->result.instructions[cur_inst_idx]["insts"]) {
                std::string inst_type = detail_inst["type"];
                inst_type = inst_type.substr(0, inst_type.find(":"));
                if (inst_type == "deactivate") {
                  atom_transfer_begin_time = std::max(detail_inst["begin_time"],
                                                      atom_transfer_begin_time);
                  break;
                }
              }
              const double tmp_begin_time =
                  atom_transfer_finish_time -
                  atom_transfer_begin_time;
              if (begin_time < tmp_begin_time) {
                begin_time = tmp_begin_time;
              }
            } else {
              if (begin_time <
                  static_cast<T*>(this)
                      ->result.instructions[inst_idx]["end_time"]) {
                begin_time = static_cast<T*>(this)
                                 ->result.instructions[inst_idx]["end_time"];
              }
            }
          }
        } else {
          for (const std::size_t inst_idx : dependency_type.value()) {
            if (begin_time < static_cast<T*>(this)
                                 ->result.instructions[inst_idx]["end_time"]) {
              begin_time = static_cast<T*>(this)
                               ->result.instructions[inst_idx]["end_time"];
            }
          }
        }
      }
    }
    return begin_time;
  }

  auto get_duration(nlohmann::json& inst) -> double {
    auto& list_detail_inst = inst["insts"];
    double duration = 0;

    for (auto& detail_inst : list_detail_inst) {
      std::string inst_type = detail_inst["type"];
      inst_type = inst_type.substr(0, inst_type.find(":"));
      detail_inst["begin_time"] = duration;
      if (inst_type == "activate" || inst_type == "deactivate") {
        duration += static_cast<T*>(this)->architecture.time_atom_transfer;
        detail_inst["end_time"] = duration;
      } else if (inst_type == "move") {
        double move_duration = 0;
        for (std::size_t r = 0; r < detail_inst["row_y_begin"].size(); ++r) {
          const auto& row_begin = detail_inst["row_y_begin"][r];
          const auto& row_end = detail_inst["row_y_end"][r];
          for (std::size_t c = 0; c < detail_inst["col_x_begin"].size(); ++c) {
            const auto& col_begin = detail_inst["col_x_begin"][c];
            const auto& col_end = detail_inst["col_x_end"][c];
            const double tmp =
                static_cast<T*>(this)->architecture.movement_duration(
                    col_begin, row_begin, col_end, row_end);
            if (move_duration < tmp) {
              move_duration = tmp;
            }
          }
        }
        detail_inst["end_time"] = move_duration + duration;
        duration += move_duration;
      } else {
        throw std::invalid_argument(
            "[ERROR] Invalid instruction type in `get_duration`, must be one "
            "of 'activate', 'deactivate', 'move'");
      }
    }

    return duration;
  }

  auto expand_arrangement(
      const nlohmann::json& inst,
      const std::vector<
          std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>&
          begin_location,
      const std::vector<
          std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>&
          end_location) -> nlohmann::json {
    nlohmann::json details{};

    // ---------------------- find out number of cols ----------------------
    std::vector<std::size_t> all_col_x{}; // all the x coord of qubits
    std::vector<nlohmann::json> coords{}; // coords of qubits
    // these coords are going to be updated as we construct the detail insts

    for (const auto& locs : begin_location) {
      std::vector<nlohmann::json> coords_row{};
      for (const auto& [q, slm, r, c] : locs) {
        const auto& [x, y] =
            static_cast<T*>(this)->architecture.exact_SLM_location(slm, c, r);
        coords_row.emplace_back(nlohmann::json{{"id", q}, {"x", x}, {"y", y}});
        all_col_x.emplace_back(q);
      }
      coords.emplace_back(std::move(coords_row));
    }
    const auto init_coords = coords;

    std::sort(all_col_x.begin(), all_col_x.end());

    // assign AOD column ids based on all x coords needed
    std::unordered_map<std::size_t, std::size_t> col_x_to_id{};
    col_x_to_id.reserve(all_col_x.size());
    for (std::size_t i = 0; i < all_col_x.size(); ++i) {
      col_x_to_id.emplace(all_col_x[i], i);
    }
    // ---------------------------------------------------------------------

    // -------------------- activation and parking -------------------------
    std::vector<std::size_t>
        all_col_idx_sofar{}; // which col has been activated
    for (std::size_t row_id = 0; row_id < begin_location.size(); ++row_id) {
      // for each row
      const auto& locs = begin_location[row_id];
      const std::size_t row_y =
          static_cast<T*>(this)
              ->architecture
              .exact_SLM_location(std::get<1>(locs.front()),
                                  std::get<2>(locs.front()),
                                  std::get<3>(locs.front()))
              .second;
      const std::pair row_loc{std::get<1>(locs.front())->id,
                              std::get<2>(locs.front())};
      // before activation, adjust column position. This is necessary
      // whenever cols are parked (the `parking` movement below).
      nlohmann::json shift_back = {
          {"type", "move"},
          {"move_type", "before"},
          {"row_id", std::vector<std::size_t>{}},
          {"row_y_begin", std::vector<std::size_t>{}},
          {"row_y_end", std::vector<std::size_t>{}},
          {"row_loc_begin", std::vector<std::size_t>{}},
          {"row_loc_end", std::vector<std::size_t>{}},
          {"col_id", std::vector<std::size_t>{}},
          {"col_x_begin", std::vector<std::size_t>{}},
          {"col_x_end", std::vector<std::size_t>{}},
          {"col_loc_begin",
           std::vector<std::pair<std::int64_t, std::int64_t>>{}},
          {"col_loc_end", std::vector<std::pair<std::int64_t, std::int64_t>>{}},
          {"begin_coord", coords},
          {"end_coord", std::vector<std::size_t>{}}};

      // activate one row and some columns
      nlohmann::json activate = {
          {"type", "activate"},
          {"row_id", std::vector{row_id}},
          {"row_y", std::vector{row_y}},
          {"row_loc", std::vector{row_loc}},
          {"col_id", std::vector<std::size_t>{}},
          {"col_x", std::vector<std::size_t>{}},
          {"col_loc", std::vector<std::pair<std::size_t, std::size_t>>{}},
      };

      for (std::size_t j = 0; j < locs.size(); ++j) {
        const auto& [q, slm, r, c] = locs[j];
        const std::size_t col_x =
            static_cast<T*>(this)
                ->architecture.exact_SLM_location(slm, r, c)
                .first;
        const std::pair col_loc{slm->id, c};
        const auto col_id = col_x_to_id[col_x];
        if (std::find(all_col_idx_sofar.cbegin(), all_col_idx_sofar.cend(),
                      col_id) == all_col_idx_sofar.end()) {
          // the col hasn't been activated, so there's no shift back
          // and we need to activate it at `col_x`.`
          all_col_idx_sofar.emplace_back(col_id);
          activate["col_id"].emplace_back(col_id);
          activate["col_x"].emplace_back(col_x);
          activate["col_loc"].emplace_back(col_loc);
        } else {
          // the col has been activated, thus parked previously and we
          // need the shift back, but we do not activate again.
          shift_back["col_id"].emplace_back(col_id);
          shift_back["col_x_begin"].emplace_back(col_x + PARKING_DIST);
          shift_back["col_x_end"].emplace_back(col_x);
          shift_back["col_loc_begin"].emplace_back(std::pair{-1, -1});
          shift_back["col_loc_end"].emplace_back(col_loc);
          // since there's a shift, update the coords of the qubit
          coords[row_id][j]["x"] = col_x;
        }
      }

      shift_back["end_coord"] = coords;

      if (!shift_back["col_id"].empty()) {
        details.emplace_back(shift_back);
      }
      details.emplace_back(activate);

      if (row_id < inst["begin_locs"].size() - 1) {
        // parking movement after the activation
        // parking is required if we have activated some col, and there is
        // some qubit we don't want to pick up at the intersection of this
        // col and some future row to activate. We just always park here.
        // the last parking is not needed since there's a big move after it.
        nlohmann::json parking = {
            {"type", "move"},
            {"move_type", "after"},
            {"row_id", std::vector{row_id}},
            {"row_y_begin", std::vector{row_y}},
            {"row_y_end", std::vector{row_y + PARKING_DIST}},
            {"row_loc_begin", std::vector{row_loc}},
            {"row_loc_end", std::vector{std::pair{-1, -1}}},
            {"col_id", std::vector<std::size_t>{}},
            {"col_x_begin", std::vector<std::size_t>{}},
            {"col_x_end", std::vector<std::size_t>{}},
            {"col_loc_begin", std::vector<std::size_t>{}},
            {"col_loc_end", std::vector<std::size_t>{}},
            {"begin_coord", coords},
            {"end_coord", std::vector<std::size_t>{}}};
        for (std::size_t j = 0; j < locs.size(); ++j) {
          const auto& [q, slm, r, c] = locs[j];
          const std::size_t col_x =
              static_cast<T*>(this)
                  ->architecture.exact_SLM_location(slm, r, c)
                  .first;
          const std::pair col_loc{slm, c};
          const auto col_id = col_x_to_id[col_x];
          // all columns used in this row are parked after the activation
          parking["col_id"].emplace_back(col_id);
          parking["col_x_begin"].emplace_back(col_x);
          parking["col_x_end"].emplace_back(col_x + PARKING_DIST);
          parking["col_loc_begin"].emplace_back(col_loc);
          parking["col_loc_end"].emplace_back(std::pair{-1, -1});
          coords[row_id][j]["x"] = parking["col_x_end"].back();
          coords[row_id][j]["y"] = parking["row_y_end"].front();
        }
        parking["end_coord"] = coords;
        details.emplace_back(parking);
      }
    }
    // ---------------------------------------------------------------------

    // ------------------------- big move ----------------------------------
    nlohmann::json big_move = {
        {"type", "move:big"},
        {"move_type", "big"},
        {"row_id", std::vector<std::size_t>{}},
        {"row_y_begin", std::vector<std::size_t>{}},
        {"row_y_end", std::vector<std::size_t>{}},
        {"row_loc_begin", std::vector<std::size_t>{}},
        {"row_loc_end", std::vector<std::size_t>{}},
        {"col_id", std::vector<std::size_t>{}},
        {"col_x_begin", std::vector<std::size_t>{}},
        {"col_x_end", std::vector<std::size_t>{}},
        {"col_loc_begin", std::vector<std::size_t>{}},
        {"col_loc_end", std::vector<std::size_t>{}},
        {"begin_coord", coords},
        {"end_coord", std::vector<std::size_t>{}},
    };

    for (std::size_t row_id = 0; row_id < begin_location.size(); ++row_id) {
      const auto& begin_locs = begin_location[row_id];
      const auto& end_locs = end_location[row_id];
      big_move["row_id"].emplace_back(row_id);
      big_move["row_y_begin"].emplace_back(coords[row_id][0]["y"]);
      if (init_coords[row_id][0]["y"] == coords[row_id][0]["y"]) {
        // AOD row is align with SLM row
        big_move["row_loc_begin"].emplace_back(std::pair{
            std::get<1>(begin_locs[0])->id, std::get<2>(begin_locs[0])});
      } else {
        big_move["row_loc_begin"].emplace_back(std::pair{-1, -1});
      }
      big_move["row_y_end"].emplace_back(
          static_cast<T*>(this)
              ->architecture
              .exact_SLM_location(std::get<1>(end_locs[0]),
                                  std::get<2>(end_locs[0]),
                                  std::get<3>(end_locs[0]))
              .second);
      big_move["row_loc_end"].emplace_back(
          std::pair{std::get<1>(end_locs[0]), std::get<2>(end_locs[0])});

      for (std::size_t j = 0; j < begin_locs.size(); ++j) {
        const auto& begin_loc = begin_locs[j];
        const auto& end_loc = end_locs[j];
        const auto col_x = static_cast<T*>(this)
                               ->architecture
                               .exact_SLM_location(std::get<1>(begin_loc),
                                                   std::get<2>(begin_loc),
                                                   std::get<3>(begin_loc))
                               .first;
        const auto col_id = col_x_to_id[col_x];

        if (!big_move["col_id"].contains(col_id)) {
          // the movement of this rol has not been recorded before
          big_move["col_id"].emplace_back(col_id);
          big_move["col_x_begin"].emplace_back(coords[row_id][j]["x"]);
          if (init_coords[row_id][j]["x"] == coords[row_id][j]["x"]) {
            // AOD col is align with SLM col
            big_move["col_loc_begin"].emplace_back(
                std::pair{std::get<1>(begin_loc)->id, std::get<3>(begin_loc)});
          } else {
            big_move["col_loc_begin"].emplace_back(std::pair{-1, -1});
          }
          big_move["col_x_end"].emplace_back(
              static_cast<T*>(this)
                  ->architecture
                  .exact_SLM_location(std::get<1>(end_loc),
                                      std::get<2>(end_loc),
                                      std::get<3>(end_loc))
                  .first);
          big_move["col_loc_end"].emplace_back(
              std::pair{std::get<1>(end_loc)->id, std::get<3>(end_loc)});
        }
        // whether or not the movement of this col has been considered
        // before, we need to update the coords of the qubit.
        coords[row_id][j]["x"] =
            static_cast<T*>(this)
                ->architecture
                .exact_SLM_location(std::get<1>(end_loc), std::get<2>(end_loc),
                                    std::get<3>(end_loc))
                .first;
        coords[row_id][j]["y"] =
            static_cast<T*>(this)
                ->architecture
                .exact_SLM_location(std::get<1>(end_locs[0]),
                                    std::get<2>(end_locs[0]),
                                    std::get<3>(end_locs[0]))
                .second;
      }
    }
    big_move["end_coord"] = coords;
    details.emplace_back(big_move);
    // ---------------------------------------------------------------------

    // --------------------------- deactivation ----------------------------
    details.emplace_back(
        nlohmann::json{{"type", "deactivate"},
                       {"row_id", std::vector<std::size_t>{}},
                       {"col_id", std::vector<std::size_t>{}}});
    for (std::size_t i = 0; i < begin_location.size(); ++i) {
      details.back()["row_id"].emplace_back(i);
    }
    for (std::size_t i = 0; i < all_col_x.size(); ++i) {
      details.back()["col_id"].emplace_back(i);
    }
    // ---------------------------------------------------------------------

    for (std::size_t inst_counter = 0; inst_counter < details.size();
         ++inst_counter) {
      auto& detail_inst = details[inst_counter];
      detail_inst["id"] = inst_counter;
    }
    return details;
  }
};

} // namespace na
