#pragma once

#include "na/NAComputation.hpp"
#include "na/azac/CompilerBase.hpp"
#include "na/operations/NALocalOperation.hpp"
#include "na/operations/NAOperation.hpp"
#include "na/operations/NAShuttlingOperation.hpp"

#include <chrono>
#include <type_traits>

namespace na {

template <typename T> class Router {
  static_assert(std::is_base_of_v<CompilerBase<T>, T>,
                "T must be a subclass of CompilerBase");

  /// constant, the distance of AOD row and col to some trap. We use 1um here.
  constexpr std::size_t PARKING_DIST = 1;
  std::vector<std::shared_ptr<Point>> current_qubit_pos;

  /// generate rearrangement layers between two Rydberg layers
  auto route_qubit() -> void {
    std::vector<std::pair<std::size_t, std::size_t>> aod_end_time{};
    aod_end_time.reserve(static_cast<T*>(this)->architecture.dict_AOD.size());
    for (std::size_t i = 0;
         i < static_cast<T*>(this)->architecture.dict_AOD.size(); ++i) {
      aod_end_time.emplace_back(0, i);
    }
    std::vector<std::size_t> aod_dependency(
        static_cast<T*>(this)->architecture.dict_AOD.size(), 0);
    std::vector<std::size_t> rydberg_dependency(
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
      ++batch;
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
        std:
          sort(remain_graph.begin(), remain_graph.end(),
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
    current_qubit_pos.clear();
    for (std::size_t i = 0; i < static_cast<T*>(this)->n_q; ++i) {
      const auto& pos = current_qubit_pos.emplace_back(std::apply(
              std::make_shared<Point>,
              std::apply(static_cast<T*>(this)->architecture.exact_SLM_location,
                         static_cast<T*>(this)->qubit_mapping[0][i])));
      static_cast<T*>(this)->result.instructions.emplaceInitialPosition(std::move(pos));
    }

    // process single-qubit gates
    const std::vector<qc::StandardOperation>& list_1q_gate =
        static_cast<T*>(this)->dict_g_1q_parent.at(nullptr);
    std::vector<NAOperation> result_gate{};
    for (const auto& gate_info : list_1q_gate) {
      // collect qubit dependency
      const auto q = gate_info.getTargets().front();
      result_gate.emplace_back(NALocalOperation{
          {gate_info.getType(), 0},
          gate_info.getParameter(),
          current_qubit_pos[q]});
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
      std::vector<const std::vector<std::size_t>&> list_aod_qubits{};
      std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>> list_end_location{};
      std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>> list_begin_location{};
      // process aod dependency
      const std::size_t inst_idx = static_cast<T*>(this)->result.instructions.size();

      for (const auto& [dict_key, dict_value] : pickup_dict) {
        // collect set of aod qubits to pick up
        list_aod_qubits.emplace_back(dict_value);
        std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>> row_begin_location{};
        std::vector<std::tuple<std::size_t, const SLM*, std::size_t, std::size_t>> row_end_location{};
        for (const auto q : pickup_dict[dict_key]) {
          // collect qubit begin location
          row_begin_location.emplace_back(q, std::get<0>(initial_mapping[q]), std::get<1>(initial_mapping[q]), std::get<2>(initial_mapping[q]));

          // collect qubit end location
          row_end_location.emplace_back(q, std::get<0>(final_mapping[q]), std::get<1>(final_mapping[q]), std::get<2>(final_mapping[q]));
        }
        list_begin_location.emplace_back(std::move(row_begin_location));
        list_end_location.emplace_back(std::move(row_end_location));
      }
      write_rearrangement_instruction(inst_idx, list_aod_qubits, list_begin_location, list_end_location);
    }

    auto write_rearrangement_instruction(const std::size_t inst_idx, const std::vector<const std::vector<size_t>&>& aod_qubits,
      const std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>& begin_location,
      const std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>& end_location) -> void {
      expand_arrangement(inst_idx, aod_qubits, begin_location, end_location);
    }
/*
        /// generate a layer for gate execution
    auto process_gate_layer(layer: int, gate_mapping: list) -> void {
      list_gate_idx = self.gate_scheduling_idx[layer]
      list_gate = self.gate_scheduling[layer]
      list_1q_gate = self.gate_1q_scheduling[layer]
      dict_gate_zone = dict()
      for i in range(len(list_gate)):
          slm_idx = gate_mapping[list_gate[i][0]][0]
          zone_idx = static_cast<T*>(this)->architecture.dict_SLM[slm_idx].entanglement_id
          if zone_idx not in dict_gate_zone:
              dict_gate_zone[zone_idx] = [i]
          else:
              dict_gate_zone[zone_idx].append(i)
      for rydberg_idx in dict_gate_zone:
          result_gate = [{"id": list_gate_idx[i], "q0": list_gate[i][0], "q1": list_gate[i][1]} for i in dict_gate_zone[rydberg_idx]]
          set_qubit_dependency = set()
          inst_idx = len(static_cast<T*>(this)->result.instructions)
          for gate_idx in dict_gate_zone[rydberg_idx]:
              gate = list_gate[gate_idx]
              // collect qubit dependency
              set_qubit_dependency.add(self.qubit_dependency[gate[0]])
              self.qubit_dependency[gate[0]] = inst_idx
              set_qubit_dependency.add(self.qubit_dependency[gate[1]])
              self.qubit_dependency[gate[1]] = inst_idx
          dependency = { "qubit": [], "rydberg": self.rydberg_dependency[rydberg_idx]}
      self.rydberg_dependency[rydberg_idx] = inst_idx
      dependency["qubit"] = list(set_qubit_dependency)
      self.write_gate_instruction(inst_idx, rydberg_idx, result_gate, dependency)

  // process single-qubit gates
  inst_idx = len(static_cast<T*>(this)->result.instructions)
  result_gate = []
  set_qubit_dependency = set()
  for gate_info in list_1q_gate:
      // collect qubit dependency
      set_qubit_dependency.add(self.qubit_dependency[gate_info[1]])
      self.qubit_dependency[gate_info[1]] = inst_idx
      result_gate.append({
          "name": gate_info[0],
          "q": gate_info[1]
      })
  dependency = { "qubit": []}
      dependency["qubit"] = list(set_qubit_dependency)
      if len(result_gate) > 0:
          self.write_1q_gate_instruction(inst_idx, result_gate, dependency, gate_mapping)
    }

    auto write_gate_instruction(inst_idx: int, rydberg_idx: int, result_gate: list, dependency: dict) -> void {
      static_cast<T*>(this)->result.instructions.append(
          {
              "type": "rydberg",
              "id": inst_idx,
              "zone_id": rydberg_idx,
              "gates": result_gate,
              "dependency": dependency
          }
      )
    }

    auto write_1q_gate_instruction(inst_idx: int, result_gate: list, dependency: dict, gate_mapping: list) -> void {
      locs = []
      for gate in result_gate:
          locs.append((gate["q"], gate_mapping[gate["q"]][0], gate_mapping[gate["q"]][1], gate_mapping[gate["q"]][2]))

      static_cast<T*>(this)->result.instructions.append(
          {
              "type": "1qGate",
              "unitary": "u3",
              "id": inst_idx,
              "locs": locs,
              "gates": result_gate,
              "dependency": dependency
          }
      )
    }

    auto construct_reverse_layer(id_layer_start: int, initial_mapping: list, final_mapping: list) ->void {
      """
      construct reverse movement layer by processing the forward movement
      """
      id_layer_end = len(static_cast<T*>(this)->result.instructions)
      for layer in range(id_layer_start, id_layer_end):
          if static_cast<T*>(this)->result.instructions[layer]["type"] == "rydberg":
              break
          else:
              // process a rearrangement layer
              inst_idx = len(static_cast<T*>(this)->result.instructions)
              dependency = {
        "qubit": [],
        "site": [],
    }
// process aod dependency
      set_qubit_dependency = set()
      set_site_dependency = set()
      list_aod_qubits = static_cast<T*>(this)->result.instructions[layer]["aod_qubits"]
      list_end_location = []
      list_begin_location = []
      for sub_list_qubits in list_aod_qubits:
          row_begin_location = []
          row_end_location = []
          for q in sub_list_qubits:
              row_begin_location.append([q, initial_mapping[q][0], initial_mapping[q][1], initial_mapping[q][2]])
              row_end_location.append([q, final_mapping[q][0], final_mapping[q][1], final_mapping[q][2]])
              // process site dependency
              site_key = (final_mapping[q][0], final_mapping[q][1], final_mapping[q][2])
              if site_key in self.site_dependency:
                  set_site_dependency.add(self.site_dependency[site_key])
              site_key = (initial_mapping[q][0], initial_mapping[q][1], initial_mapping[q][2])
              self.site_dependency[site_key] = inst_idx
              // collect qubit dependency
              set_qubit_dependency.add(self.qubit_dependency[q])
              self.qubit_dependency[q] = inst_idx

          list_begin_location.append(row_begin_location)
          list_end_location.append(row_end_location)
      dependency["qubit"] = list(set_qubit_dependency)
      dependency["site"] = list(set_site_dependency)
      self.write_rearrangement_instruction(inst_idx,
                                           list_aod_qubits,
                                           list_begin_location,
                                           list_end_location,
                                           dependency)
    }

      /// processs the aod assignment between two Rydberg stages
    auto aod_assignment(self, id_layer_start: int) -> void {
      list_instruction_duration = [[], []]
      id_layer_end = len(static_cast<T*>(this)->result.instructions)
      duration_idx = 0
      list_gate_layer_idx = []
      for idx in range(id_layer_start, id_layer_end):
          if static_cast<T*>(this)->result.instructions[idx]["type"] != "rearrangeJob":
              duration_idx = 1
              list_gate_layer_idx.append(idx)
              continue
          duration = self.get_duration(static_cast<T*>(this)->result.instructions[idx])
          list_instruction_duration[duration_idx].append((duration, idx))
      list_instruction_duration[0] = sorted(list_instruction_duration[0], reverse=True)
      list_instruction_duration[1] = sorted(list_instruction_duration[1], reverse=True)
      // assign instruction according to the duration in descending order
      for i in range(2):
          for item in list_instruction_duration[i]:
              duration = item[0]
              inst = static_cast<T*>(this)->result.instructions[item[1]]
              begin_time, aod_id = heapq.heappop(self.aod_end_time)
              begin_time = max(begin_time, self.get_begin_time(item[1], inst["dependency"]))
              end_time = begin_time + duration
              inst["dependency"]["aod"] = self.aod_dependency[aod_id]
              self.aod_dependency[aod_id] = item[1]
              inst["begin_time"] = begin_time
              inst["end_time"] = end_time
              inst["aod_id"] = aod_id
              heapq.heappush(self.aod_end_time, (end_time, aod_id))
              for detail_inst in inst["insts"]:
                  detail_inst["begin_time"] += begin_time
                  detail_inst["end_time"] += begin_time
              if self.result_json["runtime"] < end_time:
                  self.result_json["runtime"] = end_time
          if i == 0:
              for gate_layer_idx in list_gate_layer_idx:
                  // laser scheduling
                  inst = static_cast<T*>(this)->result.instructions[gate_layer_idx]
                  begin_time = self.get_begin_time(gate_layer_idx, inst["dependency"])
                  if inst["type"] == "rydberg":
                      end_time = begin_time + static_cast<T*>(this)->architecture.time_rydberg
                  else:
                      end_time = begin_time + (static_cast<T*>(this)->architecture.time_1qGate * len(inst["gates"])) + self.common_1q # for sequential gate execution
                  if self.result_json["runtime"] < end_time:
                      self.result_json["runtime"] = end_time
                  inst["begin_time"] = begin_time
                  inst["end_time"] = end_time
  }

    auto get_begin_time(cur_inst_idx: int, dependency: dict) -> void {
      begin_time = 0
      for dependency_type in dependency:
          if isinstance(dependency[dependency_type], int):
              inst_idx = dependency[dependency_type]
              if begin_time < static_cast<T*>(this)->result.instructions[inst_idx]["end_time"]:
                  begin_time = static_cast<T*>(this)->result.instructions[inst_idx]["end_time"]
          else:
              if dependency_type == "site":
                  for inst_idx in dependency[dependency_type]:
                      if static_cast<T*>(this)->result.instructions[inst_idx]["type"] == "rearrangeJob":
                          # find the time that the instruction finish atom transfer
                          # !
                          # atom_transfer_finish_time = static_cast<T*>(this)->result.instructions[inst_idx]["begin_time"] + 15
                          atom_transfer_finish_time = 0
                          for detail_inst in static_cast<T*>(this)->result.instructions[inst_idx]["insts"]:
                              inst_type = detail_inst["type"].split(":")[0]
                              if inst_type == "activate":
                                  atom_transfer_finish_time = max(detail_inst["end_time"], atom_transfer_finish_time)
                          # find the time until dropping of the qubits
                          # atom_transfer_begin_time = 15
                          atom_transfer_begin_time = 0
                          for detail_inst in static_cast<T*>(this)->result.instructions[cur_inst_idx]["insts"]:
                              # print(detail_inst["type"].split(":"))
                              inst_type = detail_inst["type"].split(":")[0]
                              if inst_type == "deactivate":
                                  atom_transfer_begin_time = max(detail_inst["begin_time"], atom_transfer_begin_time)
                                  break
                          tmp_begin_time = atom_transfer_finish_time - atom_transfer_begin_time
                          if begin_time < tmp_begin_time:
                              begin_time = tmp_begin_time
                          # print("cur_inst_idx: ", cur_inst_idx)
                          # print("atom_transfer_finish_time: ", atom_transfer_finish_time)
                          # print("atom_transfer_begin_time: ", atom_transfer_begin_time)
                          # print("begin time for site depend: ", tmp_begin_time)
                      else:
                          if begin_time < static_cast<T*>(this)->result.instructions[inst_idx]["end_time"]:
                              begin_time = static_cast<T*>(this)->result.instructions[inst_idx]["end_time"]
              else:
                  for inst_idx in dependency[dependency_type]:
                      if begin_time < static_cast<T*>(this)->result.instructions[inst_idx]["end_time"]:
                          begin_time = static_cast<T*>(this)->result.instructions[inst_idx]["end_time"]
      return begin_time
    }

    auto get_duration(self, inst: dict) -> void {
      list_detail_inst = inst["insts"]
      duration = 0

      for detail_inst in list_detail_inst:
          inst_type = detail_inst["type"].split(":")[0]
          detail_inst["begin_time"] = duration
          if inst_type == "activate" or inst_type == "deactivate":
              duration += static_cast<T*>(this)->architecture.time_atom_transfer
              detail_inst["end_time"] = duration
          elif inst_type == "move":
              move_duration = 0
              for row_begin, row_end in zip(detail_inst["row_y_begin"], detail_inst["row_y_end"]):
                  for col_begin, col_end in zip(detail_inst["col_x_begin"], detail_inst["col_x_end"]):
                      tmp = static_cast<T*>(this)->architecture.movement_duration(col_begin, row_begin, col_end, row_end)
                      if move_duration < tmp:
                          move_duration = tmp
              detail_inst["end_time"] = move_duration + duration
              duration += move_duration
          else:
              raise ValueError


      return duration
    }
*/
  auto expand_arrangement(const std::size_t inst_idx, const std::vector<const std::vector<size_t>&>& aod_qubits,
    const std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>& begin_location,
    const std::vector<std::vector<std::tuple<size_t, const SLM*, size_t, size_t>>>& end_location) -> void {
    // ---------------------- find out number of cols ----------------------
    std::vector<std::size_t> all_col_x{}; // all the x coord of qubits
    std::vector<std::vector<std::tuple<size_t, size_t, size_t>>> coords{}; // coords of qubits
    // these coords are going to be updated as we construct the detail insts

    for (const auto& locs : begin_location) {
      std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> coords_row{};
      for (const auto& [q, slm, r, c] : locs) {
        const auto& [x, y] = static_cast<T*>(this)->architecture.exact_SLM_location(slm, c, r);
        coords_row.emplace_back(q, x, y);
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
    std::vector<std::size_t> all_col_idx_sofar{}; // which col has been activated
    for (std::size_t row_id = 0; row_id < aod_qubits.size(); ++row_id) {
      // for each row
      const auto& locs = begin_location[row_id];
      const std::size_t row_y = static_cast<T*>(this)->architecture.exact_SLM_location(
          std::get<1>(locs.front()),
          std::get<2>(locs.front()),
          std::get<3>(locs.front())
      ).second;
      const std::pair row_loc{std::get<1>(locs.front()), std::get<2>(locs.front())};
      /*
                // before activation, adjust column position. This is necessary
                // whenever cols are parked (the `parking` movement below).
                shift_back = {
              "type": "move",
              "move_type": "before",
              "row_id": [],
              "row_y_begin": [],
              "row_y_end": [],
              "row_loc_begin": [],
              "row_loc_end": [],
              "col_id": [],
              "col_x_begin": [],
              "col_x_end": [],
              "col_loc_begin": [],
              "col_loc_end": [],
              "begin_coord": deepcopy(coords),
              "end_coord": [],
          }

            // activate one row and some columns
            activate = {
              "type": "activate",
              "row_id": [row_id, ],
              "row_y": [row_y, ],
              "row_loc": [row_loc, ],
              "col_id": [],
              "col_x": [],
              "col_loc": [],
          }

            for j, loc in enumerate(locs):
                col_x = static_cast<T*>(this)->architecture.exact_SLM_location(
                    loc[1],
                    loc[2],
                    loc[3],
                )[0]
                col_loc = [loc[1], loc[3]]
                col_id = col_x_to_id[col_x]
                if col_id not in all_col_idx_sofar:
                    // the col hasn't been activated, so there's no shift back
                    // and we need to activate it at `col_x`.`
                    all_col_idx_sofar.append(col_id)
                    activate["col_id"].append(col_id)
                    activate["col_x"].append(col_x)
                    activate["col_loc"].append(col_loc)
                else:
                    // the col has been activated, thus parked previously and we
                    // need the shift back, but we do not activate again.
                    shift_back["col_id"].append(col_id)
                    shift_back["col_x_begin"].append(col_x + self.PARKING_DIST)
                    shift_back["col_x_end"].append(col_x)
                    shift_back["col_loc_begin"].append([-1, -1])
                    shift_back["col_loc_end"].append(col_loc)
                    // since there's a shift, update the coords of the qubit
                    coords[row_id][j]["x"] = col_x

            shift_back["end_coord"] = deepcopy(coords)

            if len(shift_back["col_id"]) != 0:
                details.append(shift_back)
            details.append(activate)

            if row_id < len(inst["begin_locs"]) - 1:
            // parking movement after the activation
            // parking is required if we have activated some col, and there is
            // some qubit we don't want to pick up at the intersection of this
            // col and some future row to activate. We just always park here.
            // the last parking is not needed since there's a big move after it.
                parking = {
              "type": "move",
              "move_type": "after",
              "row_id": [row_id, ],
              "row_y_begin": [row_y, ],
              "row_y_end": [row_y + self.PARKING_DIST],
              "row_loc_begin": [row_loc],
              "row_loc_end": [[-1, -1]],
              "col_id": [],
              "col_x_begin": [],
              "col_x_end": [],
              "col_loc_begin": [],
              "col_loc_end": [],
              "begin_coord": deepcopy(coords),
              "end_coord": [],
          }
            for j, loc in enumerate(locs):
                col_x = static_cast<T*>(this)->architecture.exact_SLM_location(
                    loc[1],
                    loc[2],
                    loc[3],
                )[0]
                col_loc = [loc[1], loc[3]]
                col_id = col_x_to_id[col_x]
                // all columns used in this row are parked after the activation
                parking["col_id"].append(col_id)
                parking["col_x_begin"].append(col_x)
                parking["col_x_end"].append(col_x + self.PARKING_DIST)
                parking["col_loc_begin"].append(col_loc)
                parking["col_loc_end"].append([-1, -1])
                coords[row_id][j]["x"] = parking["col_x_end"][-1]
                coords[row_id][j]["y"] = parking["row_y_end"][0]
            parking["end_coord"] = deepcopy(coords)
            details.append(parking)
      // ---------------------------------------------------------------------

      // ------------------------- big move ----------------------------------
      big_move = {
              "type": "move:big",
              "move_type": "big",
              "row_id": [],
              "row_y_begin": [],
              "row_y_end": [],
              "row_loc_begin": [],
              "row_loc_end": [],
              "col_id": [],
              "col_x_begin": [],
              "col_x_end": [],
              "col_loc_begin": [],
              "col_loc_end": [],
              "begin_coord": deepcopy(coords),
              "end_coord": [],
          }

            for row_id, (begin_locs, end_locs) in enumerate(zip(
                inst["begin_locs"], inst["end_locs"],
                )):

                big_move["row_id"].append(row_id)
                big_move["row_y_begin"].append(
                    coords[row_id][0]["y"]
                )
                if init_coords[row_id][0]["y"] == coords[row_id][0]["y"]:
                    // AOD row is align with SLM row
                    big_move["row_loc_begin"].append([begin_locs[0][1], begin_locs[0][2]])
                else:
                    big_move["row_loc_begin"].append([-1, -1])

                big_move["row_y_end"].append(
                    static_cast<T*>(this)->architecture.exact_SLM_location(
                        end_locs[0][1],
                        end_locs[0][2],
                        end_locs[0][3],
                    )[1]
                )
                big_move["row_loc_end"].append([end_locs[0][1], end_locs[0][2]])

                for j, (begin_loc, end_loc) in enumerate(zip(begin_locs, end_locs)):
                    col_x = static_cast<T*>(this)->architecture.exact_SLM_location(
                                begin_loc[1],
                                begin_loc[2],
                                begin_loc[3],
                            )[0]
                    col_id = col_x_to_id[col_x]

                    if col_id not in big_move["col_id"]:
                        // the movement of this rol has not been recorded before
                        big_move["col_id"].append(col_id)
                        big_move["col_x_begin"].append(coords[row_id][j]["x"])
                        if init_coords[row_id][j]["x"] == coords[row_id][j]["x"]:
                            // AOD col is align with SLM col
                            big_move["col_loc_begin"].append([begin_loc[1], begin_loc[3]])
                        else:
                            big_move["col_loc_begin"].append([-1, -1])
                        big_move["col_x_end"].append(
                            static_cast<T*>(this)->architecture.exact_SLM_location(
                                end_loc[1],
                                end_loc[2],
                                end_loc[3],
                            )[0]
                        )
                        big_move["col_loc_end"].append([end_loc[1], end_loc[3]])

                    // whether or not the movement of this col has been considered
                    // before, we need to update the coords of the qubit.
                    coords[row_id][j]["x"] = static_cast<T*>(this)->architecture.exact_SLM_location(
                                                end_loc[1],
                                                end_loc[2],
                                                end_loc[3],
                                            )[0]
                    coords[row_id][j]["y"] = static_cast<T*>(this)->architecture.exact_SLM_location(
                        end_locs[0][1],
                        end_locs[0][2],
                        end_locs[0][3],
                    )[1]

            big_move["end_coord"] = deepcopy(coords)
            details.append(big_move)
            // ---------------------------------------------------------------------

            // --------------------------- deactivation ----------------------------
            details.append({
                    "type": "deactivate",
                    "row_id": [i for i in range(len(inst["begin_locs"]))],
                    "col_id": [i for i in range(len(all_col_x))],
                })
            // ---------------------------------------------------------------------

            for inst_counter, detail_inst in enumerate(details):
                detail_inst["id"] = inst_counter
            }
            */
    }
  }
};

} // namespace na
