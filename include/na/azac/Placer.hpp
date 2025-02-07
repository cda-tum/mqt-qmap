#pragma once

#include "na/azac/Architecture.hpp"
#include "na/azac/CompilerBase.hpp"
#include "na/azac/Utils.hpp"

#include <algorithm>
#include <cstdint>
#include <tuple>
#include <type_traits>

namespace na {

/// class to find a qubit layout
template <typename T> class Placer {
protected:
  /// generate qubit initial layout
  auto place_qubit_initial() -> void {
    const auto t_p = std::chrono::system_clock::now();
    if (static_cast<T*>(this)->get_given_initial_mapping()) {
      static_cast<T*>(this)->get_qubit_mapping().emplace_back(
          *static_cast<T*>(this)->get_given_initial_mapping());
    } else {
      if (static_cast<T*>(this)->is_trivial_placement()) {
        place_trivial();
      } else {
        throw std::invalid_argument(
            "Initial placement via simulated annealing is not implemented");
      }
    }
    static_cast<T*>(this)->get_runtime_analysis().initial_placement =
        std::chrono::system_clock::now() - t_p;
  }

  /// generate qubit initial layout
  auto place_qubit_intermediate() -> void {
    const auto t_p = std::chrono::system_clock::now();
    VertexMatchingPlacer intermediatePlacer(static_cast<T*>(this)->get_architecture(),
        static_cast<T*>(this)->get_qubit_mapping().front());
    intermediatePlacer.run(static_cast<T*>(this)->get_qubit_mapping(),
                           static_cast<T*>(this)->get_gate_scheduling(),
                           static_cast<T*>(this)->is_dynamic_placement(),
                           static_cast<T*>(this)->get_reuse_qubits());
    static_cast<T*>(this)->get_qubit_mapping() = intermediatePlacer.get_mapping();

    static_cast<T*>(this)->get_runtime_analysis().intermediate_placement =
        std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now() - t_p);
  }

private:
  auto place_trivial() -> void {
    std::vector<std::unique_ptr<SLM>>::const_iterator slm_it =
        static_cast<T*>(this)->get_architecture().storage_zone.cbegin();
    const SLM* slm = slm_it->get();
    std::size_t c = 0;
    // decide whether to begin with row 0 or row n
    const double dis1 =
        static_cast<T*>(this)->get_architecture().nearest_entanglement_site_distance(
            slm, 0, c);
    const double dis2 =
        static_cast<T*>(this)->get_architecture().nearest_entanglement_site_distance(
            slm, slm->n_r - 1, c);
    std::int64_t r = dis1 < dis2 ? 0 : static_cast<std::int64_t>(slm->n_r) - 1;
    const std::int64_t step = dis1 < dis2 ? 1 : -1;
    std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>
        listPossiblePosition{};
    for (std::size_t i = 0; i < static_cast<T*>(this)->get_n_q(); ++i) {
      listPossiblePosition.emplace_back(slm, r, c);
      ++c;
      if (c % slm->n_c == 0) {
        r += step;
        c = 0;
        if (r == static_cast<std::int64_t>(slm->n_r)) {
          ++slm_it;
          slm = slm_it->get();
          if (step > 0) {
            r = static_cast<std::int64_t>(slm->n_r) - 1;
          } else {
            r = 0;
          }
        }
      }
    }
    static_cast<T*>(this)->get_qubit_mapping().emplace_back(listPossiblePosition);
  }

  /// class to find a qubit placement via vertex matching
  class VertexMatchingPlacer {
  private:
    const Architecture& architecture;
    std::vector<std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>
        mapping{};
    bool l2 = false;
    double costAtomTransfer = 0.9999;
    std::size_t nQubit = 0;
    std::vector<std::unordered_set<std::size_t>> listReuseQubits{};

  public:
    VertexMatchingPlacer(const Architecture& architecture,
        const std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>&
            initial_mapping,
        const bool l2 = false)
        : architecture(architecture), l2(l2) {
      mapping.emplace_back(initial_mapping);
      nQubit = initial_mapping.size();
    }

    [[nodiscard]] auto get_mapping() const -> const std::vector<
        std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>& {
      return mapping;
    }

    auto run(const std::vector<
                 std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>>&
                 qubit_mapping,
             const std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>&
                 list_gate,
             const bool dynamic_placement,
             const std::vector<std::unordered_set<std::size_t>>& reuse_qubits) -> void {
      this->listReuseQubits = reuse_qubits;
      std::cout << "[INFO] ZAC: Minimum-weight-full-matching-based "
                   "intermediate placement: Start\n";
      place_gate(qubit_mapping, list_gate, 0, false);
      for (std::size_t layer = 0; layer < list_gate.size(); ++layer) {
        if (dynamic_placement) {
          place_qubit(list_gate, layer, false);
        } else {
          mapping.emplace_back(qubit_mapping[0]);
        }
        if (layer + 1 < list_gate.size()) {
          place_gate(mapping, list_gate, layer + 1, false);
        }
        if (!reuse_qubits[layer].empty()) {
          if (dynamic_placement) {
            place_qubit(list_gate, layer, true);
          } else {
            mapping.emplace_back(qubit_mapping[0]);
            for (const auto q : reuse_qubits[layer]) {
              mapping[mapping.size() - 1][q] = mapping[mapping.size() - 4][q]; // todo: check if this is correct
            }
          }
          if (layer + 1 < list_gate.size()) {
            place_gate({mapping[mapping.size() - 4], mapping[mapping.size() - 1]}, list_gate, layer + 1, true);
            filter_mapping(layer);
          }
        }
      }
      std::cout << "[INFO] ZAC: Minimum-weight-full-matching-based "
                   "intermediate placement: Finish\n";
    }

    auto filter_mapping(const std::size_t layer) -> void {
      const auto& last_gate_mapping = mapping[mapping.size() - 5];
      auto& qubit_mapping = mapping[mapping.size() - 4];
      auto& gate_mapping = mapping[mapping.size() - 3];
      double cost_no_reuse = 0;
      std::unordered_map<
          std::tuple<const SLM*, std::size_t, const SLM*, std::size_t>, double>
          movement_parallel_movement_1;
      std::unordered_map<
          std::tuple<const SLM*, std::size_t, const SLM*, std::size_t>, double>
          movement_parallel_movement_2;

      for (std::size_t q = 0; q < last_gate_mapping.size(); ++q) {
        if (last_gate_mapping[q] != gate_mapping[q]) {
          const auto slm_idx1 =
              std::get<0>(last_gate_mapping[q])->entanglement_id->front().get();
          const auto slm_idx2 =
              std::get<0>(qubit_mapping[q])->entanglement_id->front().get();
          const std::tuple key{slm_idx1, std::get<1>(last_gate_mapping[q]),
                               slm_idx2, std::get<1>(qubit_mapping[q])};
          const double dis = architecture.distance(
              std::get<0>(last_gate_mapping[q]),
              std::get<1>(last_gate_mapping[q]),
              std::get<2>(last_gate_mapping[q]), std::get<0>(qubit_mapping[q]),
              std::get<1>(qubit_mapping[q]), std::get<2>(qubit_mapping[q]));
          if (const auto it = movement_parallel_movement_1.find(key);
              it != movement_parallel_movement_1.end()) {
            movement_parallel_movement_1[key] = std::max(it->second, dis);
          } else {
            movement_parallel_movement_1[key] = dis;
          }
        }
        if (qubit_mapping[q] != gate_mapping[q]) {
          const auto slm_idx1 =
              std::get<0>(gate_mapping[q])->entanglement_id->front().get();
          const auto slm_idx2 =
              std::get<0>(qubit_mapping[q])->entanglement_id->front().get();
          const std::tuple key{slm_idx2, std::get<1>(qubit_mapping[q]),
                               slm_idx1, std::get<1>(gate_mapping[q])};
          const double dis = architecture.distance(
              std::get<0>(qubit_mapping[q]), std::get<1>(qubit_mapping[q]),
              std::get<2>(qubit_mapping[q]), std::get<0>(gate_mapping[q]),
              std::get<1>(gate_mapping[q]), std::get<2>(gate_mapping[q]));
          if (const auto it = movement_parallel_movement_2.find(key);
              it != movement_parallel_movement_2.end()) {
            movement_parallel_movement_2[key] = std::max(it->second, dis);
          } else {
            movement_parallel_movement_2[key] = dis;
          }
        }
      }
      for (const auto& [_, value] : movement_parallel_movement_1) {
        cost_no_reuse += std::sqrt(value);
      }
      for (const auto& [_, value] : movement_parallel_movement_2) {
        cost_no_reuse += std::sqrt(value);
      }
      // now calculate cost with reuse
      gate_mapping = mapping[mapping.size() - 1];
      qubit_mapping = mapping[mapping.size() - 2];
      double cost_reuse = 0;
      movement_parallel_movement_1.clear();
      movement_parallel_movement_2.clear();
      for (std::size_t q = 0; q < last_gate_mapping.size(); ++q) {
        if (last_gate_mapping[q] != gate_mapping[q]) {
          const auto slm_idx1 =
              std::get<0>(last_gate_mapping[q])->entanglement_id->front().get();
          const auto slm_idx2 =
              std::get<0>(qubit_mapping[q])->entanglement_id->front().get();
          const std::tuple key{slm_idx1, std::get<1>(last_gate_mapping[q]),
                               slm_idx2, std::get<1>(qubit_mapping[q])};
          double dis = architecture.distance(
              std::get<0>(last_gate_mapping[q]),
              std::get<1>(last_gate_mapping[q]),
              std::get<2>(last_gate_mapping[q]), std::get<0>(qubit_mapping[q]),
              std::get<1>(qubit_mapping[q]), std::get<2>(qubit_mapping[q]));
          if (const auto it = movement_parallel_movement_1.find(key);
              it != movement_parallel_movement_1.end()) {
            movement_parallel_movement_1[key] =
                std::max(movement_parallel_movement_1[key], dis);
          } else {
            movement_parallel_movement_1[key] = dis;
          }
        }
        if (qubit_mapping[q] != gate_mapping[q]) {
          const auto slm_idx1 =
              std::get<0>(gate_mapping[q])->entanglement_id->front().get();
          const auto slm_idx2 =
              std::get<0>(qubit_mapping[q])->entanglement_id->front().get();
          const std::tuple key{slm_idx2, std::get<1>(qubit_mapping[q]),
                               slm_idx1, std::get<1>(gate_mapping[q])};
          const double dis = architecture.distance(
              std::get<0>(qubit_mapping[q]), std::get<1>(qubit_mapping[q]),
              std::get<2>(qubit_mapping[q]), std::get<0>(gate_mapping[q]),
              std::get<1>(gate_mapping[q]), std::get<2>(gate_mapping[q]));
          if (const auto it = movement_parallel_movement_2.find(key);
              it != movement_parallel_movement_2.end()) {
            movement_parallel_movement_2[key] =
                std::max(movement_parallel_movement_2[key], dis);
          } else {
            movement_parallel_movement_2[key] = dis;
          }
        }
      }
      for (const auto& [_, value] : movement_parallel_movement_1) {
        cost_reuse += std::sqrt(value);
      }
      for (const auto& [_, value] : movement_parallel_movement_2) {
        cost_reuse += std::sqrt(value);
      }
      if (costAtomTransfer * pow((1 - cost_no_reuse / 1.5e6), nQubit) >
          pow((1 - cost_reuse / 1.5e6), nQubit)) {
        listReuseQubits[layer] = {};
        mapping.pop_back();
        mapping.pop_back();
      } else {
        mapping.erase(mapping.end() - 3);
        mapping.erase(mapping.end() - 3);
      }
    }
    /// generate gate mapping based on minimum weight matching
    /// @param list_qubit_mapping the initial mapping before the Rydberg stage
    /// @param list_two_gate_layer gates to be executed in the current Rydberg
    /// stage
    auto place_gate(
        const std::vector<std::vector<std::tuple<const SLM*, size_t, size_t>>>&
            list_qubit_mapping,
        const std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>&
            list_two_gate_layer,
        const std::size_t layer, const bool test_reuse) -> void {
      const auto& list_gate = list_two_gate_layer[layer];
      std::unordered_map<std::size_t, std::size_t> dict_reuse_qubit_neighbor;
      if (list_two_gate_layer.size() > layer + 1 and test_reuse) {
        for (const auto q : listReuseQubits[layer]) {
          for (const auto gate : list_two_gate_layer[layer + 1]) {
            if (q == gate->first) {
              dict_reuse_qubit_neighbor[q] = gate->second;
              break;
            }
            if (q == gate->second) {
              dict_reuse_qubit_neighbor[q] = gate->first;
              break;
            }
          }
        }
      }

      const auto& qubit_mapping = list_qubit_mapping[layer > 0 ? 1 : 0];
      const std::vector<std::tuple<const SLM*, size_t, size_t>>& gate_mapping =
          layer > 0 ? list_qubit_mapping[0]
                    : std::vector<std::tuple<const SLM*, size_t, size_t>>{};
      std::unordered_map<std::tuple<const SLM*, std::size_t, std::size_t>,
                         std::size_t>
          site_Rydberg_to_idx;
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>
          list_Rydberg;
      std::vector<std::size_t> list_row_coo;
      std::vector<std::size_t> list_col_coo;
      std::vector<std::size_t> list_data;
      const auto expand_factor =
          static_cast<std::size_t>(std::ceil(std::sqrt(list_gate.size() / 2)));
      for (std::size_t i = 0; i < list_gate.size(); ++i) {
        const auto& [q1, q2] = *list_gate[i];
        std::unordered_set<std::tuple<const SLM*, std::size_t, std::size_t>>
            set_nearby_site;

        if (test_reuse && (listReuseQubits[layer - 1].find(q1) !=
                           listReuseQubits[layer - 1].end())) {
          const auto location = gate_mapping[q1];
          const auto slm_idx =
              std::get<0>(location)->entanglement_id->front().get();
          set_nearby_site.emplace(slm_idx, std::get<1>(location),
                                  std::get<2>(location));
        } else if (test_reuse && (listReuseQubits[layer - 1].find(q2) !=
                                  listReuseQubits[layer - 1].end())) {
          const auto location = gate_mapping[q2];
          const auto slm_idx =
              std::get<0>(location)->entanglement_id->front().get();
          set_nearby_site.emplace(slm_idx, std::get<1>(location),
                                  std::get<2>(location));
        } else {
          const auto* slm = std::get<0>(qubit_mapping[q1]);
          std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>
              list_nearest_site{};
          list_nearest_site.emplace_back(architecture.nearest_entanglement_site(
              std::get<0>(qubit_mapping[q1]), std::get<1>(qubit_mapping[q1]),
              std::get<2>(qubit_mapping[q1]), std::get<0>(qubit_mapping[q2]),
              std::get<1>(qubit_mapping[q2]), std::get<2>(qubit_mapping[q2])));
          list_nearest_site.emplace_back(architecture.nearest_entanglement_site(
              std::get<0>(qubit_mapping[q1]), 0, std::get<2>(qubit_mapping[q1]),
              std::get<0>(qubit_mapping[q2]), 0,
              std::get<2>(qubit_mapping[q2])));
          list_nearest_site.emplace_back(architecture.nearest_entanglement_site(
              std::get<0>(qubit_mapping[q1]), slm->n_r - 1,
              std::get<2>(qubit_mapping[q1]), std::get<0>(qubit_mapping[q2]),
              slm->n_r - 1, std::get<2>(qubit_mapping[q2])));
          for (const auto& nearest_site : list_nearest_site) {
            set_nearby_site.emplace(nearest_site);
            const auto& [slm_idx, slm_r, slm_c] = nearest_site;
            auto low_r = slm_r > expand_factor ? slm_r - expand_factor : 0;
            auto high_r = std::min(slm->n_r, slm_r + expand_factor + 1);
            auto low_c = slm_c > expand_factor ? slm_c - expand_factor : 0;
            auto high_c = std::min(slm->n_c, slm_c + expand_factor + 1);
            if (high_c - low_c < 2 * expand_factor) {
              const auto height_gap = static_cast<std::size_t>(
                  std::ceil(list_gate.size() / (high_c - low_c)) -
                  expand_factor);
              low_r = low_r > (height_gap / 2) ? low_r - (height_gap / 2) : 0;
              high_r = std::min(slm->n_r, low_r + height_gap + expand_factor);
            }
            if (high_r - low_r < 2 * expand_factor) {
              const auto width_gap = static_cast<std::size_t>(
                  std::ceil(list_gate.size() / (high_r - low_r)) -
                  expand_factor);
              low_c = low_c > (width_gap / 2) ? low_c - (width_gap / 2) : 0;
              high_c = std::min(slm->n_c, low_c + width_gap + expand_factor);
            }
            for (std::size_t r = low_r; r < high_r; ++r) {
              for (std::size_t c = low_c; c < high_c; ++c) {
                set_nearby_site.emplace(slm_idx, r, c);
              }
            }
          }
        }
        for (const auto& site : set_nearby_site) {
          if (site_Rydberg_to_idx.find(site) == site_Rydberg_to_idx.end()) {
            site_Rydberg_to_idx[site] = list_Rydberg.size();
            list_Rydberg.emplace_back(site);
          }
          const auto idx_rydberg = site_Rydberg_to_idx[site];
          double dis1 = architecture.distance(
              std::get<0>(qubit_mapping[q1]), std::get<1>(qubit_mapping[q1]),
              std::get<2>(qubit_mapping[q1]), std::get<0>(site),
              std::get<1>(site), std::get<2>(site));
          double dis2 = architecture.distance(
              std::get<0>(qubit_mapping[q2]), std::get<1>(qubit_mapping[q2]),
              std::get<2>(qubit_mapping[q2]), std::get<0>(site),
              std::get<1>(site), std::get<2>(site));
          double dis3 = 0;
          std::optional<std::size_t> q3;
          if (dict_reuse_qubit_neighbor.find(q1) !=
              dict_reuse_qubit_neighbor.end()) {
            q3 = dict_reuse_qubit_neighbor[q1];
          } else if (dict_reuse_qubit_neighbor.find(q2) !=
                     dict_reuse_qubit_neighbor.end()) {
            q3 = dict_reuse_qubit_neighbor[q2];
          }
          if (q3) {
            dis3 = architecture.distance(std::get<0>(qubit_mapping[*q3]),
                                         std::get<1>(qubit_mapping[*q3]),
                                         std::get<2>(qubit_mapping[*q3]),
                                         std::get<0>(site), std::get<1>(site),
                                         std::get<2>(site));
          }
          list_row_coo.emplace_back(idx_rydberg);
          list_col_coo.emplace_back(i);
          if (std::get<1>(qubit_mapping[q1]) ==
                  std::get<1>(qubit_mapping[q2]) &&
              std::get<0>(qubit_mapping[q1]) ==
                  std::get<0>(qubit_mapping[q2])) {
            list_data.emplace_back(std::sqrt(std::max(dis1, dis2)) +
                                   std::sqrt(dis3));
          } else {
            list_data.emplace_back(std::sqrt(dis1) + std::sqrt(dis2) +
                                   std::sqrt(dis3));
          }
        }
      }

      if (list_Rydberg.size() < list_gate.size()) {
        std::stringstream ss;
        ss << "[Error] ZAC: Minimum-weight-full-matching-based intermediate "
              "placement: No enough sites for gates ("
           << list_Rydberg.size() << " vs " << list_gate.size() << ").";
        throw std::invalid_argument(ss.str());
      }

      std::vector matrix(
          list_Rydberg.size(),
          std::vector<std::optional<double>>(list_gate.size(), std::nullopt));
      for (std::size_t i = 0; i < list_row_coo.size(); ++i) {
        matrix[list_row_coo[i]][list_col_coo[i]] = list_data[i];
      }
      const auto& matching = minimumWeightFullBipartiteMatching(matrix);
      double cost = 0;
      for (std::size_t i = 0; i < matching.size(); ++i) {
        cost = cost +  matrix[i][matching[i]].value();
      }
      auto tmp_mapping = qubit_mapping;
      for (std::size_t idx_rydberg = 0; idx_rydberg < matching.size();
           ++idx_rydberg) {
        const auto q0 = list_gate[matching[idx_rydberg]]->first;
        const auto q1 = list_gate[matching[idx_rydberg]]->second;
        const auto& site = list_Rydberg[idx_rydberg];
        if (test_reuse && (listReuseQubits[layer - 1].find(q0) !=
                           listReuseQubits[layer - 1].end())) {
          tmp_mapping[q0] = gate_mapping[q0];
          if (site == gate_mapping[q0]) {
            tmp_mapping[q1] = {std::get<0>(site)->entanglement_id->back().get(),
                               std::get<1>(site), std::get<2>(site)};
          } else {
            tmp_mapping[q1] = site;
          }
        } else if (test_reuse && (listReuseQubits[layer - 1].find(q1) !=
                                  listReuseQubits[layer - 1].end())) {
          tmp_mapping[q1] = gate_mapping[q1];
          if (site == gate_mapping[q1]) {
            tmp_mapping[q0] = {std::get<0>(site)->entanglement_id->back().get(),
                               std::get<1>(site), std::get<2>(site)};
          } else {
            tmp_mapping[q0] = site;
          }
        } else {
          if (std::get<2>(qubit_mapping[q0]) < std::get<2>(qubit_mapping[q1])) {
            tmp_mapping[q0] = site;
            tmp_mapping[q1] = {std::get<0>(site)->entanglement_id->back().get(),
                               std::get<1>(site), std::get<2>(site)};
          } else {
            tmp_mapping[q0] = {std::get<0>(site)->entanglement_id->back().get(),
                               std::get<1>(site), std::get<2>(site)};
            tmp_mapping[q1] = site;
          }
        }
      }
      mapping.emplace_back(tmp_mapping);
    }

    /// Generate qubit mapping based on minimum weight matching.
    /// @param list_gate List of 2-qubit gates as pairs of qubits
    /// - @code list_gate[layer]@endcode: gates to be executed in the current
    ///   Rydberg stage
    /// - @code list_gate[i], i > layer@endcode: is the list of yet unexecuted
    ///   gates
    /// @param layer the current Rydberg stage
    /// @param test_reuse whether to test the reuse of qubits or not
    auto place_qubit(
        const std::vector<std::vector<const std::pair<qc::Qubit, qc::Qubit>*>>&
            list_gate,
        const std::size_t layer, const bool test_reuse) -> void {
      // the very initial placement of qubits
      const auto& qubit_mapping = mapping.front();
      // the placement of qubits after the last gate
      const auto& last_gate_mapping =
          mapping[mapping.size() - (test_reuse ? 3 : 1)];
      // for each storage SLM array, construct a matrix indicating occupancy of
      // sites
      std::unordered_map<const SLM*, std::vector<std::vector<bool>>>
          is_empty_storage_site{};
      // for each SLM array, initialize the site as empty
      for (const std::unique_ptr<SLM>& slm_id :
           architecture.storage_zone) {
        is_empty_storage_site.emplace(
            slm_id.get(),
            std::vector(slm_id->n_r, std::vector(slm_id->n_c, true)));
      }
      // qubits that need to be placed, in particular this does not include
      // qubits that are reused NOTE: the indices stored in qubit_to_place are
      // the indices of the last_gate_mapping and do not refer to the actual
      // index of a qubit
      std::vector<qc::Qubit> qubit_to_place{};
      // go through the placement of qubits after the last gate
      for (std::size_t q = 0; q < last_gate_mapping.size(); ++q) {
        const auto& last_mapping = last_gate_mapping[q];
        const auto array_id = std::get<0>(last_mapping);
        if (is_empty_storage_site.find(array_id) !=
            is_empty_storage_site.end()) {
          // the mapped qubit is in the storage zone, set the site as occupied
          is_empty_storage_site[array_id][std::get<1>(last_mapping)]
                               [std::get<2>(last_mapping)] = false;
        } else if (!test_reuse || listReuseQubits[layer].find(q) ==
                                      listReuseQubits[layer].end()) {
          // the mapped qubit is in the entangling zone and must be placed (if
          // not reused)
          qubit_to_place.emplace_back(q);
        }
      }
      // occupied sites in the initial mapping that are currently unoccupied;
      // those are used as candidate sites for qubit placement as it is
      // guaranteed that they exist
      std::unordered_set<std::tuple<const SLM*, std::size_t, std::size_t>>
          common_site{};
      for (const std::tuple<const SLM*, std::size_t, std::size_t>& m :
           mapping[0]) {
        const auto* array_id = std::get<0>(m);
        if (is_empty_storage_site.find(array_id) !=
            is_empty_storage_site.end()) {
          if (is_empty_storage_site[array_id][std::get<1>(m)][std::get<2>(m)]) {
            // add site to common sites if it is currently unoccupied
            common_site.emplace(array_id, std::get<1>(m), std::get<2>(m));
          }
        } else {
          // site is in an entangling zone and not yet added to
          // is_empty_storage_site
          is_empty_storage_site[array_id] =
              std::vector(array_id->n_r, std::vector(array_id->n_c, true));
          common_site.emplace(array_id, std::get<1>(m), std::get<2>(m));
        }
      }
      // dictionary to store the qubit interactions with other qubits
      std::unordered_map<qc::Qubit, std::vector<qc::Qubit>>
          dict_qubit_interaction{};
      for (const auto q : qubit_to_place) {
        dict_qubit_interaction.emplace(q, std::vector<qc::Qubit>{});
      }
      if (list_gate.size() > layer + 1) {
        for (const auto* gate : list_gate[layer + 1]) {
          if (dict_qubit_interaction.find(gate->first) !=
                  dict_qubit_interaction.end() &&
              ((!test_reuse) || (listReuseQubits[layer].find(gate->second) ==
                                 listReuseQubits[layer].end()))) {
            dict_qubit_interaction[gate->first].emplace_back(gate->second);
          }
          if (dict_qubit_interaction.find(gate->second) !=
                  dict_qubit_interaction.end() &&
              ((!test_reuse) || (listReuseQubits[layer].find(gate->first) ==
                                 listReuseQubits[layer].end()))) {
            dict_qubit_interaction[gate->second].emplace_back(gate->first);
          }
        }
      }
      const std::size_t expand_factor = 1;

      std::unordered_map<std::tuple<const SLM*, std::size_t, std::size_t>,
                         std::size_t>
          site_storage_to_idx{};
      std::vector<std::tuple<const SLM*, std::size_t, std::size_t>>
          list_storage{};
      std::vector<std::size_t> list_col_coo{};
      std::vector<std::size_t> list_row_coo{};
      std::vector<double> list_data{};

      for (std::size_t i = 0; i < qubit_to_place.size(); ++i) {
        const auto q = qubit_to_place[i];
        std::unordered_map<const SLM*, std::tuple<std::size_t, std::size_t,
                                                  std::size_t, std::size_t>>
            dict_bouding_box{};
        const auto* slm = std::get<0>(qubit_mapping[q]);
        auto lower_row = std::get<1>(qubit_mapping[q]);
        auto upper_row = lower_row;
        auto left_col = std::get<2>(qubit_mapping[q]);
        auto right_col = left_col;
        const std::pair<size_t, size_t>& exact_loc_q = std::apply(
            [&](auto&&... args) {
              return architecture.exact_SLM_location(
                  std::forward<decltype(args)>(args)...);
            },
            qubit_mapping[q]);
        const std::pair<size_t, size_t>& exact_loc_gate =
          std::apply([&](auto&&... args) {
            return architecture.exact_SLM_location(
                std::forward<decltype(args)>(args)...);
          }, mapping[0][q]);
        if (exact_loc_gate.second < exact_loc_q.second) {
          lower_row = 0;
        } else {
          upper_row = slm->n_r;
        }

        dict_bouding_box.emplace(
            std::get<0>(qubit_mapping[q]),
            std::tuple{lower_row, upper_row, left_col, right_col});
        for (const qc::Qubit neighbor_q : dict_qubit_interaction[q]) {
          const SLM* tmp_slm_idx = std::get<0>(last_gate_mapping[neighbor_q]);
          const std::tuple<const SLM*, size_t, size_t>& neighbor_q_location =
              tmp_slm_idx->isEntanglement()
          ? std::apply([&](auto&&... args) -> std::tuple<const SLM*, size_t, size_t> {
      return architecture.nearest_storage_site(
          std::forward<decltype(args)>(args)...);
    },
                               last_gate_mapping[neighbor_q])
                  : last_gate_mapping[neighbor_q];
          if (const auto& it =
                  dict_bouding_box.find(std::get<0>(neighbor_q_location));
              it != dict_bouding_box.end()) {
            dict_bouding_box.emplace(
                std::get<0>(neighbor_q_location),
                std::tuple{std::min(std::get<1>(neighbor_q_location),
                                    std::get<0>(it->second)),
                           std::max(std::get<1>(neighbor_q_location),
                                    std::get<1>(it->second)),
                           std::min(std::get<2>(neighbor_q_location),
                                    std::get<2>(it->second)),
                           std::max(std::get<2>(neighbor_q_location),
                                    std::get<3>(it->second))});
          } else {
            const auto* slm_id = std::get<0>(neighbor_q_location);
            lower_row = std::get<1>(neighbor_q_location);
            upper_row = std::get<1>(neighbor_q_location);
            left_col = std::get<2>(neighbor_q_location);
            right_col = std::get<2>(neighbor_q_location);
            const std::pair<std::size_t, std::size_t>& exact_loc_neightbor_q =
              std::apply([&](auto&&... args) {
    return architecture.exact_SLM_location(
        std::forward<decltype(args)>(args)...);
  },
                           neighbor_q_location);
            if (exact_loc_gate.second < exact_loc_neightbor_q.second) {
              lower_row = 0;
            } else {
              upper_row = slm_id->n_r;
            }
            dict_bouding_box.emplace(
                std::get<0>(neighbor_q_location),
                std::tuple{lower_row, upper_row, left_col, right_col});
          }
        }
        const auto& gate_location = last_gate_mapping[q];
        const std::tuple<const SLM*, std::size_t, std::size_t>
            nearest_storage_site =
              std::apply([&](auto&&... args) {
  return architecture.nearest_storage_site(
      std::forward<decltype(args)>(args)...);
}, gate_location);
        // todo: what is ratio?
        const std::size_t ratio = 3;
        if (const auto it =
                dict_bouding_box.find(std::get<0>(nearest_storage_site));
            it != dict_bouding_box.end()) {
          dict_bouding_box.emplace(
              std::get<0>(nearest_storage_site),
              std::tuple{std::min(std::get<1>(nearest_storage_site) - ratio,
                                  std::get<0>(it->second)),
                         std::max(std::get<1>(nearest_storage_site) + ratio,
                                  std::get<1>(it->second)),
                         std::min(std::get<2>(nearest_storage_site) - ratio,
                                  std::get<2>(it->second)),
                         std::max(std::get<2>(nearest_storage_site) + ratio,
                                  std::get<3>(it->second))});
        } else {
          dict_bouding_box.emplace(
              std::get<0>(nearest_storage_site),
              std::tuple{std::get<1>(nearest_storage_site) - ratio,
                         std::get<1>(nearest_storage_site) + ratio,
                         std::get<2>(nearest_storage_site) - ratio,
                         std::get<2>(nearest_storage_site) + ratio});
        }
        auto set_nearby_site = common_site;
        if (is_empty_storage_site[std::get<0>(qubit_mapping[q])][std::get<1>(
                qubit_mapping[q])][std::get<2>(qubit_mapping[q])]) {
          set_nearby_site.emplace(qubit_mapping[q]);
        }

        for (const auto& [slm_id, _] : dict_bouding_box) {
          auto& bounding_box = dict_bouding_box[slm_id];
          bounding_box = {
              std::get<0>(bounding_box) > expand_factor
                  ? static_cast<std::size_t>(std::get<0>(bounding_box) -
                                             expand_factor)
                  : 0,
              std::min(std::get<1>(bounding_box) + expand_factor + 1,
                       slm_id->n_r),
              std::get<2>(bounding_box) > expand_factor
                  ? static_cast<std::size_t>(std::get<2>(bounding_box) -
                                             expand_factor)
                  : 0,
              std::min(std::get<3>(bounding_box) + expand_factor + 1,
                       slm_id->n_c)};
          for (auto r = std::get<0>(bounding_box);
               r < std::get<1>(bounding_box); ++r) {
            for (auto c = std::get<2>(bounding_box);
                 c < std::get<3>(bounding_box); ++c) {
              if (is_empty_storage_site.find(slm_id) ==
                      is_empty_storage_site.end() ||
                  is_empty_storage_site[slm_id][r][c]) {
                set_nearby_site.emplace(slm_id, r, c);
              }
            }
          }
        }

        for (const auto& site : set_nearby_site) {
          if (site_storage_to_idx.find(site) == site_storage_to_idx.end()) {
            site_storage_to_idx.emplace(site, list_storage.size());
            list_storage.emplace_back(site);
          }
          const auto idx_storage = site_storage_to_idx[site];
          const double dis = architecture.distance(
              std::get<0>(gate_location), std::get<1>(gate_location),
              std::get<2>(gate_location), std::get<0>(site), std::get<1>(site),
              std::get<2>(site));
          double lookahead_cost = 0;
          for (const auto& neighbor_q : dict_qubit_interaction[q]) {
            const auto& site_neighbor_q = last_gate_mapping[neighbor_q];
            if (!std::get<0>(site_neighbor_q)->entanglement_id) {
              lookahead_cost += architecture.nearest_entanglement_site_distance(
                  std::get<0>(site), std::get<1>(site), std::get<2>(site),
                  std::get<0>(site_neighbor_q), std::get<1>(site_neighbor_q),
                  std::get<2>(site_neighbor_q));
            } else {
              const std::pair<std::size_t, std::size_t>& exact_loc_neightbor_q =
                std::apply([&](auto&&... args) {
  return architecture.exact_SLM_location(
      std::forward<decltype(args)>(args)...);
},
                             last_gate_mapping[neighbor_q]);
              const auto dx =
                  static_cast<std::int64_t>(exact_loc_neightbor_q.first) -
                  static_cast<std::int64_t>(exact_loc_q.first);
              const auto dy =
                  static_cast<std::int64_t>(exact_loc_neightbor_q.second) -
                  static_cast<std::int64_t>(exact_loc_q.second);
              lookahead_cost += std::sqrt(std::sqrt((dx * dx) + (dy * dy)));
            }
          }
          const double cost = std::sqrt(dis) + (0.1 * lookahead_cost);
          list_col_coo.emplace_back(idx_storage);
          list_row_coo.emplace_back(i);
          list_data.emplace_back(cost);
        }

        // to construct from three arrays:
        // - list_data the entries of the matrix, in any order
        // - list_row_coo the row indices of the matrix entries
        // - list_col_coo the column indices of the matrix entries
        // Where A[list_row_coo[k], list_col_coo[k]] = list_data[k].
        std::vector cost_matrix(qubit_to_place.size(),
                                std::vector<std::optional<double>>(
                                    list_storage.size(), std::nullopt));
        for (std::size_t k = 0; k < list_data.size(); ++k) {
          cost_matrix[list_row_coo[k]][list_col_coo[k]] = list_data[k];
        }
        const auto& matching = minimumWeightFullBipartiteMatching(cost_matrix);
        auto tmp_mapping = last_gate_mapping;
        for (std::size_t j = 0; j < matching.size(); ++j) {
          tmp_mapping[qubit_to_place[j]] = list_storage[matching[j]];
        }
        mapping.emplace_back(tmp_mapping);
      }
    }
  };
};

} // namespace na
