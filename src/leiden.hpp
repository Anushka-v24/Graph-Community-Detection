// leiden.hpp — the Leiden method for community detection.
//
// Addresses Louvain's core limitation: Louvain can strand internally
// disconnected or poorly connected communities during aggregation.
//
// Leiden solves this with a 3-phase cycle:
//   1. Local moving: Optimize modularity to obtain a partition P.
//   2. Refinement: Refine P into sub-partition P_refined where each sub-community
//      is guaranteed to be well-connected and a subset of a community in P.
//   3. Aggregation: Coarsen the graph based on P_refined, while preserving
//      the community assignments from P for the aggregated nodes.
//
// Reference: Traag, Waltman, van Eck (2019) "From Louvain to Leiden: guaranteeing
// well-connected communities", Scientific Reports 9:5233.

#pragma once

#include "csr.hpp"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <unordered_map>
#include <vector>
#include <queue>

using namespace std;

namespace gcd {

namespace detail_leiden {

// Fast local moving phase on graph g.
// Starts with initial partition 'comm' (either singletons or inherited from
// previous refinement) and sweeps until convergence.
inline bool local_moving_pass(const CSR& g, std::vector<node_t>& comm, double resolution) {
    double two_m = 2.0 * g.total_weight;
    if (two_m <= 0.0) return false;

    node_t max_c = 0;
    for (node_t v = 0; v < g.n; ++v) max_c = std::max(max_c, comm[v]);
    std::vector<double> sigma_tot(max_c + 1, 0.0);
    for (node_t v = 0; v < g.n; ++v) sigma_tot[comm[v]] += g.weighted_degree[v];

    bool moved_any = false;

    // Fast active queue: start with all nodes
    std::vector<bool> in_queue(g.n, true);
    std::queue<node_t> q;
    for (node_t v = 0; v < g.n; ++v) q.push(v);

    std::unordered_map<node_t, double> k_in;

    while (!q.empty()) {
        node_t u = q.front();
        q.pop();
        in_queue[u] = false;

        node_t cu = comm[u];
        double ku = g.weighted_degree[u];

        k_in.clear();
        for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
            node_t v = g.neighbors[i];
            if (v == u) continue;
            k_in[comm[v]] += g.weights[i];
        }

        sigma_tot[cu] -= ku;

        node_t best = cu;
        auto it_self = k_in.find(cu);
        double best_gain = (it_self != k_in.end() ? it_self->second : 0.0)
                          - resolution * sigma_tot[cu] * ku / two_m;

        for (const auto& kv : k_in) {
            node_t c = kv.first;
            if (c == cu) continue;
            if (c >= sigma_tot.size()) {
                sigma_tot.resize(c + 1, 0.0);
            }
            double gain = kv.second - resolution * sigma_tot[c] * ku / two_m;
            if (gain > best_gain) {
                best_gain = gain;
                best = c;
            }
        }

        if (best >= sigma_tot.size()) {
            sigma_tot.resize(best + 1, 0.0);
        }
        sigma_tot[best] += ku;

        if (best != cu) {
            comm[u] = best;
            moved_any = true;
            // Enqueue all neighbors not already in queue
            for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
                node_t v = g.neighbors[i];
                if (v != u && !in_queue[v]) {
                    q.push(v);
                    in_queue[v] = true;
                }
            }
        }
    }

    return moved_any;
}

// Refinement phase:
// Starts with each node in its own singleton refined community (P_refined[v] = v).
// Merges nodes within each community C of P_fast only if well-connected.
inline std::vector<node_t> refine_partition(const CSR& g,
                                            const std::vector<node_t>& p_fast,
                                            double resolution) {
    double two_m = 2.0 * g.total_weight;
    std::vector<node_t> p_refined(g.n);
    std::iota(p_refined.begin(), p_refined.end(), 0);

    if (two_m <= 0.0) return p_refined;

    // Group nodes by community in p_fast
    node_t max_p = 0;
    for (node_t v = 0; v < g.n; ++v) max_p = std::max(max_p, p_fast[v]);
    std::vector<std::vector<node_t>> comm_nodes(max_p + 1);
    for (node_t v = 0; v < g.n; ++v) {
        comm_nodes[p_fast[v]].push_back(v);
    }

    // Track total degree of refined communities
    std::vector<double> sigma_tot_refined(g.n);
    for (node_t v = 0; v < g.n; ++v) {
        sigma_tot_refined[v] = g.weighted_degree[v];
    }

    std::unordered_map<node_t, double> k_in_refined;

    for (node_t c = 0; c <= max_p; ++c) {
        const auto& nodes = comm_nodes[c];
        if (nodes.size() <= 1) continue;

        // Total degree within this P_fast community
        double total_c_degree = 0.0;
        for (node_t u : nodes) total_c_degree += g.weighted_degree[u];

        for (node_t u : nodes) {
            double ku = g.weighted_degree[u];
            node_t r_u = p_refined[u];

            k_in_refined.clear();
            for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
                node_t v = g.neighbors[i];
                if (v == u) continue;
                // Only consider neighbors belonging to the SAME community in P_fast
                if (p_fast[v] == c) {
                    k_in_refined[p_refined[v]] += g.weights[i];
                }
            }

            // Remove u from its current refined sub-community
            sigma_tot_refined[r_u] -= ku;

            node_t best_r = r_u;
            auto it_self = k_in_refined.find(r_u);
            double best_gain = (it_self != k_in_refined.end() ? it_self->second : 0.0)
                              - resolution * sigma_tot_refined[r_u] * ku / two_m;

            for (const auto& kv : k_in_refined) {
                node_t r_cand = kv.first;
                if (r_cand == r_u) continue;

                // Candidate sub-community well-connectedness check
                double k_cand = sigma_tot_refined[r_cand];
                double well_connected_thresh = resolution * (k_cand * (total_c_degree - k_cand)) / two_m;
                if (kv.second >= well_connected_thresh || kv.second > 0) {
                    double gain = kv.second - resolution * k_cand * ku / two_m;
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_r = r_cand;
                    }
                }
            }

            sigma_tot_refined[best_r] += ku;
            p_refined[u] = best_r;
        }
    }

    return p_refined;
}

// Renumber community ids into a dense 0..k-1 range
inline node_t renumber(std::vector<node_t>& comm) {
    std::unordered_map<node_t, node_t> remap;
    for (auto& c : comm) {
        auto it = remap.find(c);
        if (it == remap.end()) {
            node_t id = static_cast<node_t>(remap.size());
            remap.emplace(c, id);
            c = id;
        } else {
            c = it->second;
        }
    }
    return static_cast<node_t>(remap.size());
}

// Aggregation by refined partition P_refined, mapping each refined node to its parent P_fast
inline CSR aggregate_leiden(const CSR& g,
                            const std::vector<node_t>& p_refined,
                            node_t num_refined,
                            const std::vector<node_t>& p_fast,
                            std::vector<node_t>& comm_next) {
    std::unordered_map<uint64_t, double> acc;
    auto key_of = [](node_t a, node_t b) {
        if (a > b) std::swap(a, b);
        return (uint64_t(a) << 32) | b;
    };

    comm_next.resize(num_refined);

    for (node_t u = 0; u < g.n; ++u) {
        node_t ru = p_refined[u];
        comm_next[ru] = p_fast[u]; // parent community in P_fast
        for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
            node_t rv = p_refined[g.neighbors[i]];
            acc[key_of(ru, rv)] += g.weights[i];
        }
    }

    std::vector<std::pair<node_t, node_t>> edges;
    std::vector<double> weights;
    edges.reserve(acc.size());
    weights.reserve(acc.size());
    for (const auto& kv : acc) {
        node_t a = static_cast<node_t>(kv.first >> 32);
        node_t b = static_cast<node_t>(kv.first & 0xffffffffu);
        edges.emplace_back(a, b);
        weights.push_back(kv.second / 2.0);
    }

    return build_csr(edges, weights, num_refined);
}

}  // namespace detail_leiden

// Run Leiden to convergence.
// Returns comm[v] = final community id for each original node, dense 0..k-1.
inline std::vector<node_t> leiden(const CSR& g, double resolution = 1.0) {
    if (g.n == 0) return {};

    std::vector<node_t> membership(g.n);
    std::iota(membership.begin(), membership.end(), 0);

    CSR level_graph = g;
    std::vector<node_t> level_comm(g.n);
    std::iota(level_comm.begin(), level_comm.end(), 0);

    while (true) {
        // Phase 1: Fast local moving to find P_fast
        bool moved = detail_leiden::local_moving_pass(level_graph, level_comm, resolution);
        detail_leiden::renumber(level_comm);

        // Phase 2: Refinement to find P_refined
        std::vector<node_t> p_refined = detail_leiden::refine_partition(level_graph, level_comm, resolution);
        node_t num_refined = detail_leiden::renumber(p_refined);

        // Update overall node membership using P_refined
        for (size_t v = 0; v < membership.size(); ++v) {
            membership[v] = p_refined[membership[v]];
        }

        // Stop if no nodes moved, or graph collapsed to 1 node per community
        if (!moved || num_refined == level_graph.n) break;

        // Phase 3: Aggregate based on P_refined, mapping initial level_comm from P_fast
        std::vector<node_t> comm_next;
        level_graph = detail_leiden::aggregate_leiden(level_graph, p_refined, num_refined, level_comm, comm_next);
        level_comm = std::move(comm_next);
        detail_leiden::renumber(level_comm);
    }

    // Final dense renumbering of original nodes
    detail_leiden::renumber(membership);
    return membership;
}

}  // namespace gcd
