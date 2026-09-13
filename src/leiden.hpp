// leiden.hpp — the Leiden method for community detection.
//
// Addresses Louvain's core limitation: Louvain can strand internally
// disconnected or poorly connected communities during aggregation.
//
// Leiden solves this with a 3-phase cycle:
//   1. Local moving: Optimize modularity to obtain a partition P_fast.
//   2. Refinement: Refine P_fast into sub-partition P_refined where each sub-community
//      is guaranteed to be well-connected and a subset of a community in P_fast.
//   3. Aggregation: Coarsen the graph based on P_refined, while preserving
//      the community assignments from P_fast for the aggregated nodes.
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

namespace gcd {

namespace detail_leiden {

// Fast local moving phase on graph g.
// Sweeps until convergence using an active queue of nodes.
inline bool local_moving_pass(const CSR& g, std::vector<node_t>& comm, double resolution) {
    double two_m = 2.0 * g.total_weight;
    if (two_m <= 0.0) return false;

    node_t max_c = 0;
    for (node_t v = 0; v < g.n; ++v) max_c = std::max(max_c, comm[v]);
    std::vector<double> sigma_tot(max_c + 1, 0.0);
    for (node_t v = 0; v < g.n; ++v) sigma_tot[comm[v]] += g.weighted_degree[v];

    bool moved_any = false;

    // Active queue: start with all nodes
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
            // Enqueue all neighbors not currently in the active queue
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
// Merges nodes within each community C of P_fast only if well-connected:
//   1. Node u must be well-connected within C:
//      E(u, C \ {u}) >= theta * resolution * k_u * (K(C) - k_u) / 2m
//   2. Candidate sub-community R must be well-connected within C:
//      E(R, C \ R) >= theta * resolution * K(R) * (K(C) - K(R)) / 2m
//   3. Merge must provide strictly positive modularity gain:
//      E(u, R) - resolution * k_u * K(R) / 2m > 0
inline std::vector<node_t> refine_partition(const CSR& g,
                                            const std::vector<node_t>& p_fast,
                                            double resolution,
                                            double theta = 0.05) {
    double two_m = 2.0 * g.total_weight;
    std::vector<node_t> p_refined(g.n);
    std::iota(p_refined.begin(), p_refined.end(), 0);

    if (two_m <= 0.0 || g.n == 0) return p_refined;

    // Group nodes by parent community in p_fast
    node_t max_p = 0;
    for (node_t v = 0; v < g.n; ++v) max_p = std::max(max_p, p_fast[v]);
    std::vector<std::vector<node_t>> comm_nodes(max_p + 1);
    for (node_t v = 0; v < g.n; ++v) {
        comm_nodes[p_fast[v]].push_back(v);
    }

    // Precompute edge weight from each node to other nodes in its own P_fast community: E(u, C \ {u})
    std::vector<double> e_node_in_c(g.n, 0.0);
    for (node_t u = 0; u < g.n; ++u) {
        node_t cu = p_fast[u];
        double e_in = 0.0;
        for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
            node_t v = g.neighbors[i];
            if (v != u && p_fast[v] == cu) {
                e_in += g.weights[i];
            }
        }
        e_node_in_c[u] = e_in;
    }

    // Track total degree and external connectivity E(R, C \ R) for each refined sub-community
    std::vector<double> sigma_tot_refined = g.weighted_degree;
    std::vector<double> e_to_c_minus_r = e_node_in_c;

    std::unordered_map<node_t, double> k_in_refined;

    for (node_t c = 0; c <= max_p; ++c) {
        const auto& nodes = comm_nodes[c];
        if (nodes.size() <= 1) continue;

        // Total degree K(C) within this P_fast community
        double total_c_degree = 0.0;
        for (node_t u : nodes) total_c_degree += g.weighted_degree[u];

        for (node_t u : nodes) {
            double ku = g.weighted_degree[u];

            // 1. Node well-connectedness safety check in C:
            // Node u must have sufficient connection density to C \ {u}
            double node_wc_thresh = theta * resolution * (ku * (total_c_degree - ku)) / two_m;
            if (e_node_in_c[u] < node_wc_thresh && e_node_in_c[u] == 0.0) {
                // Node u is completely disconnected from C \ {u}; isolate as singleton
                continue;
            }

            // Collect edge weights from u to other refined sub-communities inside community C
            k_in_refined.clear();
            for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
                node_t v = g.neighbors[i];
                if (v == u) continue;
                if (p_fast[v] == c) {
                    k_in_refined[p_refined[v]] += g.weights[i];
                }
            }

            node_t best_r = u;
            double best_gain = 0.0; // Must be strictly positive gain to merge

            for (const auto& kv : k_in_refined) {
                node_t r_cand = kv.first;
                if (r_cand == p_refined[u]) continue;

                // 2. Candidate sub-community well-connectedness check in C:
                // R_cand must satisfy E(R, C \ R) >= theta * resolution * K(R) * (K(C) - K(R)) / 2m
                double k_cand = sigma_tot_refined[r_cand];
                double r_wc_thresh = theta * resolution * (k_cand * (total_c_degree - k_cand)) / two_m;
                if (e_to_c_minus_r[r_cand] < r_wc_thresh && e_to_c_minus_r[r_cand] == 0.0) {
                    // Candidate sub-community is completely disconnected in C
                    continue;
                }

                // 3. Positive modularity gain check: Delta Q(u -> R_cand) > 0
                // Requires direct connection E(u, R_cand) > 0 and Delta Q > 0
                if (kv.second > 0.0) {
                    double gain = kv.second - resolution * (ku * k_cand) / two_m;
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_r = r_cand;
                    }
                }
            }

            // Perform well-connected merge if a strictly positive candidate was found
            if (best_r != u && best_gain > 0.0) {
                node_t old_r = p_refined[u];
                p_refined[u] = best_r;

                double e_u_best = k_in_refined[best_r];

                // Update connectivity and degree of best_r:
                // E(R U {u}, C \ (R U {u})) = E(R, C \ R) + E(u, C \ {u}) - 2 * E(u, R)
                e_to_c_minus_r[best_r] += e_node_in_c[u] - 2.0 * e_u_best;
                sigma_tot_refined[best_r] += ku;

                // If old_r was a singleton {u}, its degree is now transferred
                if (old_r == u) {
                    sigma_tot_refined[u] = 0.0;
                    e_to_c_minus_r[u] = 0.0;
                }
            }
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
//
// If hierarchy_out is non-null, it receives membership[v] snapshotted
// after every round (see louvain.hpp — same contract). Note Leiden's
// refinement-based well-connectedness check has a known bug (see
// connectivity.hpp); callers wanting a real connectivity guarantee should
// post-process the result with enforce_connectivity() rather than trust
// this alone.
inline std::vector<node_t> leiden(const CSR& g, double resolution = 1.0,
                                   std::vector<std::vector<node_t>>* hierarchy_out = nullptr) {
    if (g.n == 0) return {};
    if (hierarchy_out) hierarchy_out->clear();

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

        // Stop if no nodes moved, or graph collapsed to 1 node per community
        if (!moved || num_refined == level_graph.n) {
            for (size_t v = 0; v < membership.size(); ++v) {
                membership[v] = level_comm[membership[v]];
            }
            if (hierarchy_out) hierarchy_out->push_back(membership);
            break;
        }

        // Update overall node membership to point to the coarsened nodes for next level
        for (size_t v = 0; v < membership.size(); ++v) {
            membership[v] = p_refined[membership[v]];
        }
        if (hierarchy_out) hierarchy_out->push_back(membership);

        // Phase 3: Aggregate based on P_refined, mapping initial level_comm from P_fast
        std::vector<node_t> comm_next;
        level_graph = detail_leiden::aggregate_leiden(level_graph, p_refined, num_refined, level_comm, comm_next);
        level_comm = std::move(comm_next);
        detail_leiden::renumber(level_comm);
    }

    // Final dense renumbering of original nodes
    detail_leiden::renumber(membership);
    if (hierarchy_out && !hierarchy_out->empty()) hierarchy_out->back() = membership;
    return membership;
}

}  // namespace gcd
