// louvain.hpp — the Louvain method for modularity maximization.
//
// Two nested phases, repeated on a shrinking sequence of graphs until
// modularity stops improving:
//
//   1. Local moving — visit every node, move it into whichever neighboring
//      community increases modularity the most (or leave it if nothing
//      helps). Keep sweeping until a full pass makes no move.
//   2. Aggregation — collapse each community into a single node. Edges
//      between two communities merge into one weighted edge; edges within
//      a community become a self-loop. Recurse local moving on this
//      smaller graph.
//
// The final per-node community is recovered by composing the assignment
// at every level back down to the original nodes.
//
// This is a from-scratch, correctness-first implementation: no OpenMP yet
// (see Makefile), and each level currently copies the graph rather than
// mutating in place. Fine at the sizes Phase 2 tests against; worth
// revisiting once Phase 9's scale runs make the copy show up in a profile.

#pragma once

#include "csr.hpp"

#include <cstdint>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace gcd {

namespace detail {

// One call runs local moving to convergence on graph g. comm[v] is both
// input (starting assignment — normally comm[v] == v, everyone alone) and
// output (final assignment, not yet renumbered to a dense range).
// Returns true if any node ever moved.
inline bool local_moving_pass(const CSR& g, std::vector<node_t>& comm, double resolution) {
    double two_m = 2.0 * g.total_weight;
    if (two_m <= 0.0) return false;

    std::vector<double> sigma_tot(g.n, 0.0);  // total weighted degree of each community
    for (node_t v = 0; v < g.n; ++v) sigma_tot[comm[v]] += g.weighted_degree[v];

    bool moved_any = false;
    bool moved_this_sweep = true;

    // Reused across nodes to avoid reallocating every iteration: weight
    // from the current node into each community it currently touches.
    std::unordered_map<node_t, double> k_in;

    while (moved_this_sweep) {
        moved_this_sweep = false;

        for (node_t u = 0; u < g.n; ++u) {
            node_t cu = comm[u];
            double ku = g.weighted_degree[u];

            k_in.clear();
            for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
                node_t v = g.neighbors[i];
                if (v == u) continue;  // self-loops don't bonus any community; see aggregate()
                k_in[comm[v]] += g.weights[i];
            }

            // Remove u from its own community before scoring candidates —
            // otherwise "stay put" is compared unfairly against itself.
            sigma_tot[cu] -= ku;

            node_t best = cu;
            auto it_self = k_in.find(cu);
            double best_gain = (it_self != k_in.end() ? it_self->second : 0.0)
                              - resolution * sigma_tot[cu] * ku / two_m;

            for (auto& kv : k_in) {
                node_t c = kv.first;
                if (c == cu) continue;
                double gain = kv.second - resolution * sigma_tot[c] * ku / two_m;
                if (gain > best_gain) { best_gain = gain; best = c; }
            }

            sigma_tot[best] += ku;
            if (best != cu) {
                comm[u] = best;
                moved_any = true;
                moved_this_sweep = true;
            }
        }
    }

    return moved_any;
}

// Collapse g down to one node per community (num_comm of them, ids assumed
// already dense 0..num_comm-1). Each original undirected edge is visited
// twice — once from each endpoint's adjacency list, matching CSR's
// double-stored convention — so weights are accumulated doubled and then
// halved back before rebuilding.
inline CSR aggregate(const CSR& g, const std::vector<node_t>& comm, node_t num_comm) {
    std::unordered_map<uint64_t, double> acc;
    auto key_of = [](node_t a, node_t b) {
        if (a > b) std::swap(a, b);
        return (uint64_t(a) << 32) | b;
    };

    for (node_t u = 0; u < g.n; ++u) {
        node_t cu = comm[u];
        for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
            node_t cv = comm[g.neighbors[i]];
            acc[key_of(cu, cv)] += g.weights[i];
        }
    }

    std::vector<std::pair<node_t, node_t>> edges;
    std::vector<double> weights;
    edges.reserve(acc.size());
    weights.reserve(acc.size());
    for (auto& kv : acc) {
        node_t a = node_t(kv.first >> 32), b = node_t(kv.first & 0xffffffffu);
        edges.emplace_back(a, b);
        weights.push_back(kv.second / 2.0);
    }

    return build_csr(edges, weights, num_comm);
}

// Renumber arbitrary community ids into a dense 0..k-1 range, in place.
// Returns k.
inline node_t renumber(std::vector<node_t>& comm) {
    std::unordered_map<node_t, node_t> remap;
    for (auto& c : comm) {
        auto it = remap.find(c);
        if (it == remap.end()) {
            node_t id = node_t(remap.size());
            remap.emplace(c, id);
            c = id;
        } else {
            c = it->second;
        }
    }
    return node_t(remap.size());
}

}  // namespace detail

// Run Louvain to convergence. resolution > 1 favors more, smaller
// communities; < 1 favors fewer, larger ones. 1.0 is standard modularity.
// Returns comm[v] = final community id for each original node, dense 0..k-1.
//
// If hierarchy_out is non-null, it receives membership[v] snapshotted
// after every aggregation round — hierarchy_out->back() equals the
// returned vector. That lets a caller ask "what was v's community if we
// only zoom out N rounds" instead of only the fully-collapsed answer.
// See query.hpp.
inline std::vector<node_t> louvain(const CSR& g, double resolution = 1.0,
                                    std::vector<std::vector<node_t>>* hierarchy_out = nullptr) {
    if (hierarchy_out) hierarchy_out->clear();

    // membership[v] = which current-level node original node v belongs to,
    // composed across levels as the graph shrinks.
    std::vector<node_t> membership(g.n);
    std::iota(membership.begin(), membership.end(), 0);

    CSR level_graph = g;
    std::vector<node_t> level_comm(g.n);
    std::iota(level_comm.begin(), level_comm.end(), 0);

    while (true) {
        bool moved = detail::local_moving_pass(level_graph, level_comm, resolution);
        node_t num_comm = detail::renumber(level_comm);

        for (size_t v = 0; v < membership.size(); ++v) {
            membership[v] = level_comm[membership[v]];
        }
        if (hierarchy_out) hierarchy_out->push_back(membership);

        // Stop once a pass made no move, or the graph has fully collapsed
        // to one node per community and aggregating again would be a no-op.
        if (!moved || num_comm == level_graph.n) break;

        level_graph = detail::aggregate(level_graph, level_comm, num_comm);
        level_comm.assign(num_comm, 0);
        std::iota(level_comm.begin(), level_comm.end(), 0);
    }

    return membership;
}

}  // namespace gcd
