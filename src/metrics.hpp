// metrics.hpp — modularity and NMI, the two scores every community
// detection run is judged by.
//
// Modularity measures how much more densely connected a community is
// than you'd expect by chance, given the degree sequence. It only looks
// at the graph and the partition you handed it — it has no idea what the
// "right" answer is, so a high score does not mean a correct answer.
//
// NMI (Normalized Mutual Information) compares two partitions of the same
// nodes and asks how much knowing one tells you about the other. It needs
// ground truth, which is exactly what LFR graphs provide — see README.

#pragma once

#include "csr.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

using namespace std;

namespace gcd {

// Q = (1/2m) * sum over ordered (u,v) pairs of [ w(u,v) - k_u*k_v/(2m) ]
// restricted to pairs in the same community. Computed via per-community
// running sums so it stays O(n + m) instead of comparing every pair.
inline double modularity(const CSR& g, const std::vector<node_t>& comm) {
    if (g.total_weight <= 0.0 || g.n == 0) return 0.0;
    double two_m = 2.0 * g.total_weight;

    node_t num_comm = 0;
    for (node_t v = 0; v < g.n; ++v) num_comm = std::max(num_comm, node_t(comm[v] + 1));

    std::vector<double> internal(num_comm, 0.0);   // internal edge weight (double-counted), per community
    std::vector<double> tot_degree(num_comm, 0.0); // total weighted degree, per community

    for (node_t u = 0; u < g.n; ++u) {
        tot_degree[comm[u]] += g.weighted_degree[u];
        for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
            node_t v = g.neighbors[i];
            if (comm[v] == comm[u]) internal[comm[u]] += g.weights[i];
        }
    }

    double q = 0.0;
    for (node_t c = 0; c < num_comm; ++c) {
        double frac_tot = tot_degree[c] / two_m;
        q += internal[c] / two_m - frac_tot * frac_tot;
    }
    return q;
}

// NMI between two partitions of the same n nodes, via the standard
// contingency-table formula: I(A,B) / sqrt(H(A) * H(B)).
// 1.0 = identical partitions, 0.0 = independent (or degenerate) ones.
inline double nmi(const std::vector<node_t>& a, const std::vector<node_t>& b) {
    size_t n = a.size();
    if (n == 0 || b.size() != n) return 0.0;

    std::unordered_map<uint64_t, edge_t> joint;
    std::unordered_map<node_t, edge_t> count_a, count_b;

    for (size_t i = 0; i < n; ++i) {
        count_a[a[i]]++;
        count_b[b[i]]++;
        joint[(uint64_t(a[i]) << 32) | b[i]]++;
    }

    double h_a = 0.0, h_b = 0.0, mi = 0.0;
    for (auto& kv : count_a) { double p = double(kv.second) / n; h_a -= p * std::log(p); }
    for (auto& kv : count_b) { double p = double(kv.second) / n; h_b -= p * std::log(p); }

    for (auto& kv : joint) {
        node_t ca = node_t(kv.first >> 32), cb = node_t(kv.first & 0xffffffffu);
        double p_ab = double(kv.second) / n;
        double p_a = double(count_a[ca]) / n;
        double p_b = double(count_b[cb]) / n;
        mi += p_ab * std::log(p_ab / (p_a * p_b));
    }

    // Both sides collapsed to a single community: no information to share,
    // but also nothing to disagree on. Call that trivially identical.
    if (h_a == 0.0 && h_b == 0.0) return 1.0;
    if (h_a == 0.0 || h_b == 0.0) return 0.0;
    return mi / std::sqrt(h_a * h_b);
}

}  // namespace gcd
