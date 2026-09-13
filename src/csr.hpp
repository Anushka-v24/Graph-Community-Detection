// csr.hpp — Compressed Sparse Row adjacency structure.
//
// Layout:
//   offsets[v]   = index into neighbors[] where v's adjacency list starts
//   offsets[v+1] = one past the end
//   neighbors[]  = flat array of all adjacency lists, back to back
//
// Why this and not vector<vector<uint32_t>>: one allocation instead of n,
// and neighbours of a node sit contiguously in memory, so walking them is
// a linear scan the CPU prefetcher can predict. At 10^8 edges that
// difference is roughly an order of magnitude.

#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace gcd {

using node_t = uint32_t;  // up to ~4.2B nodes
using edge_t = uint64_t;  // edge counts exceed 32 bits well before node counts

struct CSR {
    node_t n = 0;                    // number of nodes
    edge_t m = 0;                    // number of undirected edges
    std::vector<edge_t> offsets;     // size n+1
    std::vector<node_t> neighbors;   // size 2*m (each edge stored both ways)
    std::vector<double> weights;     // size 2*m, parallel to neighbors

    inline edge_t degree(node_t v) const { return offsets[v + 1] - offsets[v]; }

    inline const node_t* nbr_begin(node_t v) const { return neighbors.data() + offsets[v]; }
    inline const node_t* nbr_end(node_t v) const { return neighbors.data() + offsets[v + 1]; }

    // Sum of incident edge weights. Louvain needs this constantly, so it is
    // worth caching once built rather than recomputing per pass.
    std::vector<double> weighted_degree;

    double total_weight = 0.0;  // sum of all edge weights (counting each edge once)

    void compute_degrees() {
        weighted_degree.assign(n, 0.0);
        double sum = 0.0;
        for (node_t v = 0; v < n; ++v) {
            double d = 0.0;
            for (edge_t i = offsets[v]; i < offsets[v + 1]; ++i) d += weights[i];
            weighted_degree[v] = d;
            sum += d;
        }
        total_weight = sum / 2.0;
    }
};

// Build CSR from an edge list using a two-pass counting sort.
// Pass 1 counts degrees, prefix-sums them into offsets.
// Pass 2 places each endpoint. O(n + m), no per-node allocation.
//
// Treats the graph as undirected: every (u,v) is stored at u and at v.
inline CSR build_csr(const std::vector<std::pair<node_t, node_t>>& edges,
                     const std::vector<double>& edge_weights,
                     node_t num_nodes) {
    CSR g;
    g.n = num_nodes;
    g.m = edges.size();

    const bool weighted = !edge_weights.empty();
    if (weighted && edge_weights.size() != edges.size())
        throw std::runtime_error("weights size does not match edges size");

    // Pass 1: degree count
    g.offsets.assign(g.n + 1, 0);
    for (const auto& e : edges) {
        if (e.first >= g.n || e.second >= g.n)
            throw std::runtime_error("edge endpoint out of range");
        g.offsets[e.first + 1]++;
        g.offsets[e.second + 1]++;
    }

    // Prefix sum turns counts into start positions
    for (node_t v = 0; v < g.n; ++v) g.offsets[v + 1] += g.offsets[v];

    // Pass 2: scatter
    g.neighbors.resize(2 * g.m);
    g.weights.assign(2 * g.m, 1.0);
    std::vector<edge_t> cursor(g.offsets.begin(), g.offsets.end() - 1);

    for (size_t i = 0; i < edges.size(); ++i) {
        node_t u = edges[i].first, v = edges[i].second;
        double w = weighted ? edge_weights[i] : 1.0;
        g.neighbors[cursor[u]] = v;  g.weights[cursor[u]] = w;  cursor[u]++;
        g.neighbors[cursor[v]] = u;  g.weights[cursor[v]] = w;  cursor[v]++;
    }

    g.compute_degrees();
    return g;
}

// Read a whitespace-separated edge list: "u v" or "u v w" per line.
// Lines starting with # or % are skipped (SNAP and Matrix Market both use these).
// Node ids are assumed to be 0-based and dense. If your file is 1-based or
// sparse, remap first — a dense id space is what makes CSR indexing free.
inline CSR load_edge_list(const std::string& path, bool has_weights = false) {
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) throw std::runtime_error("cannot open " + path);

    std::vector<std::pair<node_t, node_t>> edges;
    std::vector<double> weights;
    edges.reserve(1 << 20);

    node_t max_id = 0;
    char line[512];

    while (std::fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '%' || line[0] == '\n') continue;

        unsigned long u, v;
        double w = 1.0;
        int got = has_weights ? std::sscanf(line, "%lu %lu %lf", &u, &v, &w)
                              : std::sscanf(line, "%lu %lu", &u, &v);
        if (got < 2) continue;
        if (u == v) continue;  // drop self-loops; they do not affect modularity

        edges.emplace_back(static_cast<node_t>(u), static_cast<node_t>(v));
        if (has_weights) weights.push_back(w);
        max_id = std::max<node_t>(max_id, std::max<node_t>(u, v));
    }
    std::fclose(f);

    return build_csr(edges, weights, max_id + 1);
}

}  // namespace gcd
