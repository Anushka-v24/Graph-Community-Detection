// connectivity.hpp — guarantee every reported community is actually one
// connected piece of the graph, regardless of which algorithm produced it.
//
// Louvain has no such guarantee at all. Leiden's refinement phase is
// *supposed* to guarantee it via a well-connectedness threshold — but that
// threshold is easy to get subtly wrong (this project's own leiden.hpp has
// done it twice: a check of the shape `x < threshold && x == 0.0` always
// collapses to just `x == 0.0`, silently deleting the threshold entirely).
//
// Rather than trust either algorithm's internal bookkeeping, this runs an
// explicit DFS over each community's induced subgraph after detection and
// splits any community that turns out to be more than one connected piece.
// That makes "communities are connected" a property you can prove by
// construction, not just hope the threshold math got right.

#pragma once

#include "csr.hpp"

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gcd {

// Returns a new community assignment where every community is a single
// connected component. Communities that were already connected keep their
// members (just possibly renumbered); disconnected ones are split, each
// piece getting a fresh id. If splits_out is given, it's set to the number
// of original communities that were not already connected (0 = input was
// already clean).
inline std::vector<node_t> enforce_connectivity(const CSR& g,
                                                 const std::vector<node_t>& comm,
                                                 size_t* splits_out = nullptr) {
    const node_t UNSET = node_t(-1);
    std::vector<node_t> out(g.n, UNSET);
    std::vector<node_t> stack;
    node_t next_id = 0;

    for (node_t start = 0; start < g.n; ++start) {
        if (out[start] != UNSET) continue;

        // DFS (explicit stack, so depth can't blow the call stack on a
        // large connected piece) over every unvisited node reachable from
        // `start` while staying inside its declared community.
        node_t original_comm = comm[start];
        node_t this_piece = next_id++;
        out[start] = this_piece;
        stack.push_back(start);

        while (!stack.empty()) {
            node_t u = stack.back();
            stack.pop_back();
            for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i) {
                node_t v = g.neighbors[i];
                if (v != u && comm[v] == original_comm && out[v] == UNSET) {
                    out[v] = this_piece;
                    stack.push_back(v);
                }
            }
        }
    }

    if (splits_out) {
        std::unordered_map<node_t, std::unordered_set<node_t>> pieces_per_comm;
        for (node_t v = 0; v < g.n; ++v) pieces_per_comm[comm[v]].insert(out[v]);
        size_t splits = 0;
        for (auto& kv : pieces_per_comm) if (kv.second.size() > 1) ++splits;
        *splits_out = splits;
    }

    return out;
}

}  // namespace gcd
