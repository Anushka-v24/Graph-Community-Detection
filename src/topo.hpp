// topo.hpp — Kahn's algorithm topological sort.
//
// Used to validate the multi-level community hierarchy that louvain() and
// leiden() build internally (original nodes -> level-0 supernodes ->
// level-1 supernodes -> ...) but normally discard, keeping only the final
// flattened membership. Every node in that hierarchy depends on its parent
// not existing yet when it's created, and stops depending on anything once
// it's folded into its parent — that dependency relation is a DAG (in
// fact a forest), and processing it in topological order (children before
// parents) is exactly the ordering an incremental engine (Phase 4) would
// need when only part of the hierarchy has to be recomputed after an edit.
//
// query.hpp uses this to prove a built Hierarchy has no cycles before
// trusting it for lookups — a real bug in the aggregation loop (a level
// pointing back at itself, say) would otherwise show up as a wrong answer
// somewhere downstream instead of a clear error here.

#pragma once

#include "csr.hpp"

#include <queue>
#include <stdexcept>
#include <vector>

namespace gcd {

// adj[u] lists nodes that depend on u — i.e. an edge u -> v means u must
// be processed before v (here: child supernode before the parent it folds
// into). Returns all nodes in a valid topological order. Throws if the
// graph has a cycle, which should be structurally impossible for a
// hierarchy built by repeated aggregation — if this ever throws, the
// hierarchy itself was built wrong.
inline std::vector<node_t> topo_sort(const std::vector<std::vector<node_t>>& adj) {
    node_t n = node_t(adj.size());
    std::vector<int> indeg(n, 0);
    for (node_t u = 0; u < n; ++u)
        for (node_t v : adj[u]) indeg[v]++;

    std::queue<node_t> q;
    for (node_t u = 0; u < n; ++u)
        if (indeg[u] == 0) q.push(u);

    std::vector<node_t> order;
    order.reserve(n);
    while (!q.empty()) {
        node_t u = q.front();
        q.pop();
        order.push_back(u);
        for (node_t v : adj[u])
            if (--indeg[v] == 0) q.push(v);
    }

    if (order.size() != n) throw std::runtime_error("cycle detected in hierarchy graph");
    return order;
}

}  // namespace gcd
