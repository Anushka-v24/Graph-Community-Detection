// query.hpp — a query layer over a finished community-detection result.
//
// Three kinds of query, three different tools, each picked for what it's
// actually good at rather than forced together:
//
//   same_community(u, v)     flat array compare, O(1) — the partition IS
//                             the answer once detection has run.
//   connected(u, v)          DSU (dsu.hpp) over raw graph edges, ignoring
//                             community labels entirely — "can u reach v
//                             at all". add_edge() previews Phase 4: new
//                             edges update this without a full rerun.
//   community_at_level(v, L) direct index into the aggregation hierarchy
//                             that louvain()/leiden() can optionally record
//                             (see their hierarchy_out parameter), so a
//                             caller can ask for a coarser or finer view
//                             than the fully-converged answer.
//
// The hierarchy is validated once at construction with topo_sort (topo.hpp):
// each level must depend only on the level before it, which is exactly a
// topological order over "level i communities" -> "level i+1 communities"
// they get folded into. That's a real DAG check, not decoration — if a
// future change to the aggregation loop ever produced a level that folded
// back into an earlier one, this catches it as a thrown exception right
// here instead of a silently wrong query answer downstream.

#pragma once

#include "csr.hpp"
#include "dsu.hpp"
#include "topo.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace gcd {

// One DAG node per (level, community) pair. level_offset_out[lvl] is the
// first global id belonging to level lvl, so community c at level lvl has
// global id level_offset_out[lvl] + c. Assumes every hierarchy[lvl] is
// already densely renumbered 0..k-1, which is what louvain()/leiden()
// produce.
inline std::vector<std::vector<node_t>> build_hierarchy_dag(
        const std::vector<std::vector<node_t>>& hierarchy,
        std::vector<node_t>& level_offset_out) {
    size_t num_levels = hierarchy.size();
    level_offset_out.assign(num_levels + 1, 0);
    for (size_t lvl = 0; lvl < num_levels; ++lvl) {
        node_t num_comm = 0;
        for (node_t c : hierarchy[lvl]) num_comm = std::max(num_comm, node_t(c + 1));
        level_offset_out[lvl + 1] = level_offset_out[lvl] + num_comm;
    }

    std::vector<std::vector<node_t>> adj(level_offset_out.back());
    for (size_t lvl = 0; lvl + 1 < num_levels; ++lvl) {
        std::vector<std::unordered_set<node_t>> seen(level_offset_out[lvl + 1] - level_offset_out[lvl]);
        for (size_t v = 0; v < hierarchy[lvl].size(); ++v) {
            node_t local_from = hierarchy[lvl][v];
            node_t from = level_offset_out[lvl] + local_from;
            node_t to = level_offset_out[lvl + 1] + hierarchy[lvl + 1][v];
            if (seen[local_from].insert(to).second) adj[from].push_back(to);
        }
    }
    return adj;
}

class CommunityIndex {
public:
    // Takes the graph the partition was computed on (used once, to seed
    // the connectivity DSU — not retained), the final partition, and
    // optionally the hierarchy that produced it.
    CommunityIndex(const CSR& g, std::vector<node_t> comm,
                   std::vector<std::vector<node_t>> hierarchy = {})
        : comm_(std::move(comm)), hierarchy_(std::move(hierarchy)), conn_(g.n) {
        for (node_t u = 0; u < g.n; ++u)
            for (edge_t i = g.offsets[u]; i < g.offsets[u + 1]; ++i)
                conn_.unite(u, g.neighbors[i]);

        if (!hierarchy_.empty()) {
            std::vector<node_t> level_offset;
            auto dag = build_hierarchy_dag(hierarchy_, level_offset);
            topo_sort(dag);  // throws if the hierarchy is somehow cyclic
        }
    }

    bool same_community(node_t u, node_t v) const { return comm_[u] == comm_[v]; }
    node_t community_of(node_t v) const { return comm_[v]; }

    // Raw reachability, ignoring community boundaries entirely.
    bool connected(node_t u, node_t v) { return conn_.connected(u, v); }

    // Phase-4 preview: register a new edge without a full recompute.
    // Keeps connected() correct; does NOT update same_community() /
    // community_of() for it — making detection itself incremental is the
    // actual Phase 4 work, this only keeps raw connectivity live.
    void add_edge(node_t u, node_t v) { conn_.unite(u, v); }

    // v's community if detection had stopped after `level` aggregation
    // rounds (0 = first round) instead of running to convergence.
    node_t community_at_level(node_t v, size_t level) const {
        if (level >= hierarchy_.size()) throw std::out_of_range("no such hierarchy level");
        return hierarchy_[level][v];
    }

    size_t num_levels() const { return hierarchy_.size(); }

private:
    std::vector<node_t> comm_;
    std::vector<std::vector<node_t>> hierarchy_;
    DSU conn_;
};

}  // namespace gcd
