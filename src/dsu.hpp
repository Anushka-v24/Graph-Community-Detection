// dsu.hpp — Union-Find (Disjoint Set Union) with path compression and
// union by size.
//
// This is the standard structure for "are these two things in the same
// group" queries that need to keep working cheaply as groups keep merging
// — which is exactly the shape of a graph that gets new edges over time.
// Two operations, both close to O(1) amortized:
//   find(x)   — which group is x currently in (returns a representative id)
//   unite(a,b) — merge a's group and b's group into one
//
// Used two ways in this project: as the backbone of the incremental query
// engine (query.hpp), and inside connectivity.hpp's connected-components
// check as an alternative way to think about the same DFS result.

#pragma once

#include "csr.hpp"

#include <numeric>
#include <vector>

namespace gcd {

class DSU {
public:
    explicit DSU(node_t n) : parent_(n), size_(n, 1) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    node_t find(node_t x) {
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];  // path halving
            x = parent_[x];
        }
        return x;
    }

    // Returns true if a and b were in different groups (and are now merged).
    bool unite(node_t a, node_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (size_[a] < size_[b]) std::swap(a, b);
        parent_[b] = a;
        size_[a] += size_[b];
        return true;
    }

    bool connected(node_t a, node_t b) { return find(a) == find(b); }
    node_t component_size(node_t x) { return size_[find(x)]; }
    node_t n() const { return node_t(parent_.size()); }

private:
    std::vector<node_t> parent_;
    std::vector<node_t> size_;
};

}  // namespace gcd
