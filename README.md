# gcd — graph community detection at scale

An incremental community detection engine for continuously updating graphs.

**Goal:** maintain community structure over a large, constantly changing graph,
and answer filtered community queries in bounded time — without recomputing
from scratch on every update.

Current status: **Phase 2 Complete — static Louvain & Leiden + modularity/NMI harness, a DFS-verified connectivity guarantee, and a first query layer over the result (DSU + topo-sort).** Benchmarking against NetworKit/igraph (Phase 3) and the Incremental Update Engine (Phase 4) are next.

---

## Setup

Requires g++ (or clang++) with C++17, make, and Python 3.

```bash
make                       # builds bin/stats and bin/community
pip install networkx       # for the LFR generator
python3 tools/gen_lfr.py --n 1000 --mu 0.1 --out data/lfr_1k
./bin/stats data/lfr_1k.edges
./bin/community data/lfr_1k.edges data/lfr_1k.truth --algo leiden   # runs Leiden, reports modularity + NMI
```

### VS Code

Open the `gcd/` folder directly (not a parent directory) so `.vscode/` is picked up.

Install the **C/C++** extension (`ms-vscode.cpptools`). Then:

- `Ctrl+Shift+B` — build
- `Ctrl+Shift+P` → *Tasks: Run Task* — LFR generation, stats, clean
- `F5` — debug build with sanitizers, breakpoints in the loader

---

## Layout

```
src/csr.hpp           Compressed Sparse Row adjacency — the core data structure
src/stats.cpp         Loads a graph, prints node/edge/degree/memory stats
src/louvain.hpp       Louvain method: local moving + aggregation
src/leiden.hpp        Leiden method: fast local moving + sub-community refinement + coarsening
src/metrics.hpp       Modularity and NMI scoring
src/dsu.hpp           Union-Find (Disjoint Set) — path compression + union by size
src/connectivity.hpp  DFS pass that guarantees every community is one connected piece
src/topo.hpp          Kahn's-algorithm topological sort, used to validate the hierarchy DAG
src/query.hpp         Query layer: same-community, DSU connectivity, multi-level lookups
src/community.cpp     Runs Louvain/Leiden, enforces connectivity, reports metrics, logs to bench/
tools/gen_lfr.py      LFR benchmark graphs with ground-truth communities
data/                 Graph files (gitignored)
bench/                Benchmark results and plots (results.csv from bin/community)
```

---

## Why LFR graphs

LFR graphs come with the correct community assignment baked in, controlled by a
mixing parameter `mu` — the fraction of each node's edges that leave its own
community.

| mu    | difficulty                                |
|-------|-------------------------------------------|
| 0.1   | easy, any method finds the communities    |
| 0.4   | moderate                                  |
| 0.5+  | hard, where good methods separate from bad|

Because the truth is known, correctness can be measured with NMI rather than
only modularity. Modularity tells you a partition scores well against its own
objective; it does not tell you the partition is right. Every measurement in
this project is validated against ground truth first.

---

## Roadmap

| Phase | Work                                                            |
|-------|-----------------------------------------------------------------|
| 1     | CSR loader, LFR generation, ID remapping, binary format + mmap  |
| 2     | Static Louvain & Leiden, Modularity + NMI harness               |
| 3     | Benchmark against NetworKit and igraph on identical graphs      |
| 4     | **Incremental update engine** — drift curve vs full recompute   |
| 5     | Run-to-run stability of incremental results                     |
| 6     | Time-windowed communities (sliding window)                      |
| 7     | Filtered queries: detect-on-subgraph vs filter-precomputed      |
| 8     | Fraud ring demo on a synthetic transaction graph                |
| 9     | Scale run — 10^8 local, largest affordable cloud run            |
| 10    | Multi-tenancy, UI, design doc                                   |

Phase 4 is the contribution. Everything before it exists to make its numbers
credible; everything after extends it.

---

## Remaining Phase 1 tasks

- [ ] Node ID remapper — real datasets have sparse or 1-based IDs; the loader
      assumes dense 0-based
- [ ] Binary edge list format + `mmap` loader — text parsing dominates above
      roughly 10M edges
- [x] Generate at mu = 0.1 / 0.4 / 0.6 and confirm stats look sane at each —
      done on `data/lfr_1k_mu{01,04,06}`; all load clean (0 isolated nodes,
      consistent bytes/edge)

## Phase 2 tasks (Complete)

- [x] Static Louvain (`src/louvain.hpp`)
- [x] Modularity + NMI harness (`src/metrics.hpp`, `bin/community`)
- [x] Leiden (`src/leiden.hpp`) — fast local moving with queue + sub-community refinement
- [x] Connectivity guarantee (`src/connectivity.hpp`) — Leiden's refinement phase is
      *supposed* to guarantee every community is connected via a well-connectedness
      threshold, but that threshold has twice been implemented as `x < thresh && x == 0.0`,
      which always collapses to just `x == 0.0` and silently drops the threshold
      entirely. Rather than trust that math a third time, every run now runs an
      explicit DFS over each community's induced subgraph afterward and splits
      anything that isn't actually one connected piece — a guarantee proven by
      construction instead of assumed from a formula.
- [x] Query layer (`src/query.hpp`) over a finished result: O(1) same-community
      lookups, a DSU-backed (`src/dsu.hpp`) raw connectivity index that previews
      Phase 4's incremental idea (`add_edge()` updates it without a full rerun),
      and multi-level lookups against the aggregation hierarchy — validated once
      with Kahn's-algorithm topological sort (`src/topo.hpp`) to prove it's
      actually a DAG before trusting it.
- [x] Sanity numbers on `data/lfr_1k_mu{01,04,06}` (post connectivity-fix, 0 splits
      needed on any of them — these graphs are small/easy enough not to trigger the
      failure mode, which is expected; the guarantee matters more at scale):
      - Louvain: NMI 1.00 / 0.74 / 0.16
      - Leiden:  NMI 1.00 / 0.76 / 0.17
      Tracks the expected difficulty table and provides a baseline for Phase 4
      incremental drift curves.

## Metrics tracked

Recorded for every run from Phase 2 onward, into `bench/`:

- modularity
- NMI against ground truth
- wall-clock time
- peak memory, and bytes per edge
- (Phase 4+) update latency p50/p99, drift from full recompute
