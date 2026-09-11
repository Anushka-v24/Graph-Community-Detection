#!/usr/bin/env python3
"""
gen_lfr.py -- generate LFR benchmark graphs with known ground-truth communities.

LFR graphs are the standard test bed for community detection. They are built
so that each node has a community assigned in advance, with a tunable mixing
parameter mu: the fraction of a node's edges that go OUTSIDE its community.

  mu = 0.1  -> very clear communities, any algorithm finds them
  mu = 0.5  -> hard, where good algorithms separate from bad ones
  mu = 0.7  -> most methods break down

Because you know the true answer, you can measure NMI (how close your output
is to truth). Without this you can only measure modularity, which tells you
your partition scores well by its own objective -- not that it is correct.
That distinction is the whole reason this file exists before any algorithm.

Requires: pip install networkx

Usage:
  python tools/gen_lfr.py --n 1000 --mu 0.1 --out data/lfr_1k
  python tools/gen_lfr.py --n 100000 --mu 0.4 --out data/lfr_100k

Writes:
  <out>.edges       "u v" per line, 0-based dense node ids
  <out>.truth       "node community_id" per line
"""

import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=1000, help="number of nodes")
    p.add_argument("--mu", type=float, default=0.1,
                   help="mixing parameter, 0..1 (higher = harder)")
    p.add_argument("--tau1", type=float, default=2.5,
                   help="power-law exponent for degree distribution")
    p.add_argument("--tau2", type=float, default=1.5,
                   help="power-law exponent for community size distribution")
    p.add_argument("--avg-degree", type=float, default=15.0)
    p.add_argument("--min-community", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=str, required=True, help="output path prefix")
    args = p.parse_args()

    try:
        import networkx as nx
    except ImportError:
        sys.exit("networkx not installed. Run: pip install networkx")

    print(f"generating n={args.n} mu={args.mu} avg_deg={args.avg_degree} ...")

    G = None
    max_attempts = 50
    for attempt in range(max_attempts):
        current_seed = args.seed + attempt if args.seed is not None else None
        try:
            G = nx.LFR_benchmark_graph(
                n=args.n,
                tau1=args.tau1,
                tau2=args.tau2,
                mu=args.mu,
                average_degree=args.avg_degree,
                min_community=args.min_community,
                max_iters=2000,
                seed=current_seed,
            )
            break
        except Exception as e:
            if attempt == max_attempts - 1:
                sys.exit(f"Failed to generate LFR graph after {max_attempts} attempts: {e}")

    # LFR attaches each node's community as a set on the node itself.
    communities = {}
    comm_id = {}
    for v in G.nodes():
        members = frozenset(G.nodes[v]["community"])
        if members not in comm_id:
            comm_id[members] = len(comm_id)
        communities[v] = comm_id[members]

    G.remove_edges_from(nx.selfloop_edges(G))

    edge_path = args.out + ".edges"
    truth_path = args.out + ".truth"

    with open(edge_path, "w") as f:
        for u, v in G.edges():
            f.write(f"{u} {v}\n")

    with open(truth_path, "w") as f:
        for v in sorted(communities):
            f.write(f"{v} {communities[v]}\n")

    sizes = {}
    for c in communities.values():
        sizes[c] = sizes.get(c, 0) + 1

    print(f"  nodes:       {G.number_of_nodes()}")
    print(f"  edges:       {G.number_of_edges()}")
    print(f"  communities: {len(sizes)}")
    print(f"  sizes:       min={min(sizes.values())} max={max(sizes.values())}")
    print(f"  wrote {edge_path} and {truth_path}")


if __name__ == "__main__":
    main()
