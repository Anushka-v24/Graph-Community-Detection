// stats.cpp — load a graph, print what it looks like.
//
// This is deliberately the first thing you build. Before any algorithm,
// you want to be certain the loader is correct and to have a feel for the
// degree distribution — it drives every performance decision later.
//
// Build: g++ -O2 -std=c++17 -o bin/stats src/stats.cpp
// Run:   ./bin/stats data/lfr_1k.edges

#include "csr.hpp"

#include <chrono>
#include <cstdio>
#include <numeric>

using namespace std;
using namespace gcd;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <edgelist> [--weighted]\n", argv[0]);
        return 1;
    }
    bool weighted = (argc > 2 && std::string(argv[2]) == "--weighted");

    auto t0 = std::chrono::steady_clock::now();
    CSR g = load_edge_list(argv[1], weighted);
    auto t1 = std::chrono::steady_clock::now();
    double load_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Degree distribution
    edge_t max_deg = 0;
    edge_t zero_deg = 0;
    for (node_t v = 0; v < g.n; ++v) {
        edge_t d = g.degree(v);
        max_deg = std::max(max_deg, d);
        if (d == 0) zero_deg++;
    }
    double avg_deg = g.n ? (2.0 * g.m) / g.n : 0.0;

    // Memory: this is the number to watch as you scale up.
    double bytes = static_cast<double>(g.offsets.size() * sizeof(edge_t)
                                     + g.neighbors.size() * sizeof(node_t)
                                     + g.weights.size() * sizeof(double)
                                     + g.weighted_degree.size() * sizeof(double));

    std::printf("nodes            %u\n", g.n);
    std::printf("edges            %llu\n", (unsigned long long)g.m);
    std::printf("avg degree       %.2f\n", avg_deg);
    std::printf("max degree       %llu\n", (unsigned long long)max_deg);
    std::printf("isolated nodes   %llu\n", (unsigned long long)zero_deg);
    std::printf("total weight     %.1f\n", g.total_weight);
    std::printf("load time        %.1f ms\n", load_ms);
    std::printf("memory           %.1f MB  (%.1f bytes/edge)\n",
                bytes / 1e6, g.m ? bytes / (2.0 * g.m) : 0.0);

    // Log-scale degree histogram — real graphs are heavy-tailed and you
    // want to see that before you start reasoning about load balance.
    std::printf("\ndegree histogram (log2 buckets)\n");
    std::vector<edge_t> hist(34, 0);
    for (node_t v = 0; v < g.n; ++v) {
        edge_t d = g.degree(v);
        int b = 0;
        while ((edge_t(1) << (b + 1)) <= d && b < 32) ++b;
        hist[d == 0 ? 0 : b + 1]++;
    }
    for (size_t b = 0; b < hist.size(); ++b) {
        if (!hist[b]) continue;
        if (b == 0) std::printf("  deg 0       : %llu\n", (unsigned long long)hist[0]);
        else std::printf("  deg %6llu+ : %llu\n",
                         (unsigned long long)(edge_t(1) << (b - 1)),
                         (unsigned long long)hist[b]);
    }
    return 0;
}
