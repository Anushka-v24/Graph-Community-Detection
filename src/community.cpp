// community.cpp — run Louvain or Leiden on a graph, report modularity, and (given a
// ground-truth file) NMI. This is the Phase 2 harness: every run's numbers
// get appended to bench/results.csv so runs are comparable over time.
//
// Build: g++ -O2 -std=c++17 -o bin/community src/community.cpp
// Run:   ./bin/community data/lfr_1k_mu04.edges data/lfr_1k_mu04.truth --algo leiden

#include "csr.hpp"
#include "louvain.hpp"
#include "leiden.hpp"
#include "metrics.hpp"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
using namespace gcd;

// Truth file: "node community_id" per line, as written by tools/gen_lfr.py.
static std::vector<node_t> load_truth(const std::string& path, node_t n) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);

    std::vector<node_t> truth(n, 0);
    unsigned long node, c;
    while (f >> node >> c) {
        if (node < n) truth[node] = static_cast<node_t>(c);
    }
    return truth;
}

static void append_bench_row(const std::string& algo, const std::string& edge_path, const CSR& g, node_t num_comm,
                              double modularity_score, double nmi_score, double ms) {
    bool need_header = std::ifstream("bench/results.csv").peek() == std::ifstream::traits_type::eof();
    std::ofstream out("bench/results.csv", std::ios::app);
    if (!out) return;  // bench/ missing or unwritable — don't fail the run over it

    if (need_header) out << "algo,graph,nodes,edges,communities,modularity,nmi,ms\n";
    out << algo << ',' << edge_path << ',' << g.n << ',' << g.m << ',' << num_comm << ','
        << modularity_score << ',' << (nmi_score >= 0.0 ? std::to_string(nmi_score) : "")
        << ',' << ms << '\n';
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <edgelist> [truth_file] [--algo louvain|leiden] [--weighted]\n", argv[0]);
        return 1;
    }

    std::string edge_path = argv[1];
    std::string truth_path;
    std::string algo = "louvain";
    bool weighted = false;

    for (int i = 2; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--weighted") {
            weighted = true;
        } else if (a == "--algo" && i + 1 < argc) {
            algo = argv[++i];
        } else if (a == "--leiden") {
            algo = "leiden";
        } else if (a == "--louvain") {
            algo = "louvain";
        } else if (a.rfind("--", 0) != 0 && truth_path.empty()) {
            truth_path = a;
        }
    }

    CSR g = load_edge_list(edge_path, weighted);

    auto t0 = std::chrono::steady_clock::now();
    std::vector<node_t> comm;
    if (algo == "leiden") {
        comm = leiden(g);
    } else {
        comm = louvain(g);
        algo = "louvain";
    }
    auto t1 = std::chrono::steady_clock::now();
    double runtime_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    node_t num_comm = 0;
    for (node_t c : comm) num_comm = std::max(num_comm, node_t(c + 1));
    double q = modularity(g, comm);

    std::printf("algorithm        %s\n", algo.c_str());
    std::printf("nodes            %u\n", g.n);
    std::printf("edges            %llu\n", (unsigned long long)g.m);
    std::printf("communities      %u\n", num_comm);
    std::printf("modularity       %.4f\n", q);
    std::printf("runtime          %.1f ms\n", runtime_ms);

    double score_nmi = -1.0;
    if (!truth_path.empty()) {
        std::vector<node_t> truth = load_truth(truth_path, g.n);
        score_nmi = nmi(comm, truth);
        std::printf("NMI vs truth     %.4f\n", score_nmi);
    }

    append_bench_row(algo, edge_path, g, num_comm, q, score_nmi, runtime_ms);
    return 0;
}
