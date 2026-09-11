CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Isrc
DBGFLAGS ?= -O0 -g -std=c++17 -Wall -Wextra -Isrc
LDFLAGS  ?=

# Uncomment once you add OpenMP to the local-moving loop in Phase 4.
# CXXFLAGS += -fopenmp
# LDFLAGS  += -fopenmp

BIN := bin

.PHONY: all debug clean

all: $(BIN)/stats $(BIN)/community

debug: $(BIN)/stats_debug $(BIN)/community_debug

$(BIN)/stats: src/stats.cpp src/csr.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ src/stats.cpp $(LDFLAGS)

$(BIN)/community: src/community.cpp src/csr.hpp src/louvain.hpp src/leiden.hpp src/metrics.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -o $@ src/community.cpp $(LDFLAGS)

# Address + UB sanitizers on. Slow, but they catch the kind of off-by-one in
# the CSR scatter loop that would otherwise surface as silently wrong
# communities three phases from now.
$(BIN)/stats_debug: src/stats.cpp src/csr.hpp | $(BIN)
	$(CXX) $(DBGFLAGS) -fsanitize=address,undefined -o $@ src/stats.cpp -fsanitize=address,undefined

$(BIN)/community_debug: src/community.cpp src/csr.hpp src/louvain.hpp src/leiden.hpp src/metrics.hpp | $(BIN)
	$(CXX) $(DBGFLAGS) -fsanitize=address,undefined -o $@ src/community.cpp -fsanitize=address,undefined

$(BIN):
	mkdir -p $(BIN)

clean:
	rm -rf $(BIN)
