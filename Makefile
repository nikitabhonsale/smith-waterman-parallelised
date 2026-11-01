# Makefile - Wavefront Smith-Waterman with OpenMP
#
#   macOS  :  brew install libomp      then   make
#             (Apple's clang has no built-in OpenMP, so this file points it at Homebrew's libomp.
#              If you prefer GCC instead:  brew install gcc  then  make CXX=g++-14  -- use your version)
#   Linux  :  make                     (g++ with OpenMP)
#
#   make test    build and run the correctness tests
#   make demo    align the small example sequences
#   make bench   run the benchmark on the large example sequences (writes bench_results.csv)
#   make clean

UNAME := $(shell uname -s)

CXXFLAGS := -O3 -std=c++17 -Wall -Wextra
OMP_CXXFLAGS := -fopenmp
OMP_LDFLAGS :=

ifeq ($(UNAME),Darwin)
  ifeq ($(origin CXX),default)
    CXX := clang++
  endif
  ifneq ($(findstring clang,$(CXX)),)
    LIBOMP := $(shell brew --prefix libomp 2>/dev/null)
    ifeq ($(LIBOMP),)
      $(warning libomp was not found. Install it with: brew install libomp)
    endif
    OMP_CXXFLAGS := -Xpreprocessor -fopenmp -I$(LIBOMP)/include
    OMP_LDFLAGS  := -L$(LIBOMP)/lib -lomp -Wl,-rpath,$(LIBOMP)/lib
  endif
endif

SRC_DIR := src
BUILD   := build
BIN     := bin

# Code shared by all three programs
LIB_OBJS := $(BUILD)/fasta.o $(BUILD)/sw_sequential.o $(BUILD)/sw_wavefront.o \
            $(BUILD)/sw_tiled.o $(BUILD)/alignment.o

HEADERS := $(wildcard $(SRC_DIR)/*.h)

.PHONY: all test demo bench clean

all: $(BIN)/sw_align $(BIN)/sw_bench $(BIN)/sw_tests

$(BIN)/sw_align: $(BUILD)/main.o $(LIB_OBJS) | $(BIN)
	$(CXX) $^ -o $@ $(OMP_CXXFLAGS) $(OMP_LDFLAGS)

$(BIN)/sw_bench: $(BUILD)/bench.o $(LIB_OBJS) | $(BIN)
	$(CXX) $^ -o $@ $(OMP_CXXFLAGS) $(OMP_LDFLAGS)

$(BIN)/sw_tests: $(BUILD)/tests.o $(LIB_OBJS) | $(BIN)
	$(CXX) $^ -o $@ $(OMP_CXXFLAGS) $(OMP_LDFLAGS)

$(BUILD)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(OMP_CXXFLAGS) -c $< -o $@

$(BUILD) $(BIN):
	mkdir -p $@

test: $(BIN)/sw_tests
	./$(BIN)/sw_tests

demo: $(BIN)/sw_align
	./$(BIN)/sw_align examples/small_a.fasta examples/small_b.fasta --out alignment_result.txt

bench: $(BIN)/sw_bench
	./$(BIN)/sw_bench examples/large_a.fasta examples/large_b.fasta --csv bench_results.csv

clean:
	rm -rf $(BUILD) $(BIN)
