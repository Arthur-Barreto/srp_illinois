# Compiler and flags
GPP = g++
MPICPP = mpic++
NVCC = nvcc
CXXFLAGS = -O3 -fopenmp
NVCCFLAGS = -arch=sm_86 -O3
BUILDDIR = build
UTILS_SRC = src/utils.cpp
UTILS_OBJ = $(BUILDDIR)/utils
BLS_SRC = src/bls.cpp
BLS_MPI_SRC = src/bls_mpi.cpp
BLS_CU_SRC = src/bls.cu
BLS_BIN = $(BUILDDIR)/bls
BLS_MPI_BIN = $(BUILDDIR)/bls_mpi
GPU_BLS_BIN = $(BUILDDIR)/gpu_bls

# Default target
all: $(BLS_BIN) $(BLS_MPI_BIN) $(GPU_BLS_BIN)

# Create build directory if it doesn't exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile utils
$(UTILS_OBJ): $(UTILS_SRC) | $(BUILDDIR)
	$(GPP) $(CXXFLAGS) -c $(UTILS_SRC) -o $(UTILS_OBJ)

# Build bls
$(BLS_BIN): $(UTILS_OBJ) $(BLS_SRC) | $(BUILDDIR)
	$(GPP) $(CXXFLAGS) $(UTILS_OBJ) $(BLS_SRC) -o $(BLS_BIN)

# Build bls_mpi
$(BLS_MPI_BIN): $(UTILS_OBJ) $(BLS_MPI_SRC) | $(BUILDDIR)
	$(MPICPP) $(CXXFLAGS) $(UTILS_OBJ) $(BLS_MPI_SRC) -o $(BLS_MPI_BIN)

# Build gpu_bls
$(GPU_BLS_BIN): $(BLS_CU_SRC) | $(BUILDDIR)
	$(NVCC) $(NVCCFLAGS) $(BLS_CU_SRC) -o $(GPU_BLS_BIN)

# Clean up build directory
clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean