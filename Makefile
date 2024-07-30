# Compiler and flags
GPP = g++
MPICPP = mpic++
CXXFLAGS = -O3 -Wall -Wextra -g -fsanitize=address
OPENMP_FLAGS = -fopenmp
MPI_LIB = -lmpi
BUILDDIR = build
UTILS_SRC = src/utils.cpp
UTILS_OMP_SRC = src/utils_omp.cpp
UTILS_MPI_SRC = src/utils_mpi.cpp
UTILS_OBJ = $(BUILDDIR)/utils
UTILS_OMP_OBJ = $(BUILDDIR)/utils_omp
UTILS_MPI_OBJ = $(BUILDDIR)/utils_mpi
BLS_SRC = src/bls.cpp
BLS_OMP_SRC = src/bls_omp.cpp
BLS_MPI_SRC = src/bls_mpi.cpp
BLS_BIN = $(BUILDDIR)/bls
BLS_OMP_BIN = $(BUILDDIR)/bls_omp
BLS_MPI_BIN = $(BUILDDIR)/bls_mpi

# Default target
all: $(BLS_BIN) $(BLS_OMP_BIN) $(BLS_MPI_BIN)

# Create build directory if it doesn't exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile utils
$(UTILS_OBJ): $(UTILS_SRC) | $(BUILDDIR)
	$(GPP) $(CXXFLAGS) -c $(UTILS_SRC) -o $(UTILS_OBJ)

# Compile utils_omp
$(UTILS_OMP_OBJ): $(UTILS_OMP_SRC) | $(BUILDDIR)
	$(GPP) $(CXXFLAGS) $(OPENMP_FLAGS) -c $(UTILS_OMP_SRC) -o $(UTILS_OMP_OBJ)

# Compile utils_mpi
$(UTILS_MPI_OBJ): $(UTILS_MPI_SRC) | $(BUILDDIR)
	$(MPICPP) $(CXXFLAGS) -c $(UTILS_MPI_SRC) -o $(UTILS_MPI_OBJ)

# Build bls
$(BLS_BIN): $(UTILS_OBJ) $(BLS_SRC) | $(BUILDDIR)
	$(GPP) $(CXXFLAGS) $(UTILS_OBJ) $(BLS_SRC) -o $(BLS_BIN)

# Build bls_omp
$(BLS_OMP_BIN): $(UTILS_OMP_OBJ) $(BLS_OMP_SRC) | $(BUILDDIR)
	$(GPP) $(CXXFLAGS) $(OPENMP_FLAGS) $(UTILS_OMP_OBJ) $(BLS_OMP_SRC) -o $(BLS_OMP_BIN)

# Build bls_mpi
$(BLS_MPI_BIN): $(UTILS_MPI_OBJ) $(BLS_MPI_SRC) | $(BUILDDIR)
	$(MPICPP) $(CXXFLAGS) $(UTILS_MPI_OBJ) $(BLS_MPI_SRC) -o $(BLS_MPI_BIN)

# Clean up build directory
clean:
	rm -rf $(BUILDDIR)

.PHONY: all clean
