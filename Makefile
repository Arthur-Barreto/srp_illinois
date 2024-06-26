# Compiler and flags
CXX = g++
NVCC = nvcc
CXXFLAGS = -O3
NVCCFLAGS = -arch=sm_75 -O3

# Source files
SRC_CPP = src/bls.cpp
SRC_CU = src/bls.cu

# Output files
OUT_CPP = bls
OUT_CU = gpu_bls

# Default target
all: $(OUT_CPP) $(OUT_CU)

# Compile C++ source
$(OUT_CPP): $(SRC_CPP)
	$(CXX) $(CXXFLAGS) $< -o $@

# Compile CUDA source
$(OUT_CU): $(SRC_CU)
	$(NVCC) $(NVCCFLAGS) $< -o $@

# Clean target
clean:
	rm -f $(OUT_CPP) $(OUT_CU)

.PHONY: all clean
