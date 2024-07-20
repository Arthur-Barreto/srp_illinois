# Compiler and flags
CXX = g++
NVCC = nvcc
CXXFLAGS = -O3 -fopenmp
NVCCFLAGS = -arch=sm_75 -O3

# Source files
SRC_UTILS = src/utils.cpp
SRC_CPP = src/bls.cpp
SRC_CU = src/bls.cu

# Output files
OUT_UTILS = build/utils
OUT_CPP = build/bls
OUT_CU = build/gpu_bls

# Default target
all: $(OUT_UTILS) $(OUT_CPP) $(OUT_CU) 

# COMpile uitls
$(OUT_UTILS) : $(SRC_UTILS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile C++ source
$(OUT_CPP): $(SRC_CPP)
	$(CXX) $(CXXFLAGS) $(OUT_UTILS) $< -o $@

# Compile CUDA source
$(OUT_CU): $(SRC_CU)
	$(NVCC) $(NVCCFLAGS) $< -o $@

# Clean target
clean:
	rm -f $(OUT_CPP) $(OUT_CU) $(OUT_UTILS)

.PHONY: all clean
