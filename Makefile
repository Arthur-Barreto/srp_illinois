all: bls gpu_bls

bls: src/bls.cpp
	g++ -O3 src/bls.cpp -o build/bls 

gpu_bls: src/bls.cu
	nvcc -arch=sm_75 -O3 src/bls.cu -o build/gpu_bls 
