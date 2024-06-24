#!/bin/bash
g++ -O3 src/bls.cpp -o bls
nvcc -arch=sm_75 -O3 src/bls.cu -o gpu_bls 