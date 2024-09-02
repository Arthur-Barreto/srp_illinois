#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <vector>

using namespace std;

struct BlsResult {
  double period;
  double d_value;
};

__global__ void weight_sum_kernel(const double *flux_err, double *weight,
                                  int n) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < n) {
    atomicAdd(weight, pow(flux_err[idx], -2));
  }
}

__global__ void compute_weight(const double *flux_err, double *weight,
                               double sum_w, int n) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < n) {
    weight[idx] = sum_w * pow(flux_err[idx], -2);
  }
}

__global__ void bls_kernel(const double *time, const double *flux,
                           const double *weight, int n, BlsResult *result) {
  int idx1 = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx1 >= n)
    return;

  for (int idx2 = idx1 + 1; idx2 < n; idx2++) {
    double r = 0, s = 0, aux = 0;
    for (int i = idx1; i <= idx2; i++) {
      r += weight[i];
      s += weight[i] * flux[i];
    }
    for (int i = 0; i < n; i++) {
      aux += weight[i] * pow(flux[i], 2);
    }
    double d = aux - pow(s, 2) / (r * (1 - r));
    double period = time[idx2] - time[idx1];
    if (d < result->d_value) {
      result->d_value = d;
      result->period = period;
    }
  }
}

void readCSV(const string &filename, vector<double> &time, vector<double> &flux,
             vector<double> &flux_err) {
  ifstream file(filename);
  string line;

  if (!file.is_open()) {
    cerr << "Error opening file: " << filename << endl;
    return;
  }

  // Read the header line
  if (getline(file, line)) {
    // Do nothing with the header line
  }

  while (getline(file, line)) {
    stringstream lineStream(line);
    string cell;
    vector<string> row;

    while (getline(lineStream, cell, ',')) {
      row.push_back(cell);
    }

    if (row.size() >= 3) { // Ensure there are at least 3 columns
      time.push_back(stod(row[0]));
      flux.push_back(stod(row[1]));
      flux_err.push_back(stod(row[2]));
    }
  }

  file.close();
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <filename>" << endl;
    return 1;
  }

  string filename = argv[1];
  vector<double> time, flux, flux_err;

  readCSV(filename, time, flux, flux_err);

  int n = flux.size();

  thrust::host_vector<double> h_time = time;
  thrust::host_vector<double> h_flux = flux;
  thrust::host_vector<double> h_flux_err = flux_err;

  thrust::device_vector<double> d_time = h_time;
  thrust::device_vector<double> d_flux = h_flux;
  thrust::device_vector<double> d_flux_err = h_flux_err;
  thrust::device_vector<double> d_weight(n);

  double h_weight_sum = 0;
  double *d_weight_sum;

  cudaMalloc((void **)&d_weight_sum, sizeof(double));
  cudaMemcpy(d_weight_sum, &h_weight_sum, sizeof(double),
             cudaMemcpyHostToDevice);

  int blockSize = 256;
  int numBlocks = (n + blockSize - 1) / blockSize;
  if (n % blockSize != 0) {
    numBlocks++;
  }

  dim3 dimGrid(numBlocks, 1, 1);
  dim3 dimBlock(blockSize, 1, 1);

  auto start = chrono::high_resolution_clock::now();

  weight_sum_kernel<<<dimGrid, dimBlock>>>(
      thrust::raw_pointer_cast(d_flux_err.data()), d_weight_sum, n);
  cudaMemcpy(&h_weight_sum, d_weight_sum, sizeof(double),
             cudaMemcpyDeviceToHost);

  h_weight_sum = pow(h_weight_sum, -1);
  compute_weight<<<dimGrid, dimBlock>>>(
      thrust::raw_pointer_cast(d_flux_err.data()),
      thrust::raw_pointer_cast(d_weight.data()), h_weight_sum, n);

  BlsResult h_result;
  h_result.d_value = LONG_MAX;
  BlsResult *d_result;

  cudaMalloc((void **)&d_result, sizeof(BlsResult));
  cudaMemcpy(d_result, &h_result, sizeof(BlsResult), cudaMemcpyHostToDevice);

  bls_kernel<<<dimGrid, dimBlock>>>(thrust::raw_pointer_cast(d_time.data()),
                                    thrust::raw_pointer_cast(d_flux.data()),
                                    thrust::raw_pointer_cast(d_weight.data()),
                                    n, d_result);
  cudaMemcpy(&h_result, d_result, sizeof(BlsResult), cudaMemcpyDeviceToHost);

  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::seconds>(end - start);

  cout << "Period: " << h_result.period << endl;
  cout << "D Value: " << h_result.d_value << endl;
  cout << "Time: " << duration.count() << " s" << endl;

  cudaFree(d_weight_sum);
  cudaFree(d_result);

  return 0;
}
