#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <vector>

using namespace std;

struct BlsResult {
    double period;
    double d_value;
};

struct WeightSum {
    __host__ __device__ double operator()(const double &x) const {
        return pow(x, -2);
    }
};

struct WeightCalculation {
    double sum_w;
    WeightCalculation(double _sum_w) : sum_w(_sum_w) {}

    __host__ __device__ double operator()(const double &x) const {
        return sum_w * pow(x, -2);
    }
};

struct RValue {
    double *weight;
    int i1, i2;
    RValue(double *_weight, int _i1, int _i2) : weight(_weight), i1(_i1), i2(_i2) {}

    __host__ __device__ double operator()(const int &x) const {
        double r = 0;
        for (int i = i1; i <= i2; i++) {
            r += weight[i];
        }
        return r;
    }
};

struct SValue {
    double *weight;
    double *flux;
    int i1, i2;
    SValue(double *_weight, double *_flux, int _i1, int _i2) : weight(_weight), flux(_flux), i1(_i1), i2(_i2) {}

    __host__ __device__ double operator()(const int &x) const {
        double s = 0;
        for (int i = i1; i <= i2; i++) {
            s += (weight[i] * flux[i]);
        }
        return s;
    }
};

struct DValue {
    double *weight;
    double *flux;
    double r, s;
    int size;
    DValue(double *_weight, double *_flux, double _r, double _s, int _size) : weight(_weight), flux(_flux), r(_r), s(_s), size(_size) {}

    __host__ __device__ double operator()(const int &x) const {
        double aux = 0;
        for (int i = 0; i < size; i++) {
            aux += (weight[i] * pow(flux[i], 2));
        }
        return aux - pow(s, 2) / (r * (1 - r));
    }
};

BlsResult my_bls(thrust::device_vector<double> &time, thrust::device_vector<double> &flux, thrust::device_vector<double> &flux_err) {
    BlsResult result;
    result.d_value = DBL_MAX;

    double sum_w = thrust::transform_reduce(flux_err.begin(), flux_err.end(), WeightSum(), 0.0, thrust::plus<double>());

    thrust::device_vector<double> weight(flux.size());
    thrust::transform(flux_err.begin(), flux_err.end(), weight.begin(), WeightCalculation(sum_w));

    for (int i1 = 0; i1 < flux.size(); i1++) {
        for (int i2 = i1 + 1; i2 < flux.size(); i2++) {
            double r_ = thrust::transform_reduce(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(flux.size()), RValue(thrust::raw_pointer_cast(weight.data()), i1, i2), 0.0, thrust::plus<double>());
            double s_ = thrust::transform_reduce(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(flux.size()), SValue(thrust::raw_pointer_cast(weight.data()), thrust::raw_pointer_cast(flux.data()), i1, i2), 0.0, thrust::plus<double>());
            double d_ = thrust::transform_reduce(thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(flux.size()), DValue(thrust::raw_pointer_cast(weight.data()), thrust::raw_pointer_cast(flux.data()), r_, s_, flux.size()), 0.0, thrust::plus<double>());

            double period = (time[i2] - time[i1]);

            if (d_ < result.d_value) {
                result.d_value = d_;
                result.period = period;
            }

            cout << "i1: " << i1 << " i2: " << i2 << " period: " << period << " d: " << d_ << endl;
        }
    }

    return result;
}

void readCSV(const string &filename, thrust::host_vector<double> &time, thrust::host_vector<double> &flux, thrust::host_vector<double> &flux_err) {
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

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

        if (row.size() >= 3) {
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
    thrust::host_vector<double> h_time, h_flux, h_flux_err;

    readCSV(filename, h_time, h_flux, h_flux_err);

    thrust::device_vector<double> d_time = h_time;
    thrust::device_vector<double> d_flux = h_flux;
    thrust::device_vector<double> d_flux_err = h_flux_err;

    auto start = chrono::high_resolution_clock::now();

    BlsResult result = my_bls(d_time, d_flux, d_flux_err);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);

    cout << "Period: " << result.period << endl;
    cout << "D Value: " << result.d_value << endl;
    cout << "Time: " << duration.count() << " s" << endl;

    return 0;
}
