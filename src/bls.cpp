#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <vector>

using namespace std;

struct BlsResult {
    long double period;
    long double d_value;
};

long double weight_sum(vector<long double> flux_err) {
    long double weight = 0;
    for (int i = 0; i < flux_err.size(); i++) {
        weight += pow(flux_err[i], -2);
    }
    return pow(weight, -1);
}

long double r_value(int i1, int i2, vector<long double> weight) {
    long double r = 0;
    for (int i = i1; i <= i2; i++) {
        r += weight[i];
    }
    return r;
}

long double s_value(int i1, int i2, vector<long double> weight, vector<long double> flux) {
    long double s = 0;
    for (int i = i1; i <= i2; i++) {
        s += (weight[i] * flux[i]);
    }
    return s;
}

long double d_value(vector<long double> weight, vector<long double> flux, long double r, long double s) {
    long double aux = 0;

    for (int i = 0; i < weight.size(); i++) {
        aux += (weight[i] * pow(flux[i], 2));
    }

    return aux - pow(s, 2) / (r * (1 - r));
}

BlsResult my_bls(vector<long double> time, vector<long double> flux, vector<long double> flux_err) {

    long double sum_w = weight_sum(flux_err);

    vector<long double> weight(flux.size());
    for (int i = 0; i < flux.size(); i++) {
        weight[i] = sum_w * pow(flux_err[i], -2);
    }

    BlsResult data;
    data.d_value = INT_MAX;

    long double period;

    for (int i1 = 0; i1 < flux.size(); i1++) {
        for (int i2 = i1 + 1; i2 < flux.size(); i2++) {
            long double r_ = r_value(i1, i2, weight);
            long double s_ = s_value(i1, i2, weight, flux);
            long double d_ = d_value(weight, flux, r_, s_);

            period = (time[i2] - time[i1]);

            if (d_ < data.d_value) {
                data.d_value = d_;
                data.period = period;
            }

            cout << "i1: " << i1 << " i2: " << i2 << " period: " << period << " d: " << d_ << endl;
        }
    }

    return data;
}

void readCSV(const string &filename, vector<long double> &time, vector<long double> &flux, vector<long double> &flux_err) {

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
    vector<long double> time, flux, flux_err;

    readCSV(filename, time, flux, flux_err);

    // init chonos time
    auto start = chrono::high_resolution_clock::now();

    BlsResult result = my_bls(time, flux, flux_err);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);

    cout << "Period: " << result.period << endl;
    cout << "D Value: " << result.d_value << endl;
    cout << "Time: " << duration.count() << " s" << endl;

    return 0;
}