#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <float.h>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdlib.h>
#include <vector>

using namespace std;

struct BlsResult {
    double period;
    double d_value;
};

double weight_sum(const vector<double> &flux_err) {
    double weight = accumulate(flux_err.begin(), flux_err.end(), 0.0, [](double sum, double val) {
        return sum + pow(val, -2);
    });
    return pow(weight, -1);
}

void precompute_cumulative_r_s(
    const vector<double> &weight,
    const vector<double> &flux,
    vector<double> &cumulative_r,
    vector<double> &cumulative_s) {
    cumulative_r.reserve(weight.size());
    cumulative_s.reserve(weight.size());

    double cr = weight[0];
    double cs = cr * flux[0];

    cumulative_r.push_back(cr);
    cumulative_s.push_back(cs);

    for (size_t i = 1; i < weight.size(); ++i) {
        const double &w = weight[i];
        const double &wf = w * flux[i];
        cr += w;
        cs += wf;
        cumulative_r.push_back(cr);
        cumulative_s.push_back(cs);
    }
}

double sum_wff_value(const vector<double> &weight, const vector<double> &flux) {
    return inner_product(weight.begin(), weight.end(), flux.begin(), 0.0, plus<>(), [](double w, double f) {
        return w * f * f;
    });
}

BlsResult my_bls(const vector<double> &time, const vector<double> &flux, const vector<double> &flux_err) {
    double sum_w = weight_sum(flux_err);

    vector<double> weight(flux.size());
    for (size_t i = 0; i < flux.size(); ++i) {
        weight[i] = sum_w * pow(flux_err[i], -2);
    }

    BlsResult data;
    data.d_value = DBL_MAX;

    vector<double> cumulative_r, cumulative_s;
    precompute_cumulative_r_s(weight, flux, cumulative_r, cumulative_s);

    const double aux = sum_wff_value(weight, flux);

    for (size_t i1 = 0; i1 < flux.size(); ++i1) {
        for (size_t i2 = i1 + 1; i2 < flux.size(); ++i2) {
            double r = cumulative_r[i2] - cumulative_r[i1] + weight[i1];
            double s = cumulative_s[i2] - cumulative_s[i1] + weight[i1] * flux[i1];
            double d = aux - (s * s) / (r * (1 - r));

            if (d < data.d_value) {
                data.d_value = d;
                data.period = time[i2] - time[i1];
            }
        }
    }

    return data;
}

void readCSV(const string &filename, vector<double> &time, vector<double> &flux, vector<double> &flux_err) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;
    getline(file, line); // Read the header line

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

    auto start = chrono::high_resolution_clock::now();

    BlsResult result = my_bls(time, flux, flux_err);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);

    cout << "Period: " << result.period << endl;
    cout << "D Value: " << result.d_value << endl;
    cout << "Time: " << duration.count() << " s" << endl;

    return 0;
}