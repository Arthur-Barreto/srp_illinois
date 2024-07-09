#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <float.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

struct BLSResult {
    double best_period;
    double best_phase;
    double best_duration;
};

template <typename T>
T min_value(const vector<T> &v) {
    return *min_element(v.begin(), v.end());
}

template <typename T>
T max_value(const vector<T> &v) {
    return *max_element(v.begin(), v.end());
}

template <typename T>
T ptp(const vector<T> &v) {
    T min_val = min_value(v);
    T max_val = max_value(v);
    return max_val - min_val;
}

vector<double> logspace(double start, double end, int num) {
    vector<double> result(num);
    double log_start = log10(start);
    double log_end = log10(end);
    double step = (log_end - log_start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = pow(10, log_start + i * step);
    }
    return result;
}

vector<double> auto_period(
    const vector<double> &t,
    const vector<double> &y,
    const vector<double> &dy,
    const double &duration,
    double minimum_period = -1,
    double maximum_period = -1,
    int minimum_n_transit = 3,
    double frequency_factor = 1.0) {
    // Check that the input arrays have the same length.
    assert(t.size() == y.size() && y.size() == dy.size() && "t, y, and dy must have the same length");

    // Check that all durations are positive.
    // assert(all_of(duration.begin(), duration.end(), [](double d) { return d > 0; }) && "duration must be positive");

    // Check that the minimum period is positive.
    if (minimum_period == -1) {
        minimum_period = 2.0 * duration;
    } else if (minimum_period > 0) {
        assert(minimum_period > 0 && "minimum_period must be positive");
    } else {
        throw invalid_argument("minimum_period must be positive");
    }

    // Compute the baseline.
    double baseline = ptp(t);
    double min_duration = duration;

    // Check that the maximum period is positive.
    if (maximum_period == -1) {
        maximum_period = baseline;
    } else if (maximum_period > 0) {
        assert(maximum_period > 0 && "maximum_period must be positive");
    } else {
        throw invalid_argument("maximum_period must be positive");
    }

    // Check that the minimum period is less than the maximum period.
    assert(minimum_period < maximum_period && "minimum_period must be less than maximum_period");

    // Check that the minimum number of transits is positive.
    assert(minimum_n_transit > 0 && "minimum_n_transit must be positive");

    // Check that the frequency factor is positive.
    assert(frequency_factor > 0 && "frequency_factor must be positive");

    // Estimate the required frequency spacing
    double df = frequency_factor * min_duration / (baseline * baseline);

    // Convert bounds to frequency
    double minimum_frequency = 1.0 / maximum_period;
    double maximum_frequency = 1.0 / minimum_period;

    // Compute the number of frequencies and the frequency grid
    int nf = 1 + static_cast<int>(round((maximum_frequency - minimum_frequency) / df));
    vector<double> frequencies(nf);
    for (int i = 0; i < nf; ++i) {
        frequencies[i] = 1.0 / (maximum_frequency - df * i);
    }

    return frequencies;
}

void fit(
    vector<double> &flux,
    vector<double> &weights) {

    // fit data (flux)
    double mean = accumulate(flux.begin(), flux.end(), 0.0) / flux.size();
    for (size_t i = 0; i < flux.size(); ++i) {
        flux[i] -= mean;
    }

    double std = 0.0;
    for (size_t i = 0; i < flux.size(); ++i) {
        std += flux[i] * flux[i];
    }

    std = sqrt(std / flux.size());
    for (size_t i = 0; i < flux.size(); ++i) {
        flux[i] /= std;
    }

    // fit weight vector
    for (size_t i = 0; i < flux.size(); ++i) {
        weights[i] = 1.0 / (flux[i] * flux[i]);
    }

    double sum_weights = accumulate(weights.begin(), weights.end(), 0.0);
    for (size_t i = 0; i < flux.size(); ++i) {
        weights[i] /= sum_weights;
    }

    return;
}

vector<double> auto_phase(double period, double duration) {
    double delta_t = int(period / duration);
    // create a vector that goes [0,period[, wirth delta_t steps
    vector<double> phase(delta_t);
    for (size_t i = 0; i < delta_t; ++i) {
        phase[i] = i * duration / period;
    }
    return phase;
}

double model(
    const vector<double> &time,
    const vector<double> &flux,
    const vector<double> &weights,
    const double &periods,
    const double &durations,
    const double &phases) {

    vector<double> in_transit(time.size(), 0.0);
    for (size_t i = 0; i < time.size(); ++i) {
        double t_mod_p = fmod(time[i], periods);
        double phase_diff = t_mod_p - phases;
        if (phase_diff >= 0 && phase_diff <= durations) {
            in_transit[i] = 1.0;
        }
    }

    double r = 0.0;
    double s = 0.0;
    double wx = 0.0;
    for (size_t i = 0; i < time.size(); ++i) {
        r += weights[i] * in_transit[i];
        s += weights[i] * flux[i] * in_transit[i];
        wx += weights[i] * flux[i] * flux[i];
    }

    r = r + DBL_MAX * (r * (1 - r));

    return wx - (s * s) / (r * (1 - r));
}

BLSResult bls(
    const vector<double> &time,
    vector<double> &flux,
    const vector<double> &flux_err,
    const vector<double> &durations) {

    // fit flux and create weights
    vector<double> weights(flux.size(), 0.0);
    fit(flux, weights);

    BLSResult result;
    double best_d = DBL_MAX;

    for (auto &d : durations) {
        vector<double> periods = auto_period(time, flux, flux_err, d);
        for (auto &p : periods) {
            vector<double> phases = auto_phase(p, d);
            for (auto &phi : phases) {

                double d_val = model(time, flux, weights, p, d, phi);

                if (d_val < best_d) {
                    best_d = d_val;
                    result.best_period = p;
                    result.best_phase = phi;
                    result.best_duration = d;
                }
            }
        }
    }

    return result;
}

void readCSV(const string &filename, vector<double> &time, vector<double> &flux, vector<double> &flux_err) {
    ifstream file(filename); // Change 'ifstream' to 'ifstream'
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
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    string filename = argv[1];
    vector<double> time, flux, flux_err;

    readCSV(filename, time, flux, flux_err);

    vector<double> durations(10000);
    double start_duration = 20.0 * 0.01;
    double end_duration = 20.0 * 0.05;
    double step_duration = (end_duration - start_duration) / (durations.size() - 1);
    for (size_t i = 0; i < durations.size(); ++i) {
        durations[i] = start_duration + i * step_duration;
    }

    auto start = chrono::high_resolution_clock::now();

    BLSResult result = bls(time, flux, flux_err, durations);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    cout << "Best period: " << result.best_period << endl;
    cout << "Best phase: " << result.best_phase << endl;
    cout << "Best duration: " << result.best_duration << endl;
    cout << "Execution time: " << duration << " ms" << endl;

    return 0;
}
