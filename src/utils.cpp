#include "utils.h";
using namespace std;

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

vector<double> arange(double start, double end, double step) {
    vector<double> result;
    for (double value = start; value < end; value += step) {
        result.push_back(value);
    }
    return result;
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

PERIODParameters auto_max_min_period(vector<double> &time) {
    double total_duration = ptp(time);

    // compute the difference between each pair of time values
    vector<double> dt(time.size() - 1);
    for (size_t i = 0; i < time.size() - 1; ++i) {
        dt[i] = time[i + 1] - time[i];
    }

    // compute the minimum period
    double min_period = 2.0 * min_value(dt);

    // compute the maximum period
    double max_period = 0.5 * total_duration;

    PERIODParameters params;
    params.minimum_period = min_period;
    params.maximum_period = max_period;
    params.total_duration = total_duration;

    return params;
}

vector<double> auto_period(
    double minimum_period = -1,
    double maximum_period = -1,
    double total_duration = -1) {

    // Check that the minimum and maxmimum period is positive.
    assert(minimum_period > 0 && "minimum_period must be positive");
    assert(maximum_period > 0 && "maximum_period must be positive");

    // check the minimum period is less than the maximum period.
    assert(minimum_period < maximum_period && "minimum_period must be less than maximum_period");

    // check the maximum period does not exceed the total duration.
    assert(maximum_period <= total_duration && "maximum_period must be less than or equal to total_duration");

    // Convert bounds to frequency
    double minimum_frequency = 1.0 / maximum_period;
    double maximum_frequency = 1.0 / minimum_period;
    double frequency_resolution = 1.0 / total_duration;

    // Compute the number of frequencies and the frequency grid
    size_t nf = 1 + static_cast<int>(round((maximum_frequency - minimum_frequency) / frequency_resolution));
    vector<double> periods(nf);
    for (size_t i = 0; i < nf; ++i) {
        periods[i] = 1.0 / (maximum_frequency - frequency_resolution * i);
    }

    return periods;
}

std::vector<SPECParameters> spec_generator(std::vector<double> &time) {

    PERIODParameters p_params = auto_max_min_period(time);
    vector<double> periods = auto_period(p_params.minimum_period, p_params.maximum_period, p_params.total_duration);

    vector<SPECParameters> spec_params;

    for (const auto &p : periods) {
        vector<double> durations = arange(0.01 * p, 0.05 * p, 1e-5);
        for (const auto &d : durations) {
            vector<double> phases = auto_phase(p, d);
            for (const auto &phi : phases) {
                spec_params.push_back(make_tuple(p, d, phi));
            }
        }
    }

    return spec_params;
}

void compute_trel(vector<double> &time) {
    double t0 = min_value(time);
    for (auto &t : time) {
        t -= t0;
    }
}

void normalize(vector<double> &flux) {
    double mean_flux = accumulate(flux.begin(), flux.end(), 0.0) / flux.size();
    for (auto &f : flux) {
        f /= mean_flux + numeric_limits<double>::epsilon();
    }
}

vector<double> compute_weights(vector<double> &flux_err) {
    vector<double> weights(flux_err.size());
    for (size_t i = 0; i < flux_err.size(); ++i) {
        weights[i] = 1.0 / (flux_err[i] * flux_err[i]);
    }

    double sum_weights = accumulate(weights.begin(), weights.end(), 0.0);

    assert(fabs(sum_weights) > numeric_limits<double>::epsilon() && "Sum of weights must be greater than epsilon");

    for (auto &w : weights) {
        w /= sum_weights;
    }

    return weights;
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

        if (row.size() >= 3) {
            time.push_back(stod(row[0]));
            flux.push_back(stod(row[1]));
            flux_err.push_back(stod(row[2]));
        }
    }

    file.close();
}