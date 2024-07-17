#include "utils.h"
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

vector<double> arange(double start, double end, double step) {
    vector<double> result;
    for (double value = start; value < end; value += step) {
        result.push_back(value);
    }
    return result;
}

vector<double> linspace(double start, double end, size_t num) {
    vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (size_t i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

vector<double> logarithmic_space(double start, double end, size_t num) {
    vector<double> result(num);
    double step = (log10(end) - log10(start)) / (num - 1);
    for (size_t i = 0; i < num; ++i) {
        result[i] = pow(10, log10(start) + i * step);
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
    double minimum_period,
    double maximum_period,
    double total_duration) {

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
    size_t nf = static_cast<int>(round((maximum_frequency - minimum_frequency) / frequency_resolution));
    vector<double> periods(nf);
    double step = (maximum_frequency - minimum_frequency) / (nf - 1);
    for (size_t i = 0; i < nf; ++i) {
        periods[i] = 1.0 / (minimum_frequency + i * step);
    }

    return periods;
}

vector<SPECParameters> spec_generator(vector<double> &time) {

    PERIODParameters p_params = auto_max_min_period(time);
    vector<double> periods = auto_period(p_params.minimum_period, p_params.maximum_period, p_params.total_duration);

    vector<SPECParameters> spec_params;

    double frequency_resolution = 1.0 / ptp(time);
    frequency_resolution *= frequency_resolution;

    for (const auto &p : periods) {
        vector<double> durations = logarithmic_space(0.01 * p, 0.05 * p, time.size());
        for (const auto &d : durations) {
            vector<double> phases = auto_phase(p, d);
            for (const auto &phi : phases) {
                spec_params.push_back(make_tuple(p, d, phi));
            }
        }
    }

    return spec_params;
}

vector<SPECParameters> spec_generator_gambiarra(vector<double> &time) {
    PERIODParameters p_params = auto_max_min_period(time);
    vector<double> periods = auto_period(p_params.minimum_period, p_params.maximum_period, p_params.total_duration);

    vector<SPECParameters> spec_params;

    for (const auto &p : linspace(45, 55, 101)) {
        for (const auto &d : linspace(1, 11, 11)) {
            for (const auto &phi : arange(0, p, 0.5)) {
                spec_params.push_back(make_tuple(p, d, phi));
            }
        }
    }

    return spec_params;
}

vector<double> compute_trel(vector<double> &time) {
    double t0 = min_value(time);
    vector<double> t_rel(time.size(), 0.0);

    for (size_t i = 0; i < time.size(); ++i) {
        t_rel[i] = time[i] - t0;
    }

    return t_rel;
}

vector<double> normalize(vector<double> &flux) {
    vector<double> normalized_flux(flux.size());
    double mean_flux = accumulate(flux.begin(), flux.end(), 0.0) / flux.size();

    for (size_t i = 0; i < flux.size(); ++i) {
        normalized_flux[i] = flux[i] - mean_flux;
    }

    double squared_sum = inner_product(normalized_flux.begin(), normalized_flux.end(), normalized_flux.begin(), 0.0);

    double sigma = sqrt(squared_sum / flux.size());

    for (auto &f : normalized_flux) {
        f /= (sigma + numeric_limits<double>::epsilon());
    }

    return normalized_flux;
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

double model(
    vector<double> &t_rel,
    vector<double> &flux,
    vector<double> &weights,
    double period,
    double duration,
    double phase) {
    vector<size_t> is_transit(flux.size());
    for (size_t i = 0; i < flux.size(); ++i) {
        is_transit[i] = ((fmod(t_rel[i], period) >= phase) && fmod(t_rel[i], period) <= phase + duration) ? 1 : 0;
    }

    double r = inner_product(weights.begin(), weights.end(), is_transit.begin(), 0.0);

    double s = 0.0;
    for (size_t i = 0; i < flux.size(); ++i) {
        s += weights[i] * flux[i] * is_transit[i];
    }

    double wx = inner_product(weights.begin(), weights.end(), flux.begin(), 0.0, plus<>(), [](double w, double f) {
        return w * f * f;
    });

    double d_value = wx - (s * s) / (r * (1 - r)) + numeric_limits<double>::epsilon();

    return d_value;
}

BLSResult bls(
    vector<double> &time,
    vector<double> &flux,
    vector<double> &flux_err,
    vector<SPECParameters> &s_params) {

    vector<double> t_rel = compute_trel(time);
    vector<double> normalized_flux = normalize(flux);
    vector<double> weights = compute_weights(flux_err);

    BLSResult result;
    result.best_d_value = DBL_MAX;

    for (const auto &params : s_params) {
        double period = get<0>(params);
        double duration = get<1>(params);
        double phase = get<2>(params);

        // if period is equal to best period, go to the next iteration
        if (fabs(period - result.best_period) < numeric_limits<double>::epsilon()) {
            continue;
        }

        double d_value = model(t_rel, normalized_flux, weights, period, duration, phase);

        if (d_value < result.best_d_value) {
            result.best_d_value = d_value;
            result.best_period = period;
            result.best_duration = duration;
            result.best_phase = phase;
        }
    }

    return result;
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