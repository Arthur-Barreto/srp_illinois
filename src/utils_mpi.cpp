#include "utils_mpi.h"

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

vector<double> auto_phase(double period, double duration) {
    vector<double> phase(int(ceil(period / duration)) + 1);
    phase = linspace(0, period, int(ceil(period / duration)) + 1);
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
    frequency_resolution *= frequency_resolution;

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

    double step = ptp(time) / time.size();

    for (const auto &p : periods) {
        vector<double> durations = arange(0.01 * p, 0.05 * p, step);
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

double model_mpi(
    vector<double> &t_rel,
    vector<double> &flux,
    vector<double> &weights,
    double period,
    double duration,
    double phase,
    int rank,
    int size) {

    size_t n = flux.size();
    size_t start = ((rank * n) / size);
    size_t end = (((rank + 1) * n) / size) - 1;

    vector<double> local_t_rel(t_rel.begin() + start, t_rel.begin() + end);
    vector<double> local_flux(flux.begin() + start, flux.begin() + end);
    vector<double> local_weights(weights.begin() + start, weights.begin() + end);
    vector<size_t> local_is_transit(local_t_rel.size());

    // Compute local is_transit
    for (size_t i = 0; i < local_t_rel.size(); ++i) {
        local_is_transit[i] = ((fmod(local_t_rel[i], period) >= phase) && fmod(local_t_rel[i], period) <= phase + duration) ? 1 : 0;
    }

    // Compute local r, s, and wx
    double local_r = inner_product(local_weights.begin(), local_weights.end(), local_is_transit.begin(), 0.0);
    double local_s = 0.0;
    for (size_t i = 0; i < local_flux.size(); ++i) {
        local_s += local_weights[i] * local_flux[i] * local_is_transit[i];
    }
    double local_wx = inner_product(local_weights.begin(), local_weights.end(), local_flux.begin(), 0.0, plus<double>(), [](double w, double x) { return w * x * x; });

    // Debugging: Print local_r, local_s, local_wx
    stringstream ss;
    ss << "Rank " << rank << " [" << start << "," << end << "]" << endl
       << "  t_rel[" << start << "] = " << t_rel[start] << endl
       << "  t_rel[" << end << "] = " << t_rel[end] << endl
       << "  next: ";
    if (end < t_rel.size()) {
        ss << t_rel[end + 1];
    } else {
        ss << "None";
    }
    ss << endl;
    ss << "  local_is_transit[" << 0 << "] = " << local_is_transit[0] << endl
       << "  local_is_transit[" << local_is_transit.size() << "] = " << local_is_transit[local_is_transit.size() - 1] << endl
       << "  next: ";
    if (local_is_transit.size() - 1 < local_is_transit.size()) {
        ss << local_is_transit[local_is_transit.size() - 1];
    } else {
        ss << "None";
    }
    ss << endl;
    ss << "    local_r: " << local_r << ", local_s: " << local_s << ", local_wx: " << local_wx << endl;
    cout << ss.str();
    // ss << "";
    // ss << "RANK " << rank << " size = " << end - start + 1 << endl;
    // ss << "t_rel[" << start << "]= " << t_rel[start] << endl;
    // ss << "t_rel[" << start + 1 << "]= " << t_rel[start + 1] << endl;
    // ss << "t_rel[" << end - 1 << "]= " << t_rel[end - 1] << endl;
    // ss << "t_rel[" << end << "]= " << t_rel[end] << endl;
    // ss << endl;
    // cout << ss.str();

    // Synchronize before reduction
    MPI_Barrier(MPI_COMM_WORLD);

    // Reduce the results to the root process
    double r = 0.0, s = 0.0, wx = 0.0;

    MPI_Reduce(&local_r, &r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_s, &s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_wx, &wx, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double d_value = 0.0;
    if (rank == 0) {
        d_value = wx - (s * s) / (r * (1 - r)) + numeric_limits<double>::epsilon();
        cout << "Global r: " << r << ", s: " << s << ", wx: " << wx << endl;
    }

    // Broadcast the result back to all processes
    MPI_Bcast(&d_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return d_value;
}

void minloc_reduction(void *in, void *inout, int *len, MPI_Datatype *dptr) {
    BLSResult *in_vals = static_cast<BLSResult *>(in);
    BLSResult *inout_vals = static_cast<BLSResult *>(inout);

    for (int i = 0; i < *len; ++i) {
        if (in_vals[i].best_d_value < inout_vals[i].best_d_value) {
            inout_vals[i] = in_vals[i];
        }
    }
}

BLSResult bls_mpi(
    vector<double> &time,
    vector<double> &flux,
    vector<double> &flux_err,
    vector<SPECParameters> &s_params,
    int rank,
    int size) {

    vector<double> t_rel = compute_trel(time);
    vector<double> normalized_flux = normalize(flux);
    vector<double> weights = compute_weights(flux_err);

    BLSResult local_result;
    local_result.best_d_value = DBL_MAX;

    size_t num_params = s_params.size();
    size_t start_idx = (rank * num_params) / size;
    size_t end_idx = ((rank + 1) * num_params) / size;

    double local_best_d_value = DBL_MAX;
    double local_best_period = 0.0;
    double local_best_duration = 0.0;
    double local_best_phase = 0.0;

    for (size_t i = start_idx; i < end_idx; ++i) {
        double period = get<0>(s_params[i]);
        double duration = get<1>(s_params[i]);
        double phase = get<2>(s_params[i]);

        double d_value = model(t_rel, normalized_flux, weights, period, duration, phase);

        if (d_value < local_best_d_value) {
            local_best_d_value = d_value;
            local_best_period = period;
            local_best_duration = duration;
            local_best_phase = phase;
        }
    }

    local_result.best_d_value = local_best_d_value;
    local_result.best_period = local_best_period;
    local_result.best_duration = local_best_duration;
    local_result.best_phase = local_best_phase;

    // Define a custom MPI datatype for the BLSResult struct
    MPI_Datatype MPI_BLSResult;
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int block_lengths[4] = {1, 1, 1, 1};
    MPI_Aint offsets[4];

    offsets[0] = offsetof(BLSResult, best_d_value);
    offsets[1] = offsetof(BLSResult, best_period);
    offsets[2] = offsetof(BLSResult, best_duration);
    offsets[3] = offsetof(BLSResult, best_phase);

    MPI_Type_create_struct(4, block_lengths, offsets, types, &MPI_BLSResult);
    MPI_Type_commit(&MPI_BLSResult);

    // Define a custom reduction operation
    MPI_Op MPI_BLSResult_minloc;
    MPI_Op_create((MPI_User_function *)minloc_reduction, true, &MPI_BLSResult_minloc);

    BLSResult global_result;
    MPI_Reduce(&local_result, &global_result, 1, MPI_BLSResult, MPI_BLSResult_minloc, 0, MPI_COMM_WORLD);

    MPI_Type_free(&MPI_BLSResult);
    MPI_Op_free(&MPI_BLSResult_minloc);

    return global_result;
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