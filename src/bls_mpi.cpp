#include "utils_mpi.h"
using namespace std;

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) {
            cout << "Usage: " << argv[0] << " <n_points> <samples>" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    size_t n_points = stoul(argv[1]);
    size_t samples = stoul(argv[2]);

    double real_period = 50.0;
    double real_phase = 5.0;
    double real_duration = 0.1 * real_period;
    double real_diff = 0.05;

    double threshold = cosf64x(M_PI * real_duration / real_period);
    vector<double> time = linspace(0, 400, n_points);
    vector<double> flux(time.size());

    for (size_t i = 0; i < time.size(); ++i) {
        flux[i] = cosf64x(2.0 * M_PI * (time[i] - real_phase - real_duration / 2.0) / real_period) > threshold ? 1.0 : 0.0;
    }

    for (size_t i = 0; i < time.size(); ++i) {
        flux[i] = 1.0 - real_diff * flux[i];
    }
    vector<double> flux_err(time.size(), 0.01);

    cout << fixed << setprecision(10);

    double total_duration = 0.0;
    BLSResult result;
    // auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < samples; ++i) {
        auto iteration_start = chrono::high_resolution_clock::now();
        vector<SPECParameters> s_params = spec_generator_gambiarra(time);
        result = bls_mpi(time, flux, flux_err, s_params, rank, size);
        auto iteration_end = chrono::high_resolution_clock::now();
        total_duration += chrono::duration_cast<chrono::milliseconds>(iteration_end - iteration_start).count();
    }

    // auto end = chrono::high_resolution_clock::now();
    auto average_duration = total_duration / samples;

    if (rank == 0) {
        cout << "Best period: " << result.best_period << endl;
        cout << "Best duration: " << result.best_duration << endl;
        cout << "Best phase: " << result.best_phase << endl;
        cout << "Best d_value: " << result.best_d_value << endl;
        cout << "Average execution time: " << average_duration << " ms" << endl;
    }

    MPI_Finalize();
    return 0;
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
