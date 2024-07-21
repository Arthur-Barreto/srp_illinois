#include "utils.h"
#include <mpi.h>
using namespace std;

BLSResult bls_mpi(
    vector<double> &time,
    vector<double> &flux,
    vector<double> &flux_err,
    vector<SPECParameters> &s_params,
    int rank,
    int size);

void minloc_reduction(void *in, void *inout, int *len, MPI_Datatype *dptr);

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) {
            cout << "Usage: " << argv[0] << " <filename>" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    double real_period = 50.0;
    double real_phase = 5.0;
    double real_duration = 0.1 * real_period;
    double real_diff = 0.05;

    double threshold = cosf64x(M_PI * real_duration / real_period);
    vector<double> time = linspace(0, 400, 5000);
    vector<double> flux(time.size());

    for (size_t i = 0; i < time.size(); ++i) {
        flux[i] = cosf64x(2.0 * M_PI * (time[i] - real_phase - real_duration / 2.0) / real_period) > threshold ? 1.0 : 0.0;
    }

    for (size_t i = 0; i < time.size(); ++i) {
        flux[i] = 1.0 - real_diff * flux[i];
    }
    vector<double> flux_err(time.size(), 0.01);

    auto start = chrono::high_resolution_clock::now();

    vector<SPECParameters> s_params = spec_generator_gambiarra_omp(time);
    BLSResult result = bls_mpi(time, flux, flux_err, s_params, rank, size);

    MPI_Finalize();

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

    if (rank == 0) {
        cout << "Best period: " << result.best_period << endl;
        cout << "Best duration: " << result.best_duration << endl;
        cout << "Best phase: " << result.best_phase << endl;
        cout << "Best d_value: " << result.best_d_value << endl;
        cout << "Execution time: " << duration.count() << " ms" << endl;
    }

    return 0;
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

    vector<double> t_rel = compute_trel_omp(time);
    vector<double> normalized_flux = normalize_omp(flux);
    vector<double> weights = compute_weights_omp(flux_err);

    BLSResult local_result;
    local_result.best_d_value = DBL_MAX;

    // Distribute the workload among MPI processes
    size_t num_params = s_params.size();
    size_t start_idx = (rank * num_params) / size;
    size_t end_idx = ((rank + 1) * num_params) / size;

    double local_best_d_value = DBL_MAX;
    double local_best_period = 0.0;
    double local_best_duration = 0.0;
    double local_best_phase = 0.0;

#pragma omp parallel for
    for (size_t i = start_idx; i < end_idx; ++i) {
        double period = get<0>(s_params[i]);
        double duration = get<1>(s_params[i]);
        double phase = get<2>(s_params[i]);

        double d_value = model_omp(t_rel, normalized_flux, weights, period, duration, phase);

#pragma omp critical
        {
            if (d_value < local_best_d_value) {
                local_best_d_value = d_value;
                local_best_period = period;
                local_best_duration = duration;
                local_best_phase = phase;
            }
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