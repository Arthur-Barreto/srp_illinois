#include "utils.h"
using namespace std;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    cout << std::fixed << setprecision(16);

    string filename = argv[1];
    // vector<double> time, flux, flux_err;

    // readCSV(filename, time, flux, flux_err);
    // create fake data
    double real_period = 50.0;
    double real_phase = 5.0;
    double real_duration = 0.1 * real_period;
    double real_diff = 0.05;

    double threshold = cosf64x(M_PI * real_duration / real_period);
    vector<double> time = linspace(0, 400, 4000);
    vector<double> flux(time.size());

    for (size_t i = 0; i < time.size(); ++i) {
        flux[i] = cosf64x(2.0 * M_PI * (time[i] - real_phase - real_duration / 2.0) / real_period) > threshold ? 1.0 : 0.0;
    }

    for (size_t i = 0; i < time.size(); ++i) {
        flux[i] = 1.0 - real_diff * flux[i];
    }
    vector<double> flux_err(time.size(), 0.01);
    // end of fake data

    auto start = chrono::high_resolution_clock::now();

    vector<SPECParameters> s_params = spec_generator_gambiarra(time);

    BLSResult result = bls(time, flux, flux_err, s_params);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << "Best period: " << result.best_period << endl;
    cout << "Best duration: " << result.best_duration << endl;
    cout << "Best phase: " << result.best_phase << endl;
    cout << "Best d_value: " << result.best_d_value << endl;

    cout << endl << "Execution time: " << duration.count() << " ms" << endl;

    return 0;
}