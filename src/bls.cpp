#include "utils.h"
using namespace std;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " <filename>" << endl;
        return 1;
    }

    string filename = argv[1];
    vector<double> time, flux, flux_err;

    readCSV(filename, time, flux, flux_err);

    auto start = chrono::high_resolution_clock::now();

    vector<SPECParameters> s_params = spec_generator(time);
    BLSResult result = bls(time, flux, flux_err, s_params);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << "Best period: " << result.best_period << endl;
    cout << "Best duration: " << result.best_duration << endl;
    cout << "Best phase: " << result.best_phase << endl;
    cout << "Best d_value: " << result.best_d_value << endl;

    return 0;
}