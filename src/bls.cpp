#include "utils.h"
using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " <n_points> <n_samples>" << endl;
    return 1;
  }

  cout << std::fixed << setprecision(16);

  size_t n_points = stoul(argv[1]);
  size_t n_samples = stoul(argv[2]);

  // fake data
  double real_period = 50.0;
  double real_phase = 5.0;
  double real_duration = 0.1 * real_period;
  double real_diff = 0.05;

  double threshold = cosf64x(M_PI * real_duration / real_period);
  vector<double> time = linspace(0, 400, n_points);
  vector<double> flux(time.size());

  for (size_t i = 0; i < time.size(); ++i) {
    flux[i] =
        cosf64x(2.0 * M_PI * (time[i] - real_phase - real_duration / 2.0) /
                real_period) > threshold
            ? 1.0
            : 0.0;
  }

  for (size_t i = 0; i < time.size(); ++i) {
    flux[i] = 1.0 - real_diff * flux[i];
  }
  vector<double> flux_err(time.size(), 0.01);
  // end of fake data

  auto start = chrono::high_resolution_clock::now();
  BLSResult result;
  double total_duration = 0.0;

  for (int i = 0; i < n_samples; ++i) {
    auto iteration_start = chrono::high_resolution_clock::now();
    vector<SPECParameters> s_params = spec_generator_gambiarra(time);
    result = bls(time, flux, flux_err, s_params);
    auto iteration_end = chrono::high_resolution_clock::now();
    total_duration += chrono::duration_cast<chrono::milliseconds>(
                          iteration_end - iteration_start)
                          .count();
  }

  auto average_duration = total_duration / n_samples;

  cout << "Best period: " << result.best_period << endl;
  cout << "Best duration: " << result.best_duration << endl;
  cout << "Best phase: " << result.best_phase << endl;
  cout << "Best d_value: " << result.best_d_value << endl;
  cout << "Average execution time: " << average_duration << " ms" << endl;

  return 0;
}