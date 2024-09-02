#include "utils_omp.h"
using namespace std;

int main(int argc, char *argv[]) {
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " <n_points> <n_threads> <n_samples>"
         << endl;
    return 1;
  }

  cout << std::fixed << setprecision(16);

  // create fake data
  size_t n_points = stoul(argv[1]);
  size_t n_threads = stoul(argv[2]);
  size_t n_samples = stoul(argv[3]);

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

  omp_set_num_threads(n_threads);

  auto start = chrono::high_resolution_clock::now();
  BLSResult result;
  vector<double> exec_time = vector<double>(n_samples);

  for (int i = 0; i < n_samples; ++i) {
    auto iteration_start = chrono::high_resolution_clock::now();
    vector<SPECParameters> s_params = spec_generator_gambiarra(time);
    result = bls_omp(time, flux, flux_err, s_params);
    auto iteration_end = chrono::high_resolution_clock::now();
    exec_time[i] = chrono::duration_cast<chrono::milliseconds>(iteration_end -
                                                               iteration_start)
                       .count();
  }

  auto average_duration =
      accumulate(exec_time.begin(), exec_time.end(), 0.0) / n_samples;
  auto min_duration = *min_element(exec_time.begin(), exec_time.end());
  auto max_duration = *max_element(exec_time.begin(), exec_time.end());
  auto std_dev =
      sqrt(accumulate(exec_time.begin(), exec_time.end(), 0.0,
                      [average_duration](double accumulator, double value) {
                        return accumulator + pow(value - average_duration, 2);
                      }) /
           n_samples);

  cout << "Best period: " << result.best_period << endl;
  cout << "Best duration: " << result.best_duration << endl;
  cout << "Best phase: " << result.best_phase << endl;
  cout << "Best d_value: " << result.best_d_value << endl;
  cout << "Average execution time: " << average_duration << " ms" << endl;
  cout << "Min execution time: " << min_duration << " ms" << endl;
  cout << "Max execution time: " << max_duration << " ms" << endl;
  cout << "Standard deviation: " << std_dev << " ms" << endl;

  return 0;
}