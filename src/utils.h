#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cfloat>      // Use <cfloat> for C++ compatibility
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// Struct for period parameters
struct PERIODParameters {
  double minimum_period;
  double maximum_period;
  double total_duration;
};

// Struct for the BLS result
struct BLSResult {
  double best_period;
  double best_duration;
  double best_phase;
  double best_d_value;
};

// Typedef for SPEC parameters (period, duration, phase)
typedef std::tuple<double, double, double> SPECParameters;

// Function declarations with const correctness for input vectors that are not modified
template <typename T>
T min_value(const std::vector<T> &v);

template <typename T>
T max_value(const std::vector<T> &v);

template <typename T>
T ptp(const std::vector<T> &v);

std::vector<double> arange(double start, double end, double step);

std::vector<double> linspace(double start, double end, size_t num);

std::vector<double> auto_phase(double period, double duration);

PERIODParameters auto_max_min_period(const std::vector<double> &time);

std::vector<double> auto_period(double minimum_period = -1,
                                double maximum_period = -1,
                                double total_duration = -1);

std::vector<SPECParameters> spec_generator(const std::vector<double> &time);

std::vector<SPECParameters> spec_generator_gambiarra(const std::vector<double> &time);

std::vector<double> compute_trel(const std::vector<double> &time);

std::vector<double> normalize(const std::vector<double> &flux);

std::vector<double> compute_weights(const std::vector<double> &flux_err);

double model(const std::vector<double> &t_rel, const std::vector<double> &flux,
             const std::vector<double> &weights, double period, double duration,
             double phase);

BLSResult bls(const std::vector<double> &time, const std::vector<double> &flux,
              const std::vector<double> &flux_err,
              const std::vector<SPECParameters> &s_params);

void readCSV(const std::string &filename, const std::vector<double> &time,
             const std::vector<double> &flux, const std::vector<double> &flux_err);

#endif  // SRC_UTILS_H_