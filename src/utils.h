#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <float.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct BLSResult {
    double best_period;
    double best_phase;
    double best_duration;
};

struct PERIODParameters {
    double minimum_period;
    double maximum_period;
    double total_duration;
};

typedef std::tuple<double, double, double> SPECParameters;

template <typename T>
T min_value(const std::vector<T> &v);

template <typename T>
T max_value(const std::vector<T> &v);

template <typename T>
T ptp(const std::vector<T> &v);

std::vector<double> logspace(double start, double end, int num);

std::vector<double> arange(double start, double end, double step);

std::vector<double> auto_phase(double period, double duration);

PERIODParameters auto_max_min_period(std::vector<double> &time);

std::vector<double> auto_period(
    double minimum_period = -1,
    double maximum_period = -1,
    double total_duration = -1);

std::vector<SPECParameters> spec_generator(std::vector<double> &time);

void compute_trel(std::vector<double> &time);

void normalize(std::vector<double> &flux);

std::vector<double> compute_weights(std::vector<double> &flux_err);

void readCSV(
    const string &filename,
    std::vector<double> &time,
    std::vector<double> &flux,
    std::vector<double> &flux_err);

#endif // UTILS_H