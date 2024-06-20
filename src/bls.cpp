#include <algorithm>
#include <cmath>
#include <map>
#include <stdlib.h>
#include <vector>

using namespace std;

double weight_sum(vector<double> flux_err) {
    double weight = 0;
    for (int i = 0; i < flux_err.size(); i++) {
        weight += 1.0 / pow(flux_err[i], -2);
    }
    return weight;
}

double r_value(int i1, int i2, vector<double> weight) {
    double r = 0;
    for (int i = i1; i <= i2; i++) {
        r += weight[i];
    }
    return r;
}

double s_value(int i1, int i2, vector<double> weight, vector<double> flux) {
    double s = 0;
    for (int i = i1; i <= i2; i++) {
        s += (weight[i] * flux[i]);
    }
    return s;
}

double d_value(vector<double> weight, vector<double> flux, double r, double s) {
    double aux = 0;

    for (int i = 0; i < weight.size(); i++) {
        aux += (weight[i] * pow(flux[i], 2));
    }

    return aux - pow(s, 2) / (r * (1 - r));
}

map<double, double> my_bls(vector<double> time, vector<double> flux, vector<double> flux_err) {

    double sum_w = weight_sum(flux_err);

    vector<double> weight(flux.size());
    for (int i = 0; i < flux.size(); i++) {
        weight[i] = sum_w * pow(flux_err[i], -2);
    }

    map<double, double> periodogram;
    double period;

    for (int i1 = 0; i1 < flux.size(); i1++) {
        for (int i2 = i1 + 1; i2 < flux.size(); i2++) {
            double r_ = r_value(i1, i2, weight);
            double s_ = s_value(i1, i2, weight, flux);
            double d_ = d_value(weight, flux, r_, s_);

            period = (time[i2] - time[i1]);
            periodogram.insert({d_, period});
        }
    }

    return periodogram;
}