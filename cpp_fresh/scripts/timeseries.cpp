#include "timeseries.h"
#include <iostream>
#include <random>
#include <cmath>

std::vector<double> create_period(int p_length) {
    std::vector<double> period(p_length);
    for (int i = 0; i < p_length; ++i) {
        period[i] = std::sin((i * 2 * M_PI) / p_length);
    }
    return period;
}

std::vector<double> create_ts(int n_points, int p_length, double sigma) {
    std::vector<double> pts(n_points);
    std::vector<double> pattern = create_period(p_length);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, sigma);
    for (int i = 0; i < n_points / p_length; ++i) {
        int j = p_length * i;
        for (int k = 0; k < p_length; ++k) {
            double error = distribution(generator);
            pts[j + k] = pattern[k] + error;
        }
    }
    return pts;
}

std::vector<std::vector<double>> shingle(const std::vector<double>& series, int dim) {
    int height = series.size() - dim + 1;
    std::vector<std::vector<double>> shingled(dim, std::vector<double>(height));
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < height; ++j) {
            shingled[i][j] = series[j + i];
        }
    }
    return shingled;
}

std::pair<std::vector<double>, std::vector<int>> inject_anomalies(std::vector<double> pts, int l_anom, int n_anom) {
    std::vector<int> loc(n_anom);
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distribution(0, pts.size() - l_anom - 1);
    for (int i = 0; i < n_anom; ++i) {
        loc[i] = distribution(generator);
        for (int j = 1; j < l_anom; ++j) {
            pts[loc[i] + j] = pts[loc[i]];
        }
    }
    return std::make_pair(pts, loc);
}