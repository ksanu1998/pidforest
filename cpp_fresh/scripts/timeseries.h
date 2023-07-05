#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>

std::vector<double> create_period(int p_length);
std::vector<double> create_ts(int n_points, int p_length, double sigma);
std::vector<std::vector<double>> shingle(const std::vector<double>& series, int dim);
std::pair<std::vector<double>, std::vector<int>> inject_anomalies(std::vector<double> pts, int l_anom, int n_anom);
void print_anomalies(const std::vector<double>& pts, const std::vector<int>& indices);

#endif  // TIMESERIES_H
