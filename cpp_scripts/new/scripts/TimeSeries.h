#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>

class TimeSeries {
public:
    static std::vector<std::vector<double>> shingle(const std::vector<double>& series, int dim);
    static std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix);
};

#endif // TIMESERIES_H
