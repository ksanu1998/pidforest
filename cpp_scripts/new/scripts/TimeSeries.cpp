#include "TimeSeries.h"

std::vector<std::vector<double>> TimeSeries::shingle(const std::vector<double>& series, int dim) {
    int height = series.size() - dim + 1;
    std::vector<std::vector<double>> shingled(dim, std::vector<double>(height));

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < height; ++j) {
            shingled[i][j] = series[j + i];
        }
    }

    return shingled;
}

std::vector<std::vector<double>> TimeSeries::transpose(const std::vector<std::vector<double>>& matrix) {
    if (matrix.empty())
        return {};

    int rows = matrix.size();
    int cols = matrix[0].size();

    std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}
