#include "Histogram.h"

Histogram::Histogram(const std::vector<int>& values, const std::vector<int>& counts, int maxBuckets, double epsilon)
    : num(values.size()), max_buckets(maxBuckets), val(values), count(counts), eps(epsilon) {
    std::tie(err, b_values) = approx_buckets(val, count, max_buckets, eps);
}

std::string Histogram::toString() const {
    std::string str_val;
    for (int i = 0; i < max_buckets; i++) {
        str_val += "Level " + std::to_string(i) + ":";
        for (auto& b_val : b_values[i]) {
            int a_val = std::get<0>(b_val.second);
            double err_a = std::get<1>(b_val.second);
            str_val += " (" + std::to_string(a_val) + ", " + std::to_string(err_a) + "), ";
        }
        str_val += "\n";
    }
    return str_val;
}

Histogram::BestSplit Histogram::best_split() const {
    if (err[0] == 0) {
        BestSplit bestsplit;
        bestsplit.opt = 0;
        bestsplit.var_red = 0;
        bestsplit.buckets = std::vector<int>();
        return bestsplit;
    }
    std::vector<double> err_red(max_buckets - 1);
    for (int i = 1; i < max_buckets; i++) {
        err_red[i - 1] = err[0] - err[i];
    }
    double max_err_red = *std::max_element(err_red.begin(), err_red.end());
    double var_red = max_err_red / err[0];
    if (var_red < 0) {
        std::cout << "error: var_red is " << var_red << std::endl;
        var_red = 0;
    }
    int opt = std::distance(err_red.begin(), std::max_element(err_red.begin(), err_red.end())) + 2;
    std::vector<int> buckets = compute_buckets(opt);
    BestSplit bestsplit;
    bestsplit.opt = opt;
    bestsplit.var_red = var_red;
    bestsplit.buckets = buckets;
    return bestsplit;
}

std::vector<int> Histogram::compute_buckets(int num_buckets) const {
    std::vector<int> buckets;
    int end = num - 1;
    int k = num_buckets - 1;
    while (end >= 0) {
        int start = std::get<0>(b_values[k].at(end));
        if (start <= end) {
            buckets.push_back(start);
        }
        end = start - 1;
        k -= 1;
    }
    std::reverse(buckets.begin(), buckets.end());
    return buckets;
}

std::pair<std::vector<double>, std::vector<std::unordered_map<int, std::tuple<int, double, int, int, int>>>> Histogram::approx_buckets(
    const std::vector<int>& arr, const std::vector<int>& count, int max_buckets, double eps) const {
    std::vector<double> err_a(max_buckets, -1);
    std::vector<double> cur_err(max_buckets);
    std::vector<std::unordered_map<int, std::tuple<int, double, int, int, int>>> b_values(max_buckets);
    double cur_sum = 0;
    double cur_sq = 0;
    int cur_pts = 0;
    for (int j = 0; j < arr.size(); j++) {
        cur_sum += arr[j] * count[j];
        cur_sq += std::pow(arr[j], 2) * count[j];
        cur_pts += count[j];
        cur_err[0] = cur_sq - std::pow(cur_sum, 2) / cur_pts;
        if (cur_err[0] > (1 + eps) * err_a[0]) {
            err_a[0] = cur_err[0];
        } else {
            b_values[0].erase(j - 1);
        }
        b_values[0][j] = std::make_tuple(0, cur_err[0], cur_sum, cur_sq, cur_pts);
        for (int k = 1; k < max_buckets; k++) {
            cur_err[k] = cur_err[k - 1];
            int a_val = j + 1;
            for (const auto& b_val : b_values[k - 1]) {
                if (b_val.first >= j) {
                    break;
                }
                double b_err = std::get<1>(b_val.second);
                double b_sum = std::get<2>(b_val.second);
                double b_sq = std::get<3>(b_val.second);
                int b_pts = std::get<4>(b_val.second);
                double tmp_error = b_err + cur_sq - b_sq - std::pow((cur_sum - b_sum), 2) / (cur_pts - b_pts);
                if (tmp_error < cur_err[k]) {
                    cur_err[k] = tmp_error;
                    a_val = b_val.first + 1;
                }
            }
            b_values[k][j] = std::make_tuple(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
            if (cur_err[k] > (1 + eps) * err_a[k]) {
                err_a[k] = cur_err[k];
            } else {
                b_values[k].erase(j - 1);
            }
        }
    }
    return std::make_pair(cur_err, b_values);
}
