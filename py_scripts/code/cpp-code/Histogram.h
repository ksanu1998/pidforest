// #include "Forest.h"
// #include "Node.h"
// #include "PointSet.h"
// #include "Cube.h"

class Node;
class PointSet;
class Cube;
class Forest;

class Histogram {
private:
    int num;
    int max_buckets;
    std::vector<double> val;
    std::vector<double> count;
    double eps;
    std::vector<double> err;
    std::vector<std::map<int, std::tuple<double, double, double, double, double>>> b_values;

    public:
    Histogram(std::vector<double>& val, std::vector<double>& count, int max_buckets, double eps)
        : num(val.size()), max_buckets(max_buckets), val(val), count(count), eps(eps)
    {
        std::tie(err, b_values) = approx_buckets(val, count, max_buckets, eps);
    }

    std::pair<std::vector<double>, std::vector<std::map<int, std::tuple<double, double, double, double, double>>>> approx_buckets(std::vector<double>& arr, std::vector<double>& counts, int max_buckets, double eps)
    {
        std::vector<double> err_a(max_buckets, -1);
        std::vector<double> cur_err(max_buckets, 0);
        std::vector<std::map<int, std::tuple<double, double, double, double, double>>> b_values(max_buckets);

        double cur_sum = 0;
        double cur_sq = 0;
        double cur_pts = 0;

        for (int j = 0; j < num; j++) {
            cur_sum += arr[j] * counts[j];
            cur_sq += std::pow(arr[j], 2) * counts[j];
            cur_pts += counts[j];

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

                for (auto& b : b_values[k - 1]) {
                    int b_val = b.first;
                    auto [b_a, b_err, b_sum, b_sq, b_pts] = b.second;

                    if (b_val < j) {
                        double tmp_error = b_err + cur_sq - b_sq - std::pow((cur_sum - b_sum), 2) / (cur_pts - b_pts);

                        if (tmp_error < cur_err[k]) {
                            cur_err[k] = tmp_error;
                            a_val = b_val + 1;
                        }
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

    std::string toString() {
        std::string str_val = "";
        for (int i = 0; i < max_buckets; i++) {
            str_val += "Level " + std::to_string(i) + ":";
            for (auto const& [b_val, val_tuple] : b_values[i]) {
                double a_val, err_a;
                std::tie(a_val, err_a, std::ignore, std::ignore, std::ignore) = val_tuple;
                str_val += " (" + std::to_string(a_val) + ", " + std::to_string(b_val) + "): " + std::to_string(err_a) + ", ";
            }
            str_val += "\n";
        }
        return str_val;
    }

    

    void test() {
        std::cout << toString();
        // for (int i = 1; i < max_buckets; i++) {
        //     std::cout << compute_buckets(i) << std::endl;
        // }
        std::cout << "Best Buckets: " << std::endl;
        auto [opt, var_red, buckets] = best_split();
        std::cout << opt << ", " << var_red << ", ";
        for (auto const& bucket : buckets) {
            std::cout << bucket << " ";
        }
        std::cout << std::endl;
    }

    std::tuple<int, double, std::vector<int>> best_split() {
        if (err[0] == 0) {
            return {0, 0, {}};
        }
        std::vector<double> err_red;
        for (int i = 1; i < max_buckets; i++) {
            err_red.push_back(err[0] - err[i]);
        }
        double max_err_red = *std::max_element(err_red.begin(), err_red.end());
        double var_red = max_err_red / err[0];
        if (var_red < 0) {
            std::cout << "error: var_red is " << var_red << std::endl;
            var_red = 0;
        }
        int opt = std::distance(err_red.begin(), std::max_element(err_red.begin(), err_red.end())) + 2;
        auto buckets = compute_buckets(opt);
        return {opt, var_red, buckets};
    }

    std::vector<int> compute_buckets(int num_buckets) {
        std::vector<int> buckets;
        int end = num - 1;
        int k = num_buckets - 1;
        while (end >= 0) {
            int start = std::get<0>(b_values[k][end]);
            if (start <= end) {
                buckets.push_back(start);
            }
            end = start - 1;
            k -= 1;
        }
        std::reverse(buckets.begin(), buckets.end());
      return buckets;
	}
};