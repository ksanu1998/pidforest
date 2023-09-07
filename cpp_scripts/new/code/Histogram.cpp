#include "Histogram.h"

Histogram::Histogram(const std::vector<double>& values, const std::vector<int>& counts, int maxBuckets, double epsilon)
    : num(values.size()), max_buckets(maxBuckets), val(values), count(counts), eps(epsilon) {
    std::tie(err, b_values) = approx_buckets(val, count, max_buckets, eps);
    /*
    std::cout << "err, [";
    for (int e = 0; e < err.size(); e++) {
            std::cout << err[e] << " ";
        }
    std::cout << "]" << std::endl;
    std::cout << "num, " << num << std::endl;
    // Print the elements of b_values for each level
    for (int k = 0; k < b_values.size(); k++) {
        std::cout << "b_values[" << k << "]: ";
        for (const auto& entry : b_values[k]) {
            int a = entry.second.a;
            double err_a = entry.second.err_a;
            double cur_sum = entry.second.cur_sum;
            double cur_sq = entry.second.cur_sq;
            int cur_pts = entry.second.cur_pts;
            std::cout << "(" << entry.first << ": " << a << ", " << err_a << ", "
                      << cur_sum << ", " << cur_sq << ", " << cur_pts << ") ";
        }
        std::cout << std::endl;
    }
    */
}

std::string Histogram::toString() const {
    std::string str_val;
    for (int i = 0; i < max_buckets; i++) {
        str_val += "Level " + std::to_string(i) + ":";
        for (const auto& b_val : b_values[i]) {
            // int a_val = b_val.second.a_val;
            int a = b_val.second.a;
            double err_a = b_val.second.err_a;
            str_val += " (" + std::to_string(a) + ", " + std::to_string(err_a) + "), ";
        }
        str_val += "\n";
    }
    return str_val;
}

void Histogram::test() const {
    std::cout << toString();
    /*
    for (int i = 1; i < max_buckets; i++) {
        std::vector<int> buckets = compute_buckets(i);
        std::cout << "Buckets for " << i << " buckets: ";
        for (int bucket : buckets) {
            std::cout << bucket << " ";
        }
        std::cout << std::endl;
    }
    */
    BestSplit bestsplit = best_split();
    /*
    std::cout << "Best Split: opt = " << bestsplit.opt << ", var_red = " << bestsplit.var_red << ", buckets = ";
    for (int bucket : bestsplit.buckets) {
        std::cout << bucket << " ";
    }
    std::cout << std::endl;
    */
}

Histogram::BestSplit Histogram::best_split() const {
    // std::cout << "***** Histogram::best_split() *****" << std::endl;
    if (err[0] == 0) {
        BestSplit bestsplit;
        bestsplit.opt = 0;
        bestsplit.var_red = 0;
        bestsplit.buckets = std::vector<int>();
        // std::cout << "NULL" << std::endl;
        return bestsplit;
    }
    // std::cout << "NOT NULL" << std::endl;
    std::vector<double> err_red(max_buckets - 1);
    
    // std::cout << "err_red, var_red: [";
    for (int i = 1; i < max_buckets; i++) {
        // std::cout << "err[0], err[i], err[0] - err[i]: " << err[0] << ", " << err[i] << ", " << err[0] - err[i] << std::endl;
        err_red[i - 1] = err[0] - err[i];
        // std::cout << err_red[i - 1] << " ";
    }
    // std::cout << "], ";
    
    double max_err_red = *std::max_element(err_red.begin(), err_red.end());
    // std::cout << "max_err_red, " << max_err_red << std::endl;
    double var_red = max_err_red / err[0];
    // std::cout << "var_red, " << var_red << std::endl;
    
    // std::cout <<  var_red << std::endl;
    if (var_red < 0) {
        var_red = 0;
        // std::cout << "var_red < 0" << std::endl;
    }
    
    int opt = std::distance(err_red.begin(), std::max_element(err_red.begin(), err_red.end())) + 2;
    // std::cout << "opt, " << opt << std::endl;
    // std::cout << "compute_buckets" << std::endl;
    std::vector<int> buckets = compute_buckets(opt);
    // std::cout << "compute_buckets [DONE]" << std::endl;
    BestSplit bestsplit;
    bestsplit.opt = opt;
    bestsplit.var_red = var_red;
    bestsplit.buckets = buckets;
    // std::cout << "var_red, " << var_red << std::endl;
    // exit(1);
    /*
    std::cout << "buckets, [";
    for (int b = 0; b < buckets.size(); b++) {
        std::cout << buckets[b] << " ";
    }
    std::cout << "]" << std::endl;
    */
    return bestsplit;
    // std::cout << "***** Histogram::best_split() [DONE]*****" << std::endl;
    
}

// std::vector<int> Histogram::compute_buckets(int num_buckets) const {
//     std::vector<int> buckets;
//     int end = num - 1;
//     int k = num_buckets - 1;
//     std::cout << "compute_buckets.while" << std::endl;
//     while (end >= 0) {
//         int start = b_values[k].at(end).a;
//         if (start <= end) {
//             buckets.push_back(start);
//         }
//         end = start - 1;
//         k -= 1;
//     }
//     std::cout << "compute_buckets.while [DONE]" << std::endl;
//     std::cout << "compute_buckets.reverse" << std::endl;
//     std::reverse(buckets.begin(), buckets.end());
//     std::cout << "compute_buckets.reverse [DONE]" << std::endl;
//     return buckets;
// }

std::vector<int> Histogram::compute_buckets(int num_buckets) const {
    std::vector<int> buckets;
    int end = num - 1;
    int k = num_buckets - 1;
    // std::cout << "compute_buckets.while" << std::endl;
    while (end >= 0) {
        if (b_values[k].find(end) != b_values[k].end()) {
            int start = b_values[k].at(end).a;
            if (start <= end) {
                buckets.push_back(start);
            }
            end = start - 1;
            k -= 1;
        } else {
            // Handle the case when key 'end' is not found in the map
            // This might involve an error in the algorithm or data processing
            // Print relevant debug information to identify the issue
            // std::cerr << "Key " << end << " not found in b_values[" << k << "]" << std::endl;
            break; // Exit the loop to prevent an infinite loop
        }
    }
    // std::cout << "compute_buckets.while [DONE]" << std::endl;
    // std::cout << "compute_buckets.reverse" << std::endl;
    std::reverse(buckets.begin(), buckets.end());
    // std::cout << "compute_buckets.reverse [DONE]" << std::endl;
    return buckets;
}


/*
std::pair<std::vector<double>, std::vector<std::unordered_map<int, Histogram::BucketValues>>> Histogram::approx_buckets(
    const std::vector<double>& arr, const std::vector<int>& count, int max_buckets, double eps) const {
    std::vector<double> err_a(max_buckets, -1);
    std::vector<double> cur_err(max_buckets);
    std::vector<std::unordered_map<int, BucketValues>> b_values(max_buckets);

    // std::cout << "arr: ";
    // for (const double& val : arr) {
    //     std::cout << val << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "count: ";
    // for (const int& c : count) {
    //     std::cout << c << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "max_buckets: " << max_buckets << std::endl;
    // std::cout << "eps: " << eps << std::endl;

    // Initialize cur_sum, cur_sq, cur_pts
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
        b_values[0][j] = BucketValues(0, cur_err[0], cur_sum, cur_sq, cur_pts);

        for (int k = 1; k < max_buckets; k++) {
            cur_err[k] = cur_err[k - 1];
            int a_val = j + 1;
            for (const auto& b_val : b_values[k - 1]) {
                if (b_val.first >= j) {
                    break;
                }
                const BucketValues& b_data = b_val.second;
                double b_err = b_data.err_a;
                double b_sum = b_data.cur_sum;
                double b_sq = b_data.cur_sq;
                int b_pts = b_data.cur_pts;
                double tmp_error = b_err + cur_sq - b_sq - std::pow((cur_sum - b_sum), 2) / (cur_pts - b_pts);
                if (tmp_error < cur_err[k]) {
                    cur_err[k] = tmp_error;
                    a_val = b_val.first + 1;
                }
            }
            b_values[k][j] = BucketValues(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
            if (cur_err[k] > (1 + eps) * err_a[k]) {
                err_a[k] = cur_err[k];
            } else {
                b_values[k].erase(j - 1);
            }
        }

        // for (int k = 1; k < max_buckets; k++) {
        //     cur_err[k] = cur_err[k - 1];
        //     int a_val = j + 1;
        //     for (const auto& b_val : b_values[k - 1]) {
        //         if (b_val.first >= j) {
        //             break;
        //         }
        //         double b_err = b_val.second.err;
        //         double b_sum = b_val.second.sum;
        //         double b_sq = b_val.second.sq;
        //         int b_pts = b_val.second.pts;
        //         double tmp_error = b_err + cur_sq - b_sq - std::pow((cur_sum - b_sum), 2) / (cur_pts - b_pts);
        //         if (tmp_error < cur_err[k]) {
        //             cur_err[k] = tmp_error;
        //             a_val = b_val.first + 1;
        //         }
        //     }
        //     b_values[k][j] = BucketValues(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
        //     if (cur_err[k] > (1 + eps) * err_a[k]) {
        //         err_a[k] = cur_err[k];
        //     } else {
        //         b_values[k].erase(j - 1);
        //     }
        // }
    }
    std::cout << "cur_err, [";
    for (int e = 0; e < cur_err.size(); e++) {
        std::cout << cur_err[e] << " ";
    }
    std::cout << "]" << std::endl;
    return std::make_pair(cur_err, b_values);
}
*/
/*
std::pair<std::vector<double>, std::vector<std::unordered_map<int, Histogram::BucketValues>>> Histogram::approx_buckets(
    const std::vector<double>& arr, const std::vector<int>& count, int max_buckets, double eps) const {
    std::vector<double> err_a(max_buckets, -1);
    std::vector<double> cur_err(max_buckets);
    std::vector<std::unordered_map<int, BucketValues>> b_values(max_buckets);

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
            b_values[0][j] = BucketValues(0, cur_err[0], cur_sum, cur_sq, cur_pts);
        }

        for (int k = 1; k < max_buckets; k++) {
            cur_err[k] = cur_err[k - 1];
            int a_val = j + 1;
            for (const auto& b_val : b_values[k - 1]) {
                if (b_val.first >= j) {
                    break;
                }
                const BucketValues& b_data = b_val.second;
                double b_err = b_data.err_a;
                double b_sum = b_data.cur_sum;
                double b_sq = b_data.cur_sq;
                int b_pts = b_data.cur_pts;
                
                double numerator = cur_sq - b_sq - std::pow(cur_sum - b_sum, 2) / (cur_pts - b_pts);
                double denominator = cur_pts - b_pts;
                double tmp_error = b_err + (numerator / denominator);

                if (tmp_error < cur_err[k] && tmp_error > err_a[k] * (1 - eps)) {
                    cur_err[k] = tmp_error;
                    a_val = b_val.first + 1;
                }
            }
            if (cur_err[k] > (1 + eps) * err_a[k]) {
                err_a[k] = cur_err[k];
                b_values[k][j] = BucketValues(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
            } else {
                b_values[k].erase(j);
            }
        }
    }
    return std::make_pair(cur_err, b_values);
}
*/
/*
std::pair<std::vector<double>, std::vector<std::unordered_map<int, Histogram::BucketValues>>> Histogram::approx_buckets(
    const std::vector<double>& arr, const std::vector<int>& count, int max_buckets, double eps) const {
    std::vector<double> err_a(max_buckets, -1);
    std::vector<double> cur_err(max_buckets);
    std::vector<std::unordered_map<int, BucketValues>> b_values(max_buckets);

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
            b_values[0][j] = BucketValues(0, cur_err[0], cur_sum, cur_sq, cur_pts);
        }

        for (int k = 1; k < max_buckets; k++) {
            cur_err[k] = cur_err[k - 1];
            int a_val = j + 1;
            for (const auto& b_val : b_values[k - 1]) {
                if (b_val.first >= j) {
                    break;
                }
                const BucketValues& b_data = b_val.second;
                double b_err = b_data.err_a;
                double b_sum = b_data.cur_sum;
                double b_sq = b_data.cur_sq;
                int b_pts = b_data.cur_pts;
                double tmp_error = b_err + cur_sq - b_sq - std::pow((cur_sum - b_sum), 2) / (cur_pts - b_pts);
                if (tmp_error < cur_err[k]) {
                    cur_err[k] = tmp_error;
                    a_val = b_val.first + 1;
                }
            }
            // if (cur_err[k] > (1 + eps) * err_a[k]) {
            //     err_a[k] = cur_err[k];
            //     b_values[k][j] = BucketValues(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
            // } else {
            //     b_values[k].erase(j);
            // }

            b_values[k][j] = BucketValues(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
            if (cur_err[k] > (1 + eps) * err_a[k]) {
                err_a[k] = cur_err[k];
            } else {
                b_values[k].erase(j - 1);
            }

        }
    }
    std::cout << "cur_err, [";
    for (int e = 0; e < cur_err.size(); e++) {
        std::cout << cur_err[e] << " ";
    }
    std::cout << "]" << std::endl;
    return std::make_pair(cur_err, b_values);
}
*/
/*
std::pair<std::vector<double>, std::vector<std::unordered_map<int, Histogram::BucketValues>>> Histogram::approx_buckets(
    const std::vector<double>& arr, const std::vector<int>& count, int max_buckets, double eps) const {
    
    std::vector<double> err_a(max_buckets, -1);
    std::vector<double> cur_err(max_buckets);
    std::vector<std::unordered_map<int, BucketValues>> b_values(max_buckets);
    
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
        b_values[0][j] = BucketValues(0, cur_err[0], cur_sum, cur_sq, cur_pts);
        
        for (int k = 1; k < max_buckets; k++) {
            cur_err[k] = cur_err[k - 1];
            int a_val = j + 1;
            for (const auto& b_val : b_values[k - 1]) {
                if (b_val.first >= j) {
                    break;
                }
                const BucketValues& b_data = b_val.second;
                double b_err = b_data.err_a;
                double b_sum = b_data.cur_sum;
                double b_sq = b_data.cur_sq;
                int b_pts = b_data.cur_pts;
                double tmp_error = b_err + cur_sq - b_sq - std::pow((cur_sum - b_sum), 2) / (cur_pts - b_pts);
                if (tmp_error < cur_err[k]) {
                    cur_err[k] = tmp_error;
                    a_val = b_val.second.a + 1;  // Use the 'a' value from the previous map
                }
            }
            b_values[k][j] = BucketValues(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
            if (cur_err[k] > (1 + eps) * err_a[k]) {
                err_a[k] = cur_err[k];
            } else {
                b_values[k].erase(j - 1);
            }
        }
    }
    std::cout << "cur_err, [";
    for (int e = 0; e < cur_err.size(); e++) {
        std::cout << cur_err[e] << " ";
    }
    std::cout << "]" << std::endl;
    return std::make_pair(cur_err, b_values);
    return std::make_pair(cur_err, b_values);
}
*/
std::pair<std::vector<double>, std::vector<std::unordered_map<int, Histogram::BucketValues>>> Histogram::approx_buckets(
    const std::vector<double>& arr, const std::vector<int>& count, int max_buckets, double eps) const {
    
    std::vector<double> err_a(max_buckets, -1);
    std::vector<double> cur_err(max_buckets);
    std::vector<std::unordered_map<int, BucketValues>> b_values(max_buckets);
    /*
    std::cout << "arr: [";
    for (int a = 0; a < arr.size(); a++) {
        std::cout << arr[a] << " ";
    }
    std::cout << "]" << std::endl;
    */
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
        b_values[0][j] = BucketValues(0, cur_err[0], cur_sum, cur_sq, cur_pts);
        
        for (int k = 1; k < max_buckets; k++) {
            cur_err[k] = cur_err[k - 1];
            int a_val = j + 1;
            for (const auto& b_val_data : b_values[k - 1]) {
                int b_val = b_val_data.first;
                if (b_val >= j) {
                    break;
                }
                const BucketValues& b_data = b_val_data.second;
                double b_err = b_data.err_a;
                double b_sum = b_data.cur_sum;
                double b_sq = b_data.cur_sq;
                int b_pts = b_data.cur_pts;
                double tmp_error = b_err + cur_sq - b_sq - std::pow((cur_sum - b_sum), 2) / (cur_pts - b_pts);
                if (tmp_error < cur_err[k]) {
                    cur_err[k] = tmp_error;
                    // a_val = b_data.a + 1;  // Use the 'a' value from the previous map
                    a_val = b_val + 1;
                }
            }
            b_values[k][j] = BucketValues(a_val, cur_err[k], cur_sum, cur_sq, cur_pts);
            if (cur_err[k] > (1 + eps) * err_a[k]) {
                err_a[k] = cur_err[k];
            } else {
                b_values[k].erase(j - 1);
            }
        }
    }
    return std::make_pair(cur_err, b_values);
}


