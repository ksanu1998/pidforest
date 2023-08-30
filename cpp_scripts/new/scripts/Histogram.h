#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>

class Histogram {
public:
    struct BestSplit {
        int opt;
        double var_red;
        std::vector<int> buckets;
    };

    Histogram(const std::vector<double>& val, const std::vector<int>& count, int max_buckets, double eps);

    std::string toString() const;
    void test() const;
    BestSplit best_split() const;
    std::vector<int> compute_buckets(int num_buckets) const;

private:
        struct BucketValues {
        int a;
        double err_a;
        double cur_sum;
        double cur_sq;
        int cur_pts;

        BucketValues(int a_, double err, double sum, double sq, int pts)
            : a(a_), err_a(err), cur_sum(sum), cur_sq(sq), cur_pts(pts) {}

        // Default constructor for map operations
        BucketValues() : a(0), err_a(0), cur_sum(0), cur_sq(0), cur_pts(0) {}
    };

    int num;
    int max_buckets;
    std::vector<double> val;
    std::vector<int> count;
    double eps;
    std::vector<double> err;
    std::vector<std::unordered_map<int, BucketValues>> b_values;

    std::pair<std::vector<double>, std::vector<std::unordered_map<int, BucketValues>>> approx_buckets(
        const std::vector<double>& arr, const std::vector<int>& count, int max_buckets, double eps) const;
};

#endif  // HISTOGRAM_H
