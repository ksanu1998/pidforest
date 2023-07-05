/*
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>

class Histogram {
private:
    int num;
    int max_buckets;
    std::vector<int> val;
    std::vector<int> count;
    double eps;
    std::vector<double> err;
    std::vector<std::unordered_map<int, std::tuple<int, double, int, int, int>>> b_values;

public:
    Histogram(const std::vector<int>& values, const std::vector<int>& counts, int maxBuckets, double epsilon);

    std::string toString() const;

    struct BestSplit {
        int opt;
        int var_red;
        std::vector<int> buckets;
    };

    BestSplit best_split() const;

private:
    std::vector<int> compute_buckets(int num_buckets) const;

    std::pair<std::vector<double>, std::vector<std::unordered_map<int, std::tuple<int, double, int, int, int>>>> approx_buckets(
        const std::vector<int>& arr, const std::vector<int>& count, int max_buckets, double eps) const;
};

#endif  // HISTOGRAM_H
*/
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <cmath>

class Histogram {
public:
    struct BestSplit {
        int opt;
        double var_red;
        std::vector<int> buckets;
    };

    Histogram(const std::vector<int>& values, const std::vector<int>& counts, int maxBuckets, double epsilon);

    std::string toString() const;
    BestSplit best_split() const;
    std::vector<int> compute_buckets(int num_buckets) const;

private:
    std::pair<std::vector<double>, std::vector<std::unordered_map<int, std::tuple<int, double, int, int, int>>>> approx_buckets(
        const std::vector<int>& arr, const std::vector<int>& count, int max_buckets, double eps) const;

    int num;
    int max_buckets;
    std::vector<int> val;
    std::vector<int> count;
    double eps;
    std::vector<double> err;
    std::vector<std::unordered_map<int, std::tuple<int, double, int, int, int>>> b_values;
};

#endif  // HISTOGRAM_H

