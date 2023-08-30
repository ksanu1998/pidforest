#ifndef FOREST_H
#define FOREST_H

#include <unordered_map>
#include <vector>
#include <tuple>
#include <numeric>
#include <random>
class Node; // Forward declaration of the Node class

class Forest {
public:
    int n_trees;
    int max_depth;
    int max_samples;
    int max_buckets;
    double epsilon;
    double sample_axis;
    double threshold;
    // std::vector<std::vector<double>> points;
    std::vector<double> start;
    std::vector<double> end;
    std::vector<Node> tree;
    std::vector<double> n_leaves;

// public:
    int dim;
    int size;
    std::vector<std::vector<double>> points;

    Forest(const std::unordered_map<std::string, double>& kwargs);

    void fit(const std::vector<std::vector<double>>& pts);
    std::tuple<std::vector<int>, std::unordered_map<int, std::vector<double>>,
            std::unordered_map<int, std::vector<double>>, std::unordered_map<int, double>, std::vector<double>>
    predict(const std::vector<std::vector<double>>& pts, double err = 0.1, double pct = 50);
};

#endif // FOREST_H
