// Node.cpp

#include "Node.h"
#include "Cube.h"
#include "PointSet.h"
#include "Histogram.h"
#include <random>
#include <algorithm>

Node::Node(int depth, Forest* forest, const std::unordered_map<std::string, std::variant<std::vector<int>, std::array<double, 3>>>& kwargs)
    : depth(depth), forest(forest),
      point_set(this, std::unordered_set<int>()), 
      cube(this, std::vector<double>(), std::vector<double>())
{
    if (depth == 0) {
        id_string = {0};
        cube = Cube(this, forest->start, forest->end);
        // std::cout << "FINISHED BUILDING CUBE" << std::endl;
        point_set = PointSet(this, std::unordered_set<int>(std::get<std::vector<int>>(kwargs.at("indices")).begin(), std::get<std::vector<int>>(kwargs.at("indices")).end()));
        // std::cout << "FINISHED BUILDING POINTSET" << std::endl;
    }
    else {
        id_string = std::get<std::vector<int>>(kwargs.at("id"));
        std::array<double, 3> startArray = std::get<std::array<double, 3>>(kwargs.at("start"));
        std::array<double, 3> endArray = std::get<std::array<double, 3>>(kwargs.at("end"));
        std::vector<double> startVector(startArray.begin(), startArray.end());
        std::vector<double> endVector(endArray.begin(), endArray.end());
        cube = Cube(this, startVector, endVector);
        // point_set = PointSet(this, cube.filter_indices(std::get<std::vector<int>>(kwargs.at("indices"))));
        point_set = PointSet(this, std::unordered_set<int>(std::get<std::vector<int>>(kwargs.at("indices")).begin(), std::get<std::vector<int>>(kwargs.at("indices")).end()));
    }
    density = -1;
    child = std::vector<Node>();
    if ((depth < forest->max_depth) && (point_set.indices.size() > 1)) {
        find_split();
    }
}

void Node::find_split()
{
    std::vector<int> imp_axis;
    for (int axis = 0; axis < cube.dim; ++axis) {
        if (point_set.val[axis].size() > 1) {
            imp_axis.push_back(axis);
        }
    }

    if (imp_axis.empty()) {
        return;
    }

    int max_axes = std::min(static_cast<int>(imp_axis.size()), static_cast<int>(forest->sample_axis * cube.dim));
    std::vector<int> s_axes;
    s_axes.reserve(max_axes);

    std::sample(imp_axis.begin(), imp_axis.end(), std::back_inserter(s_axes), max_axes, std::mt19937{std::random_device{}()});

    std::unordered_map<int, std::vector<int>> buckets;
    std::unordered_map<int, double> var_red;

    for (int axis : s_axes) {
        // Histogram hist(point_set.gap[axis] / point_set.count[axis], point_set.count[axis], forest->max_buckets, forest->epsilon);
        std::vector<int> int_gap(point_set.gap[axis].begin(), point_set.gap[axis].end());
        std::vector<int> int_count(point_set.count[axis].begin(), point_set.count[axis].end());
        Histogram hist(int_gap, int_count, forest->max_buckets, forest->epsilon);// hist.best_split();
        // var_red[axis] = hist.var_red;
        // buckets[axis] = hist.buckets;
        Histogram::BestSplit best_split = hist.best_split();
        var_red[axis] = best_split.var_red;
        buckets[axis] = best_split.buckets;
    }

    double max_var_red = std::max_element(var_red.begin(), var_red.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    })->second;

    if (max_var_red <= forest->threshold) {
        return;
    }

    // int split_axis = std::sample(s_axes.begin(), s_axes.end(), std::initializer_list<double>{var_red.begin()->second / std::accumulate(var_red.begin(), var_red.end(), 0.0, [](double sum, const auto& pair) {
    //                                  return sum + pair.second;
    //                              })})
    //                      .front();
    std::vector<int> sample_indices(s_axes.size());
    std::iota(sample_indices.begin(), sample_indices.end(), 0);
    std::shuffle(sample_indices.begin(), sample_indices.end(), std::mt19937{std::random_device{}()});
    int split_axis = s_axes[sample_indices[0]];

    cube.split_axis = split_axis;
    cube.split_vals.resize(buckets[split_axis].size());
    for (std::size_t i = 0; i < buckets[split_axis].size(); ++i) {
        cube.split_vals[i] = (point_set.val[split_axis][buckets[split_axis][i] - 1] + point_set.val[split_axis][buckets[split_axis][i]]) / 2;
    }

    for (std::size_t i = 0; i < cube.split_vals.size() + 1; ++i) {
        // std::array<double, 3> new_start = cube.start;
        std::array<double, 3> new_start = {cube.start[0], cube.start[1], cube.start[2]};
        // std::array<double, 3> new_end = cube.end;
        std::array<double, 3> new_end = {cube.end[0], cube.end[1], cube.end[2]};


        if ((i > 0) && (i < cube.split_vals.size())) {
            new_start[split_axis] = cube.split_vals[i - 1];
            new_end[split_axis] = cube.split_vals[i];
        }
        else if (i == 0) {
            new_end[split_axis] = cube.split_vals[0];
        }
        else {  // i == cube.split_vals.size()
            new_start[split_axis] = cube.split_vals.back();
        }

        std::unordered_map<std::string, std::variant<std::vector<int>, std::array<double, 3>>> kwargs;
        kwargs["start"] = new_start;
        kwargs["end"] = new_end;
        // kwargs["indices"] = point_set.indices;
        kwargs["indices"] = std::vector<int>(point_set.indices.begin(), point_set.indices.end());
        kwargs["id"] = id_string;
        child.emplace_back(depth + 1, forest, kwargs);
        cube.child.push_back(&child.back().cube);
    }
}

void Node::compute_density(const std::vector<int>& indices)
{
    int num = indices.size();
    if (num == 0) {
        density = 0;
        child.clear();
        cube.child.clear();
        cube.split_axis = -1;
        return;
    }
    density = log(num) - cube.vol;
    if (!child.empty()) {
        std::vector<std::vector<int>> index_split = cube.split_indices(forest->points, indices);
        for (std::size_t i = 0; i < child.size(); ++i) {
            child[i].compute_density(index_split[i]);
        }
    }
}

int Node::compute_leaf_num()
{
    if (!child.empty()) {
        int leaves = 0;
        for (std::size_t i = 0; i < child.size(); ++i) {
            leaves += child[i].compute_leaf_num();
        }
        return leaves;
    }
    else {
        return 1;
    }
}

void Node::compute_split(const std::vector<std::vector<double>>& pts, const std::vector<int>& indices, std::vector<double>& scores)
{
    if (!child.empty()) {
        std::vector<std::vector<int>> index_split = cube.split_indices(pts, indices);
        for (std::size_t i = 0; i < child.size(); ++i) {
            if (!index_split[i].empty()) {
                child[i].compute_split(pts, index_split[i], scores);
            }
        }
    }
    else {
        for (int index : indices) {
            scores[index] = density;
        }
    }
}

