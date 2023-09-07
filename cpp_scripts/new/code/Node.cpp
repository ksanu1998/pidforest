// Node.cpp

#include "Node.h"
#include "Cube.h"
#include "PointSet.h"
#include "Histogram.h"
#include <random>
#include <algorithm>

Node::Node(int depth, Forest* forest, const std::unordered_map<std::string, std::variant<std::vector<int>, std::vector<double>>>& kwargs)
    : depth(depth), forest(forest), 
      point_set(this, std::unordered_set<int>()), 
      cube(this, int(), std::vector<double>(), std::vector<double>())
{
    if (depth == 0) {
        id_string = {0};
        cube = Cube(this, forest->dim, forest->start, forest->end);
        point_set = PointSet(this, std::unordered_set<int>(std::get<std::vector<int>>(kwargs.at("indices")).begin(), std::get<std::vector<int>>(kwargs.at("indices")).end()));
    }
    else {
        // std::cout << ">>>>>>>>>>>>>>>>Nodehere" << std::endl;

        id_string = std::get<std::vector<int>>(kwargs.at("id"));
        std::vector<double> startArray = std::get<std::vector<double>>(kwargs.at("start"));
        std::vector<double> endArray = std::get<std::vector<double>>(kwargs.at("end"));
        std::vector<double> startVector(startArray.begin(), startArray.end());
        std::vector<double> endVector(endArray.begin(), endArray.end());
        cube = Cube(this, forest->dim, startVector, endVector);
        std::vector<int> filtered_indices = cube.filter_indices(std::vector<int>(std::get<std::vector<int>>(kwargs.at("indices")).begin(), std::get<std::vector<int>>(kwargs.at("indices")).end()));
        std::unordered_set<int> filtered_indices_set(filtered_indices.begin(), filtered_indices.end());
        // std::cout << ">>>>>>>>>>>>>>>>Nodehere[BET]" << std::endl;
        point_set = PointSet(this, filtered_indices_set); 
        // point_set = PointSet(this, std::unordered_set<int>(std::get<std::vector<int>>(kwargs.at("indices")).begin(), std::get<std::vector<int>>(kwargs.at("indices")).end())); 
        // std::cout << ">>>>>>>>>>>>>>>>Nodehere[DONE]" << std::endl;
    }
    density = -1;
    child = std::vector<Node>();
    if ((depth < forest->max_depth) && (point_set.indices.size() > 1)) {
        // std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>depth, point_set.indices.size(): " << depth << ", " << point_set.indices.size() << std::endl;
        find_split();
    }    
}

void Node::find_split()
{
    // std::cout << "cube.dim, " << cube.dim << std::endl;
    std::vector<int> imp_axis;
    for (int axis = 0; axis < cube.dim; ++axis) {
        if (point_set.val[axis].size() > 1) {
            imp_axis.push_back(axis);
        }
    }
    /*
    std::cout << "imp_axis, [";
    for (int i = 0; i < imp_axis.size(); i++) {
        if(i < imp_axis.size() - 1) {
            std::cout << imp_axis[i] << ", ";
        }
        else {
            std::cout << imp_axis[i];
        }
    } 
    */
    // std::cout << "]" << std::endl;
    if (imp_axis.empty()) {
        // std::cout << "if not imp_axis" << std::endl;
        return;
    }

    int max_axes = std::min(static_cast<int>(imp_axis.size()), static_cast<int>(forest->sample_axis * cube.dim));
    // std::cout << "max_axes, " << max_axes << std::endl;
    std::vector<int> s_axes;
    s_axes.reserve(imp_axis.size());  // Reserve space for all axes

    // Copy the elements of imp_axis to s_axes
    s_axes.insert(s_axes.end(), imp_axis.begin(), imp_axis.end());

    // Shuffle the elements in s_axes
    std::shuffle(s_axes.begin(), s_axes.end(), std::mt19937{std::random_device{}()});
    s_axes.resize(max_axes);  // Keep only the first max_axes elements
    /*
    std::cout << "s_axes, [";
    for (int i = 0; i < s_axes.size(); i++) {
        if(i < s_axes.size() - 1) {
            std::cout << s_axes[i] << ", ";
        }
        else {
            std::cout << s_axes[i];
        }
    } 
    std::cout << "]" << std::endl;
    */
    std::unordered_map<int, std::vector<int>> buckets;
    std::unordered_map<int, double> var_red;
    
    for (int axis : s_axes) {
        std::vector<double> int_gap(point_set.gap[axis].size());
        for (std::size_t i = 0; i < point_set.gap[axis].size(); ++i) {
            int_gap[i] = static_cast<double>(point_set.gap[axis][i]) / point_set.count[axis][i];
        }
        /*
        std::cout << "point_set.gap[axis], [";
        for (int e = 0; e < int_gap.size(); e++) {
                std::cout << int_gap[e] << " ";
        }
        */
        std::vector<int> int_count(point_set.count[axis].begin(), point_set.count[axis].end());
        // std::cout << "int_gap.size(), " << int_gap.size() << std::endl;
        Histogram hist(int_gap, int_count, forest->max_buckets, forest->epsilon);
        Histogram::BestSplit best_split = hist.best_split();
        var_red[axis] = best_split.var_red;
        const std::vector<int>& bucketValues = best_split.buckets;
        buckets[axis] = bucketValues; // Just assign the vector, no need to copy element by element
        
        var_red[axis] = best_split.var_red;
        buckets[axis] = bucketValues;
        /*
        std::cout << "axis, var_red[axis], buckets[axis]: " << axis << ", " << var_red[axis] << ", [";
        for (int b : buckets[axis]) {
            std::cout << b << ", ";
        }
        std::cout << "]" << std::endl;
        */

    }
    
    double max_var_red = std::max_element(var_red.begin(), var_red.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
    })->second;

    if (max_var_red <= forest->threshold) {
        return;
    }
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
        // std::cout << "cube.split_vals.size() + 1 :  " << cube.split_vals.size() + 1 << std::endl;
        std::vector<double> new_start = cube.start; // = {cube.start[0], cube.start[1], cube.start[2]};
        // for (std::size_t i = 0; i < cube.start.size(); ++i) {
        //     new_start[i] = cube.start[i];
        // }
        // std::cout << "new_start.size: " << new_start.size() << std::endl;
        std::vector<double> new_end = cube.end; // = {cube.end[0], cube.end[1], cube.end[2]};
        // for (std::size_t i = 0; i < cube.end.size(); ++i) {
        //     new_end[i] = cube.end[i];
        // }
        // std::cout << "new_end.size: " << new_end.size() << std::endl;
        // std::cout << ">>>>>>>>>>>>>>>>here" << std::endl;
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
        // std::cout << ">>>>>>>>>>>>>>>>here[DONE]" << std::endl;
        

        std::unordered_map<std::string, std::variant<std::vector<int>, std::vector<double>>> kwargs;
        kwargs["start"] = new_start;
        kwargs["end"] = new_end;
        kwargs["indices"] = std::vector<int>(point_set.indices.begin(), point_set.indices.end());
        kwargs["id"] = id_string;
        
        Node child_node(depth + 1, forest, kwargs);
        child.push_back(child_node);
        cube.child.push_back(&(child.back().cube));
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
    std::cout << "density, " << density << std::endl;
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

