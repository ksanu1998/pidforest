#include "Cube.h"
#include "Node.h"
#include "Forest.h"

Cube::Cube(Node* n, int d, const std::vector<double>& s, const std::vector<double>& e)
    : node(n), dim(d), start(s), end(e) {
    // dim = s.size();
    split_axis = -1;
    vol = 0;
    for (int i = 0; i < dim; ++i) {
        vol += std::log(end[i] - start[i]);
    }
}

std::vector<int> Cube::filter_indices(const std::vector<int>& indices) {
    std::vector<int> filtered_indices;
    for (int i : indices) {
        bool in_lb = true;
        bool in_ub = true;
        for (int j = 0; j < dim; ++j) {
            if (node->forest->points[j][i] < start[j] || node->forest->points[j][i] >= end[j]) {
                in_lb = false;
                in_ub = false;
                break;
            }
        }
        if (in_lb && in_ub) {
            filtered_indices.push_back(i);
        }
    }
    return filtered_indices;
}

/*
std::vector<std::vector<int>> Cube::split_indices(const std::vector<std::vector<double>>& pts, const std::vector<int>& indices) {
    int n_child = child.size();
    if (n_child == 0) {
        return { indices };
    }

    int n_arr = indices.size();
    if (n_arr == 0) {
        return std::vector<std::vector<int>>(n_child);
    }

    //std::vector<int> s_arr = pts[split_axis];
    std::vector<int> s_arr(pts[split_axis].begin(), pts[split_axis].end());
    double s_start = start[split_axis];
    double s_end = end[split_axis];

    std::vector<std::vector<int>> index_split(n_child);
    index_split[0].reserve(n_arr);
    index_split[n_child - 1].reserve(n_arr);
    for (int i : indices) {
        if (s_arr[i] >= s_start && s_arr[i] < split_vals[0]) {
            index_split[0].push_back(i);
        }
        else if (s_arr[i] >= split_vals.back() && s_arr[i] < s_end) {
            index_split[n_child - 1].push_back(i);
        }
        else {
            for (int k = 1; k < n_child - 1; ++k) {
                if (s_arr[i] >= split_vals[k - 1] && s_arr[i] < split_vals[k]) {
                    index_split[k].push_back(i);
                    break;
                }
            }
        }
    }
    return index_split;
}
*/
std::vector<std::vector<int>> Cube::split_indices(const std::vector<std::vector<double>>& pts, const std::vector<int>& indices) {
    int n_child = child.size();
    if (n_child == 0) {
        return { indices };
    }

    int n_arr = indices.size();
    if (n_arr == 0) {
        return std::vector<std::vector<int>>(n_child);
    }

    //std::vector<int> s_arr = pts[split_axis];
    std::vector<int> s_arr(pts[split_axis].begin(), pts[split_axis].end());
    double s_start = start[split_axis];
    double s_end = end[split_axis];

    std::vector<std::vector<int>> index_split(n_child);
    index_split[0].reserve(n_arr);
    index_split[n_child - 1].reserve(n_arr);
    for (int i : indices) {
        if (s_arr[i] >= s_start && s_arr[i] < split_vals[0]) {
            index_split[0].push_back(i);
        }
        else if (s_arr[i] >= split_vals.back() && s_arr[i] < s_end) {
            index_split[n_child - 1].push_back(i);
        }
        else {
            for (int k = 1; k < n_child - 1; ++k) {
                if (s_arr[i] >= split_vals[k - 1] && s_arr[i] < split_vals[k]) {
                    index_split[k].push_back(i);
                    break;
                }
            }
        }
    }
    return index_split;
}
