#include "Cube.h"
#include "Node.h"
#include "Forest.h"
#include "TimeSeries.h"

Cube::Cube(Node* n, int d, const std::vector<double>& s, const std::vector<double>& e)
    : node(n), dim(d), start(s), end(e) {
    dim = s.size();
    split_axis = -1;
    vol = 0;
    for (int i = 0; i < dim; ++i) {
        vol += std::log(end[i] - start[i]);
        if (vol > 100) {
            std::cout << ">>>>>>>>>>> Anomalous volume: " << vol << " " << end[i] << " " << start[i] << std::endl;
            // exit(1);
        }
        else {
            std::cout << ">>>>>>>>>>> Correct volume: " << vol << " " << end[i] << " " << start[i] << std::endl;
        }
    }
}

std::vector<int> Cube::filter_indices(const std::vector<int>& indices) {
    std::vector<int> filtered_indices;
    for (int i : indices) {
        bool in_lb = true;
        bool in_ub = true;
        // std::cout << ">>>>>>>>>>>>>>>>>>dim: " << dim << std::endl;
            for (int j = 0; j < dim; ++j) {
            // std::cout << ">>>>>>>>>>>>>>>>>>bool: " << (node->forest->points[j][i] < start[j] || node->forest->points[j][i] >= end[j]) << std::endl;
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
    std::vector<std::vector<double>> pts_t = pts; // TimeSeries::transpose(pts);
    //std::vector<int> s_arr = pts[split_axis];
    // std::vector<int> s_arr(pts[split_axis].begin(), pts[split_axis].end()); // wrong when split_axis = -1, this is not python!
    // Initialize s_arr based on split_axis
    // double s_start = start[split_axis]; // wrong when split_axis = -1, this is not python!
    // double s_end = end[split_axis]; // wrong when split_axis = -1, this is not python!
    std::vector<double> s_arr;
    double s_start, s_end;
    if (split_axis == -1) {
        // Access the last row when split_axis is -1
        s_arr = pts_t[pts_t.size() - 1]; //.back();
        s_start = start[start.size() - 1]; //.back(); // Access the last element of start
        s_end = end[end.size() - 1]; //.back();     // Access the last element of end
    } else {
        s_arr = pts_t[split_axis];
        s_start = start[split_axis];
        s_end = end[split_axis];
    }

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
                    // break;
                }
            }
        }
    }
    
    /*
    // Calculate indices for the first child
    index_split[0].reserve(indices.size());
    for (int ind : indices) {
        if (s_arr[ind] >= s_start && s_arr[ind] < split_vals[0]) {
            index_split[0].push_back(ind);
        }
    }

    // Calculate indices for the last child
    // index_split[n_child - 1].reserve(indices.size());
    index_split.back().reserve(indices.size());
    for (int ind : indices) {
        // if (s_arr[ind] >= split_vals[n_child - 2] && s_arr[ind] < s_end) {
        if (s_arr[ind] >= split_vals.back() && s_arr[ind] < s_end) {
            // index_split[n_child - 1].push_back(ind);
            index_split.back().push_back(ind);
        }
    }

    // Calculate indices for the remaining children
    for (int k = 1; k < n_child - 1; ++k) {
        index_split[k].reserve(indices.size());
        for (int ind : indices) {
            if (s_arr[ind] >= split_vals[k - 1] && s_arr[ind] < split_vals[k]) {
                index_split[k].push_back(ind);
            }
        }
    }
    */
    /*
    // Calculate indices for the first child
    for (int ind : indices) {
        if (s_arr[ind] >= s_start && s_arr[ind] < split_vals[0]) {
            index_split[0].push_back(ind);
        }
    }

    // Calculate indices for the last child
    for (int ind : indices) {
        if (s_arr[ind] >= split_vals[n_child - 2] && s_arr[ind] < s_end) {
            index_split[n_child - 1].push_back(ind);
        }
    }

    // Calculate indices for the remaining children
    for (int k = 1; k < n_child - 1; ++k) {
        for (int ind : indices) {
            if (s_arr[ind] >= split_vals[k - 1] && s_arr[ind] < split_vals[k]) {
                index_split[k].push_back(ind);
            }
        }
    }
    */
    /*
    // Calculate indices for the first child
    std::vector<int> temp_indices;
    for (int ind : indices) {
        if (s_arr[ind] >= s_start && s_arr[ind] < split_vals[0]) {
            temp_indices.push_back(ind);
        }
    }
    index_split[0].swap(temp_indices);

    // Calculate indices for the last child
    temp_indices.clear();
    for (int ind : indices) {
        if (s_arr[ind] >= split_vals[n_child - 2] && s_arr[ind] < s_end) {
            temp_indices.push_back(ind);
        }
    }
    index_split[n_child - 1].swap(temp_indices);

    // Calculate indices for the remaining children
    for (int k = 1; k < n_child - 1; ++k) {
        temp_indices.clear();
        for (int ind : indices) {
            if (s_arr[ind] >= split_vals[k - 1] && s_arr[ind] < split_vals[k]) {
                temp_indices.push_back(ind);
            }
        }
        index_split[k].swap(temp_indices);
    }
    */
    std::cout << "index_split (sizes) [CUBE]: [";
    for (int i = 0; i < index_split.size(); i++) {
        std::cout << index_split[i].size() << " ";
    }
    std::cout << "]" << std::endl;
    return index_split;
}
