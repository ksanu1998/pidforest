#include "Forest.h"
#include "Node.h"
#include <set>

Forest::Forest(const std::unordered_map<std::string, double>& kwargs) {
    n_trees = static_cast<int>(kwargs.at("n_trees"));
    max_depth = static_cast<int>(kwargs.at("max_depth"));
    max_samples = static_cast<int>(kwargs.at("max_samples"));
    max_buckets = static_cast<int>(kwargs.at("max_buckets"));
    epsilon = kwargs.at("epsilon");
    sample_axis = kwargs.at("sample_axis");
    threshold = kwargs.at("threshold");
    tree.reserve(n_trees);
    n_leaves.resize(n_trees);
}

void Forest::fit(const std::vector<std::vector<double>>& pts) {
    points = pts;
    size = static_cast<int>(points.size());
    dim = static_cast<int>(points[0].size());
    std::cout << "Dim " << dim << ", Size " << size << std::endl;
    if (sample_axis * dim < 1) {
        std::cout << "sample_axis is too low" << std::endl;
        return;
    }
    start.resize(dim);
    end.resize(dim);
    
    std::vector<std::vector<double>> transposed_points(dim, std::vector<double>(size));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < dim; ++j) {
            transposed_points[j][i] = points[i][j];
        }
    }

    for (int axis = 0; axis < dim; ++axis) {
        std::vector<double> val(transposed_points[axis].begin(), transposed_points[axis].end());
        std::sort(val.begin(), val.end());
        val.erase(std::unique(val.begin(), val.end()), val.end());
        start[axis] = (3 * val[0] - val[1]) / 2;
        end[axis] = (3 * val.back() - val[val.size() - 2]) / 2;
    }
    start.resize(dim);
    end.resize(dim);
    std::vector<int> sample(std::min(size, max_depth * 200));
    std::iota(sample.begin(), sample.end(), 0);
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(sample.begin(), sample.end(), g);
    std::cout << "\n >> Building trees\n";
    for (int i = 0; i < n_trees; ++i) {
        
        std::vector<int> indices(max_samples);
        std::iota(indices.begin(), indices.end(), 0);
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(indices.begin(), indices.end(), g);
        std::cout << "\n >> Building tree #" << i + 1<< "/" << n_trees <<"\n";
        Node root_node(0, this, {{"indices", indices}});
        root_node.compute_density(sample);
        tree.push_back(root_node);
        n_leaves[i] = root_node.compute_leaf_num();
        std::cout << "\n >> Building tree #" << i + 1 << "/" << n_trees <<" [DONE]\n";
    }
    std::cout << "\n >> Building trees [DONE]\n";    
}

std::tuple<std::vector<int>, std::unordered_map<int, std::vector<double>>,
        std::unordered_map<int, std::vector<double>>, std::unordered_map<int, double>, std::vector<double>>
Forest::predict(const std::vector<std::vector<double>>& pts, double err, double pct) {
    
    int n_pts = static_cast<int>(pts.size());
    std::vector<std::vector<double>> scores(n_trees, std::vector<double>(n_pts));
    std::vector<int> indices;
    indices.reserve(n_pts);
    for (int i = 0; i < n_pts; i++) {
        indices.push_back(i);
    }

    for (int i = 0; i < n_trees; ++i) {
        tree[i].compute_split(pts, indices, scores[i]);
    }
    int n_err = static_cast<int>(err * n_pts);
    
    std::vector<double> min_score(n_pts);
    for (int i = 0; i < n_pts; ++i) {
        std::vector<double> column_scores(n_trees);
        bool all_zeroes = 1;
        for (int j = 0; j < n_trees; ++j) {
            column_scores[j] = scores[j][i];
            if (column_scores[j] > 0 || column_scores[j] < 0) {
                all_zeroes = 0;
            }
        }
        std::sort(column_scores.begin(), column_scores.end());
        int index = static_cast<int>(pct / 100.0 * (n_trees - 1)); // since pct = 0, this always sets index to 0 (min)
        min_score[i] = column_scores[index];
    }
    
    std::vector<int> top_indices(n_pts);
    std::iota(top_indices.begin(), top_indices.end(), 0);

    std::sort(top_indices.begin(), top_indices.end(), [&](int a, int b) {
        return min_score[a] < min_score[b];
    });

    top_indices.resize(n_err);

    std::unordered_map<int, std::vector<double>> anom_pts;
    std::unordered_map<int, std::vector<double>> anom_scores;
    std::unordered_map<int, double> anom_pct;

    for (int i = 0; i < n_err; ++i) {
        int index = top_indices[i];
        anom_pts[index] = std::vector<double>(n_trees);
        anom_scores[index] = std::vector<double>(n_trees);
        for (int j = 0; j < n_trees; ++j) {
            anom_pts[index][j] = pts[j][index];
            anom_scores[index][j] = scores[j][index];
        }
        anom_pct[index] = min_score[index];
    }
    std::cout << "min_score" << std::endl;
    for (int i = 0; i < min_score.size(); i++) {
        if (min_score[i] < -500) {
            std::cout << "xxx " << min_score[i] << std::endl;
        }
        else {
                std::cout << min_score[i] << std::endl;
        }
    }
    std::cout << std::endl << std::endl;
    return std::make_tuple(top_indices, anom_pts, anom_scores, anom_pct, min_score);
}
