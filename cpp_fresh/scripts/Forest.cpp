#include "Forest.h"
#include "Node.h" // Include the Node class header file if available

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
    dim = static_cast<int>(points.size());
    size = static_cast<int>(points[0].size());
    std::cout << "FOREST FIT REACHED" << std::endl;
    if (sample_axis * dim < 1) {
        std::cout << "sample_axis is too low" << std::endl;
        return;
    }

    start.resize(dim);
    end.resize(dim);

    for (int axis = 0; axis < dim; ++axis) {
        std::vector<double> val(points[axis].begin(), points[axis].end());
        std::sort(val.begin(), val.end());
        val.erase(std::unique(val.begin(), val.end()), val.end());

        if (val.size() <= 1) {
            std::cout << "No entropy in dimension: " << axis << std::endl;
            return;
        }

        start[axis] = (3 * val[0] - val[1]) / 2;
        end[axis] = (3 * val.back() - val[val.size() - 2]) / 2;
    }

    std::vector<int> sample(max_depth * 200);
    std::iota(sample.begin(), sample.end(), 0);
    std::shuffle(sample.begin(), sample.end(), std::default_random_engine(17));
    std::cout << "STARTING TO BUILD TREE" << std::endl;
    
    for (int i = 0; i < n_trees; ++i) {
        std::vector<int> indices(max_samples);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), std::default_random_engine());
        std::cout << "BUILDING TREE: " << i << std::endl;
        Node root_node(0, this, {{"indices", indices}});
        std::cout << "FINISHED BUILDING ROOT NODE OF TREE: " << i << std::endl;
        // Node root_node(0, this);
        root_node.compute_density(sample);
        tree.push_back(root_node);
        n_leaves[i] = root_node.compute_leaf_num();
    }
    std::cout << "FINISHED BUILDING TREE" << std::endl;
    
}

std::tuple<std::vector<int>, std::unordered_map<int, std::vector<double>>,
        std::unordered_map<int, std::vector<double>>, std::unordered_map<int, double>, std::vector<double>>
Forest::predict(const std::vector<std::vector<double>>& pts, double err, double pct) {
    int n_pts = static_cast<int>(pts[0].size());
    std::vector<std::vector<double>> scores(n_trees, std::vector<double>(n_pts));
    std::vector<int> indices(n_pts);
    std::iota(indices.begin(), indices.end(), 0);
    std::cout << "COMPUTING SPLIT" << std::endl;
    for (int i = 0; i < n_trees; ++i) {
        tree[i].compute_split(pts, indices, scores[i]);
    }
    std::cout << "FINISHED COMPUTING SPLIT" << std::endl;
    
    int n_err = static_cast<int>(err * n_pts);
    std::vector<double> min_score(n_pts);
    std::vector<int> top_indices(n_err);
    std::iota(top_indices.begin(), top_indices.end(), 0);
    std::cout << "BEGIN partial_sort" << std::endl;
    // std::partial_sort(top_indices.begin(), top_indices.begin() + n_err, top_indices.end(),
    // [&scores, n_pts](int a, int b) {
    //     return scores[a][b] < scores[a][n_pts];
    // });
    std::partial_sort(top_indices.begin(), top_indices.begin() + n_err, top_indices.end(),
    [&scores, n_pts](int a, int b) {
        return scores[a][n_pts] < scores[b][n_pts];
    });

    std::cout << "FINISHED partial_sort" << std::endl;
    
    std::unordered_map<int, std::vector<double>> anom_pts;
    std::unordered_map<int, std::vector<double>> anom_scores;
    std::unordered_map<int, double> anom_pct;

    for (int i = 0; i < n_err; ++i) {
        int index = top_indices[i];
        anom_pts[index] = pts[index];
        anom_scores[index] = scores[index];
        anom_pct[index] = min_score[index];
    }
    
    return std::make_tuple(top_indices, anom_pts, anom_scores, anom_pct, min_score);
}
