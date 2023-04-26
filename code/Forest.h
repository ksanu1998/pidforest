// #include "Node.h"
// #include "PointSet.h"
// #include "Cube.h"
// #include "Histogram.h"

class Node;
class PointSet;
class Cube;
class Histogram;

class Forest
{
	public:
	Forest(std::map<std::string, int> kwargs)
	{
		n_trees = kwargs["n_trees"];
        max_depth = kwargs["max_depth"];
        max_samples = kwargs["max_samples"];
        max_buckets = kwargs["max_buckets"];
        epsilon = kwargs["epsilon"];
        sample_axis = kwargs["sample_axis"];
        threshold = kwargs["threshold"];
        tree = std::vector<int>{};
        n_leaves = std::vector<int>(n_trees, 0);
	}

	private:
    int n_trees;
    int max_depth;
    int max_samples;
    int max_buckets;
    double epsilon;
    int sample_axis;
    double threshold;
    std::vector<int> tree;
    std::vector<int> n_leaves;

    void fit(std::vector<std::vector<double>> pts) 
    {
	    std::vector<std::vector<double> > points = pts;
	    int dim = points.size();
	    int size = points[0].size();
	    if (sample_axis * dim == 0) 
	    {
	        std::cout << "sample_axis is too low" << std::endl;
	        return;
    	}
	    std::vector<double> start = std::vector<double>(dim, 0.0);
	    std::vector<double> end = std::vector<double>(dim, 0.0);
	    for (int axis = 0; axis < dim; axis++) 
	    {
	        std::vector<double> val = points[axis];
	        std::sort(val.begin(), val.end());
	        auto last = std::unique(val.begin(), val.end());
	        val.erase(last, val.end());
	        if (val.size() <= 1) 
	        {
	            std::cout << "No entropy in dimension :" << axis << std::endl;
	            return;
	        }
        start[axis] = (3 * val[0] - val[1]) / 2.0;
        end[axis] = (3 * val.back() - val[val.size()-2]) / 2.0;
    	}
	    Node k_args{0, this};
	    int max_sample_size = std::min(size, max_depth*200);
	    std::vector<int> sample(max_sample_size);
	    std::iota(sample.begin(), sample.end(), 0);
	    std::shuffle(sample.begin(), sample.end(), std::default_random_engine{});
	    for (int i = 0; i < n_trees; i++) 
	    {
	        std::vector<int> indices(max_samples);
	        std::iota(indices.begin(), indices.end(), 0);
	        std::shuffle(indices.begin(), indices.end(), std::default_random_engine{});
	        std::transform(indices.begin(), indices.end(), indices.begin(),
	                       [sample](int i) { return sample[i]; });
	        k_args.indices = indices;
	        auto root_node = new Node{k_args};
	        root_node->compute_density(sample);
	        tree.push_back(root_node);
	        n_leaves[i] = root_node->compute_leaf_num();
	    }
	}

	std::vector<int> predict(std::vector<std::vector<double>>& pts, double err=0.1, double pct=50) 
	{
	    int n_trees = tree.size();
	    int n_pts = pts[0].size();
	    std::vector<std::vector<double>> scores(n_trees, std::vector<double>(n_pts, 0));
	    std::vector<int> indices(n_pts);
	    iota(indices.begin(), indices.end(), 0);
	    for (int i = 0; i < n_trees; i++) 
	    {
	        tree[i].compute_split(pts, indices, scores[i]);
	    }
	    int n_err = floor(err * n_pts);
	    std::vector<double> min_score(n_pts);
	    for (int i = 0; i < n_pts; i++) 
	    {
	        std::vector<double> col(n_trees);
	        for (int j = 0; j < n_trees; j++) {
	            col[j] = scores[j][i];
	        }
	        sort(col.begin(), col.end());
	        min_score[i] = col[floor(pct / 100.0 * n_trees)];
	    }
	    std::vector<int> top_indices(n_err);
	    std::vector<double> min_scores_copy = min_score;
	    sort(min_scores_copy.begin(), min_scores_copy.end());
	    for (int i = 0; i < n_err; i++) 
	    {
	        auto it = find(min_score.begin(), min_score.end(), min_scores_copy[i]);
	        int index = distance(min_score.begin(), it);
	        top_indices[i] = index;
	        min_score[index] = NAN;
	    }
	    std::unordered_map<int, std::vector<double>> anom_pts;
	    std::unordered_map<int, std::vector<double>> anom_scores;
	    std::unordered_map<int, double> anom_pct;
	    for (int i = 0; i < n_err; i++) 
	    {
	        int index = top_indices[i];
	        std::vector<double> pt(pts.size());
	        for (int j = 0; j < pts.size(); j++) 
	        {
	            pt[j] = pts[j][index];
	        }
	        anom_pts[index] = pt;
	        anom_scores[index] = scores[index];
	        anom_pct[index] = min_score[index];
	    }
    	return top_indices;
	}

	std::tuple<std::vector<int>, std::unordered_map<int, std::vector<double>>,
           std::unordered_map<int, std::vector<double>>, std::unordered_map<int, double>, std::vector<double>>
	predict_depth(const Eigen::MatrixXd& pts, double err, double pct) 
	{
	    int n_pts = pts.cols();
	    Eigen::MatrixXd scores(n_trees, n_pts);
	    std::vector<int> indices(n_pts);
	    std::iota(indices.begin(), indices.end(), 0);

	    for (int i = 0; i < n_trees; ++i) {
	        tree[i]->compute_depth(pts, indices, scores.row(i));
	    }

	    int n_err = static_cast<int>(err * n_pts);
	    int k = std::floor(pct / 100.0 * n_trees);
	    Eigen::VectorXd min_score = scores.colwise().percentile(k);

	    std::vector<int> top_indices(n_err);
	    std::iota(top_indices.begin(), top_indices.end(), 0);
	    std::partial_sort(top_indices.begin(), top_indices.begin() + n_err, top_indices.end(),
	                      [&](int i1, int i2) { return min_score[i1] < min_score[i2]; });

	    std::unordered_map<int, std::vector<double>> anom_pts;
	    std::unordered_map<int, std::vector<double>> anom_scores;
	    std::unordered_map<int, double> anom_pct;
	    for (int i = 0; i < n_err; ++i) {
	        int idx = top_indices[i];
	        std::vector<double> pts_col(pts.rows());
	        std::vector<double> scores_col(n_trees);
	        for (int j = 0; j < pts.rows(); ++j) {
	            pts_col[j] = pts(j, idx);
	        }
	        for (int j = 0; j < n_trees; ++j) {
	            scores_col[j] = scores(j, idx);
	        }
	        anom_pts.emplace(idx, pts_col);
	        anom_scores.emplace(idx, scores_col);
	        anom_pct.emplace(idx, min_score(idx));
	    }

	    return std::make_tuple(top_indices, anom_pts, anom_scores, anom_pct, std::vector<double>(min_score.data(),
                                                                                             min_score.data() + n_pts));
	}

};