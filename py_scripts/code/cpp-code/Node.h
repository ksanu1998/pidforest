// #include "Forest.h"
// #include "PointSet.h"
// #include "Cube.h"
// #include "Histogram.h"

class PointSet;
class Cube;
class Histogram;
class Forest;

class Node {
public:
    int depth;
    Forest* forest;
    std::vector<int> id_string;
    Cube* cube;
    PointSet* point_set;
    int density;
    std::vector<Node*> child;

    Node(int depth, Forest* forest, std::map<std::string, int> kwargs) {
        this->depth = depth;
        this->forest = forest;
        if (this->depth == 0) {
            this->id_string = {0};
            this->cube = new Cube(this, this->forest->start, this->forest->end);
            this->point_set = new PointSet(this, kwargs["indices"]);
        } else {
            this->id_string = kwargs["id"];
            this->cube = new Cube(this, kwargs["start"], kwargs["end"]);
            this->point_set = new PointSet(this, this->cube->filter_indices(kwargs["indices"]));
        }
        this->density = -1;
        if (this->depth < this->forest->max_depth && this->point_set->indices.size() > 1) {
            this->find_split();
        }
    }

    void find_split() 
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
	    int max_axes = std::min(static_cast<int>(imp_axis.size()), static_cast<int>(forest.sample_axis * cube.dim));
	    std::vector<int> s_axes = random_choice(imp_axis, max_axes, false);
	    std::unordered_map<int, std::vector<double>> buckets;
	    std::unordered_map<int, double> var_red;
	    for (int axis : s_axes) {
	        auto hist = Histogram(point_set.gap[axis] / point_set.count[axis], point_set.count[axis], forest.max_buckets, forest.epsilon);
	        auto [_, var, bkt] = hist.best_split();
	        var_red[axis] = var;
	        buckets[axis] = bkt;
	    }
	    if (*std::max_element(var_red.begin(), var_red.end(), [](const auto& p1, const auto& p2) {
	            return p1.second < p2.second;
	        }) <= forest.threshold) {
	        return;
	    }
	    int split_axis = random_choice(s_axes, get_probs(var_red));
	    cube.split_axis = split_axis;
	    std::vector<double> split_vals;
	    for (int i = 0; i < buckets[split_axis].size(); ++i) {
	        int idx = static_cast<int>(buckets[split_axis][i]);
	        double val = (point_set.val[split_axis][idx - 1] + point_set.val[split_axis][idx]) / 2;
	        split_vals.push_back(val);
	    }
	    cube.split_vals = std::move(split_vals);
	    for (int i = 0; i < static_cast<int>(cube.split_vals.size()) + 1; ++i) {
	        std::vector<double> new_start(cube.start);
	        std::vector<double> new_end(cube.end);
	        if (i > 0 && i < static_cast<int>(cube.split_vals.size())) {
	            new_start[split_axis] = cube.split_vals[i - 1];
	            new_end[split_axis] = cube.split_vals[i];
	        } else if (i == 0) {
	            new_end[split_axis] = cube.split_vals[0];
	        } else {  // i == cube.split_vals.size()
	            new_start[split_axis] = cube.split_vals.back();
	        }
	        std::vector<int> new_id(id_string);
	        new_id.push_back(i);
	        Node child_node(depth + 1, forest, new_start, new_end, point_set.indices, new_id);
	        child.emplace_back(child_node);
	        cube.child.emplace_back(child_node.cube);
	    }
	}

	void compute_density(std::vector<int> indices) 
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
	        auto index_split = cube.split_indices(forest->points, indices);
	        for (size_t i = 0; i < child.size(); i++) {
	            child[i].compute_density(index_split[i]);
	        }
	    }
	}

	int compute_leaf_num() 
	{
	    if (!child.empty()) {
	        int leaves = 0;
	        for (int i = 0; i < child.size(); i++) {
	            leaves += child[i].compute_leaf_num();
	        }
	        return leaves;
	    } else {
	        return 1;
	    }
	}

	void compute_split(const std::vector<std::vector<float>>& pts, const std::vector<int>& indices, std::vector<float>& scores) 
	{
	    if (!child.empty()) {
	        std::vector<std::vector<int>> index_split = cube.split_indices(pts, indices);
	        for (int i = 0; i < child.size(); i++) {
	            if (index_split[i].size() > 0) {
	                child[i].compute_split(pts, index_split[i], scores);
	            }
	        }
	    } else {
	        for (int i = 0; i < indices.size(); i++) {
	            scores[indices[i]] = density;
	        }
	    }
	}

	Cube compute_box(const PointSet &pts) const 
	{
	    if (child.size()) {
	        const int n_child = child.size();
	        const double *s_arr = pts.val[cube.split_axis].data();
	        const double s_start = cube.start[cube.split_axis];
	        const double s_end = cube.end[cube.split_axis];
	        if ((s_arr >= s_start) && (s_arr < cube.split_vals[0])) {
	            return child[0].compute_box(pts);
	        }
	        else if ((s_arr >= cube.split_vals[n_child - 2]) && (s_arr < s_end)) {
	            return child[n_child - 1].compute_box(pts);
	        }
	        else {
	            for (int k = 1; k < n_child - 1; k++) {
	                if ((s_arr >= cube.split_vals[k - 1]) && (s_arr < cube.split_vals[k])) {
	                    return child[k].compute_box(pts);
	                }
	            }
	        }
	    }
    	return cube;
	}

	void compute_depth(const PointSet &pts, const std::vector<int> &indices, std::vector<double> &scores) const 
	{
	    if (child.size()) {
	        std::vector<std::vector<int>> index_split = cube.split_indices(pts, indices);
	        for (int i = 0; i < child.size(); i++) {
	            if (index_split[i].size() > 0) {
	                child[i].compute_depth(pts, index_split[i], scores);
	            }
	        }
	    }
	    else {
	        for (int i = 0; i < indices.size(); i++) {
	            scores[indices[i]] = depth;
	        }
	    }
	}


	std::string __str__() 
	{
	    std::stringstream ss;
	    ss << "Id: ";
	    for (auto id : this->id_string) {
	        ss << id << ", ";
	    }
	    ss << std::endl;
	    ss << "Boundary: ";
	    for (int i = 0; i < this->cube.dim; ++i) {
	        ss << " [" << this->cube.start[i] << ", " << this->cube.end[i] << "]";
	        if (i < this->cube.dim - 1) {
	            ss << " x";
	        } else {
	            ss << std::endl;
	        }
	    }
	    ss << "Points:\n " << this->point_set.points.transpose() << std::endl;
	    ss << "Indices: " << this->point_set.indices << std::endl;
	    return ss.str();
	}

	void print_node() 
	{
	    std::vector<Node*> print_list = {this};
	    while (!print_list.empty()) {
	        Node* node = print_list[0];
	        print_list.erase(print_list.begin());
	        std::cout << node->__str__();
	        for (auto child : node->child) {
	            print_list.push_back(child);
	        }
	    }
	}


};
