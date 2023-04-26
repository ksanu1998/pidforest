// #include "Forest.h"
// #include "Node.h"
// #include "PointSet.h"
// #include "Histogram.h"

class Node;
class PointSet;
class Histogram;
class Forest;

class Cube 
{
	public:
    Cube(Node* node, std::vector<float> start, std::vector<float> end) 
    {
        assert(node != nullptr);
        this->node = node;
        this->start = start;
        this->end = end;
        this->dim = start.size();
        this->split_axis = -1;
        this->vol = 0;
        for (int i = 0; i < this->dim; ++i) {
            this->vol += std::log(this->end[i] - this->start[i]);
        }
    }

    std::vector<int> filter_indices(std::vector<int> indices) 
    {
        std::vector<std::vector<bool>> in_lb(this->node->forest->points.rows(), std::vector<bool>(indices.size()));
        std::vector<std::vector<bool>> in_ub(this->node->forest->points.rows(), std::vector<bool>(indices.size()));
        for (int i = 0; i < this->dim; ++i) {
            float lb = this->start[i];
            float ub = this->end[i];
            for (int j = 0; j < indices.size(); ++j) {
                in_lb[indices[j]][j] = this->node->forest->points(i, indices[j]) >= lb;
                in_ub[indices[j]][j] = this->node->forest->points(i, indices[j]) < ub;
            }
        }
        std::vector<int> result;
        for (int i = 0; i < indices.size(); ++i) {
            bool all_in_lb = true;
            bool all_in_ub = true;
            for (int j = 0; j < this->node->forest->points.rows(); ++j) {
                all_in_lb &= in_lb[j][i];
                all_in_ub &= in_ub[j][i];
            }
            if (all_in_lb && all_in_ub) {
                result.push_back(indices[i]);
            }
        }
        return result;
    }

    std::vector<std::vector<int>> split_indices(np::ndarray pts, std::vector<int> indices) 
    {
	    int n_child = child.size();
	    if (n_child == 0) 
	    {
	        std::vector<std::vector<int>> ret;
	        ret.push_back(indices);
	        return ret;
	    }
	    int n_arr = indices.size();
	    if (n_arr == 0) 
	    {
	        std::vector<std::vector<int>> ret(n_child);
	        return ret;
	    }
	    np::ndarray s_arr = pts[split_axis];
	    int s_start = start[split_axis];
	    int s_end = end[split_axis];
	    std::vector<std::vector<int>> index_split(n_child);
	    index_split[0] = {};
	    index_split[n_child - 1] = {};
	    for (int ind : indices) 
	    {
	        if (s_arr[ind] >= s_start && s_arr[ind] < split_vals[0]) 
	        {
	            index_split[0].push_back(ind);
	        } 
	        else if (s_arr[ind] >= split_vals[n_child - 2] && s_arr[ind] < s_end) 
	        {
	            index_split[n_child - 1].push_back(ind);
	        } 
	        else 
	        {
	            for (int k = 1; k < n_child - 1; k++) 
	            {
	                if (s_arr[ind] >= split_vals[k - 1] && s_arr[ind] < split_vals[k]) 
	                {
	                    index_split[k].push_back(ind);
	                    break;
	                }
	            }
	        }
	    }
	    return index_split;
	}

};
