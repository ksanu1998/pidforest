// #include "Forest.h"
// #include "Node.h"
// #include "Cube.h"
// #include "Histogram.h"

class Node;
class Cube;
class Histogram;
class Forest;

class PointSet {
public:
    Node* node;
    std::vector<int> indices;
    std::vector<std::vector<double>> val;
    std::vector<std::vector<int>> count;
    std::vector<std::vector<double>> gap;
    
    PointSet(Node* node, std::vector<int> indices) {
        this->node = node;
        this->indices = indices;
        this->val = std::vector<std::vector<double>>(this->node->forest.dim);
        this->count = std::vector<std::vector<int>>(this->node->forest.dim);
        this->gap = std::vector<std::vector<double>>(this->node->forest.dim);
        for (int axis = 0; axis < this->node->forest.dim; ++axis) {
            std::vector<double> val_vec, gap_vec;
            std::vector<int> count_vec;
            std::unordered_map<double, int> val_to_count;
            for (int i = 0; i < indices.size(); ++i) {
                double pt = this->node->forest.points[axis][indices[i]];
                if (val_to_count.find(pt) == val_to_count.end()) {
                    val_to_count[pt] = 1;
                } else {
                    val_to_count[pt]++;
                }
            }
            for (auto it = val_to_count.begin(); it != val_to_count.end(); ++it) {
                val_vec.push_back(it->first);
                count_vec.push_back(it->second);
            }
            this->val[axis] = val_vec;
            this->count[axis] = count_vec;
            int n_val = val_vec.size();
            gap_vec.resize(n_val);
            if (n_val <= 1) {
                gap_vec[0] = 0;
            } else {
                gap_vec[0] = (val_vec[0] + val_vec[1]) / 2 - this->node->cube.start[axis];
                gap_vec[n_val - 1] = this->node->cube.end[axis] - (val_vec[n_val - 1] + val_vec[n_val - 2]) / 2;
                for (int i = 1; i < n_val - 1; ++i) {
                    gap_vec[i] = (val_vec[i + 1] - val_vec[i - 1]) / 2;
                }
            }
            this->gap[axis] = gap_vec;
        }
    }
};