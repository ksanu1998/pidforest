#include "PointSet.h"
#include "Node.h"
#include "Forest.h"
#include <algorithm>

PointSet::PointSet(Node* n, const std::unordered_set<int>& inds) : node(n), indices(inds) {
    int dim = node->forest->dim;
    val.resize(dim);
    count.resize(dim);
    gap.resize(dim);
    // std::cout << "node->forest->dim, " << dim << std::endl;
    for (int axis = 0; axis < dim; ++axis) {
        std::vector<double> values;
        for (int index : indices) {
            double value = node->forest->points[index][axis];
            values.push_back(value);
        }
        std::sort(values.begin(), values.end());
        auto unique_end = std::unique(values.begin(), values.end());
        values.erase(unique_end, values.end());
        
        val[axis] = values;
        count[axis].resize(values.size(), 1);
        for (int index : indices) {
            double value = node->forest->points[index][axis];
            auto it = std::lower_bound(values.begin(), values.end(), value);
            int idx = std::distance(values.begin(), it);
            ++count[axis][idx];
        }
        
        /*
        std::cout << "count:" << std::endl;
        std::cout << "[" << std::endl;
        for (const auto& axis_count : count) {
            std::cout << "    [";
            for (size_t i = 0; i < axis_count.size(); ++i) {
                std::cout << axis_count[i];
                if (i < axis_count.size() - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "]" << std::endl;
        }
        std::cout << "]" << std::endl;
        */
        if (values.size() <= 1) {
            gap[axis] = {0.0};
        } else {
            // std::cout << ">>>>>>>>>>>>>>>>PointSet[BET1]" << std::endl;
            // std::cout << gap[axis].size() << ", " << values.size() << std::endl;
            gap[axis].resize(values.size());
            // std::cout << ">>>>>>>>>>>>>>>>PointSet[BET2]" << std::endl;
            gap[axis][0] = (values[1] - values[0]) / 2;
            // std::cout << ">>>>>>>>>>>>>>>>PointSet[BET3]" << std::endl;
            gap[axis][values.size() - 1] = (values[values.size() - 1] - values[values.size() - 2]) / 2;
            // std::cout << ">>>>>>>>>>>>>>>>PointSet[BET4]" << std::endl;
            for (int i = 1; i < values.size() - 1; ++i) {
                gap[axis][i] = (values[i + 1] - values[i - 1]) / 2;
            }
            // std::cout << ">>>>>>>>>>>>>>>>PointSet[BET][DONE]" << std::endl;
        }
    }
    // std::cout << ">>>>>>>>>>>>>>>>PointSet[DONE]" << std::endl;
}

int PointSet::size() const {
    return indices.size();
}