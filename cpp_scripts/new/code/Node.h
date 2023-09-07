// Node.h

#ifndef NODE_H
#define NODE_H

#include <vector>
#include <unordered_map>
#include <variant>
#include <array>

#include "Forest.h"
#include "Cube.h"
#include "PointSet.h"
#include "Histogram.h"

class Node {
public:
    Node(int depth, Forest* forest, const std::unordered_map<std::string, std::variant<std::vector<int>, std::vector<double>>>& kwargs);
    Forest* forest;
    void compute_split(const std::vector<std::vector<double>>& pts, const std::vector<int>& indices, std::vector<double>& scores);
    void compute_density(const std::vector<int>& indices);
    int compute_leaf_num();
    
private:
    int depth;
    int dim;
    std::vector<int> id_string;
    Cube cube;
    PointSet point_set;
    double density;
    std::vector<Node> child;

    void find_split();
    
};

#endif // NODE_H