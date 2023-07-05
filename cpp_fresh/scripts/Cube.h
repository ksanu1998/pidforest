#ifndef CUBE_H
#define CUBE_H

#include <vector>

class Node;  // Forward declaration of the Node class

class Cube {
public:
    Node* node;
    std::vector<Cube*> child;
    std::vector<double> start;
    std::vector<double> end;
    int dim;
    int split_axis;
    std::vector<double> split_vals;
    double vol;

    Cube(Node* n, const std::vector<double>& s, const std::vector<double>& e);

    std::vector<int> filter_indices(const std::vector<int>& indices);
    std::vector<std::vector<int>> split_indices(const std::vector<std::vector<double>>& pts, const std::vector<int>& indices);
};

#endif  // CUBE_H
