#ifndef POINT_SET_H
#define POINT_SET_H

#include <unordered_set>
#include <vector>

class Node;  // Forward declaration of the Node class

class PointSet {
public:
    Node* node;
    std::unordered_set<int> indices;
    std::vector<std::vector<double>> val;
    std::vector<std::vector<double>> count;
    std::vector<std::vector<double>> gap;

    PointSet(Node* n, const std::unordered_set<int>& inds);


    int size() const;
};

#endif  // POINT_SET_H
