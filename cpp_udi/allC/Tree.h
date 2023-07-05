/**
 * The class Tree gets a pointer to the data set and the parameters from the forest and computes a PID Tree.
 * It has an inner class Node. The tree has a member which is the root Node. This node would span
 * the rest of the Nodes until there is no data or the depth was reached.
 * Tree is intended to be single thread.
 * */

#ifndef ALLC_TREE_H
#define ALLC_TREE_H

#include <vector>
#include "ForestParam.h"
#include "DataSet.h"
#include <random>

using namespace std;

class Tree {
    class Node {
        public:
            //Node(Tree & mTree, int depth);
            Node(int depth);
            ~Node();
            void computeTree(vector <int> & indices, const vector <double> & cubeLower, const vector <double> & cubeUpper);
            void predict(const DataSet * data, vector<double> & scores, const vector<int> & indices);
            void setTree( Tree * parent);
            void printCube(); //for debugging purposes
            void printSubTree(); //prints the subtree rooted at this node.

        private:
            //const Tree & mTree_; // a pointer back to the tree.
            Tree * mTree_;
            vector<double> cubeLB_; // holds the beginning of the cube the node holds, it is a partition of space
            vector<double> cubeUB_; // holds the end points of the cube. Will have length [dim]
            vector<int> fitIndices_; // holds the indices of the data that was fit into that node.
            double volume_;
            double density_;
            int depth_; // the depth of the node. The root is depth 0.
            //vector<Node> children_;
            vector<Node *> m_children_;

            void filterIndices(const DataSet * data, vector<int> & newIndices, const vector<int> & oldIndices);
            void filterIndices(vector<int> & oldIndices);
            void computeVolume();
            void createGaps(const vector<int> & indices, vector< vector<double> > & gaps, vector< vector<int> > & counts); // computes the vectors of gaps
    };

public:
    Tree(int id, ForestParam & ForestParam);
    //~Tree();
    void computeTree(const DataSet * data);
    mt19937 mPrg_;  // This is the PRG that supplies randomnness to all the nodes and the histograms in the tree
    const DataSet * dataSet_;
    void printTree(); // prints the tree

    void predict(const DataSet * data, vector<double> & scores);

private:
    Node root_;
    ForestParam forestParam_;
    int id_;
    void getSample(vector <int> & indices);
};

#endif //ALLC_TREE_H