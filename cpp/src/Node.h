#ifndef C_ISOLATION_NODE_H
#define C_ISOLATION_NODE_H

#include "Tree_Param.h"
#include <cstdlib>
#include <random>
#include "Split.h"
#include "Histogram.h"
#include "HistogramMemory.h"

class Node {

    public:
        Node(double **arr, int dim, int size, int m_depth, Tree_Param * tParam);
        void findSplit();
        void getBuckets(int ** buckets, double * bucketErrors);
        void createChildren();

    public:
        virtual ~Node();

    public:
        int randomChoice(double  * arr, int size);

    private:
        Histogram histogram_;
        int dim_;
        double ** arr_;
        double density_;
        int size_;
        Tree_Param * t_param_;
        std::mt19937 rng_;
    
    private:
        int m_depth_;
        Split split_;
        int split_dim;
        double * split_buckets;
        Node * children;
        HistogramMemory * hist_mem_;
};


#endif