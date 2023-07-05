#ifndef C_ISOLATION_TREE_PARAM_H
#define C_ISOLATION_TREE_PARAM_H


class Tree_Param {
public:
    Tree_Param(const double epsilon, const int depth, const int nBucket, const int nSamples);

    const double epsilon_;
    const int depth_;
    const int n_bucket_;
    const int n_samples;
};


#endif
