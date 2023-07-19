// /*************************
//  * This class holds the global parameters of the forest
//  * **************************/
// #ifndef ALLC_FORESTPARAM_H
// #define ALLC_FORESTPARAM_H
// #include <vector>

// using namespace std;

// class ForestParam {
// public:
//     ForestParam(const int ntrees, const double epsilon, const int depth, const int nbucket, const int nsamples);
//     const int n_trees_; //number of trees in the forest
//     const double epsilon_; // approximation factor for the histograms
//     const int depth_; // maximum depth of each tree
//     const int n_bucket_; // the number of partitions in each histogram
//     const int n_samples_; // number of samples taken for each tree
// };


// #endif //ALLC_FORESTPARAM_H

#ifndef FORESTPARAM_H
#define FORESTPARAM_H

#include <unordered_map>

class ForestParam {
public:
    ForestParam(int ntrees, double epsilon, int depth, int nbucket, int nsamples);

    int n_trees_;
    double epsilon_;
    int depth_;
    int n_bucket_;
    int n_samples_;
};

#endif // FORESTPARAM_H
