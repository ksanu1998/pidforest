// #include "ForestParam.h"

// ForestParam::ForestParam(const int ntrees, const double epsilon, const int depth, const int nBucket, const int nSamples) : 
//                                                                                                      n_trees_(ntrees),
//                                                                                                      epsilon_(epsilon), 
//                                                                                                      depth_(depth),
//                                                                                                      n_bucket_(nBucket),
//                                                                                                      n_samples_(nSamples) { }

#include "ForestParam.h"

ForestParam::ForestParam(int ntrees, double epsilon, int depth, int nbucket, int nsamples)
    : n_trees_(ntrees), epsilon_(epsilon), depth_(depth), n_bucket_(nbucket), n_samples_(nsamples) {}
