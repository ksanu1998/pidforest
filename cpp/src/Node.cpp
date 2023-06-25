#include <random>
#include <math.h>
#include <iostream>
#include "Node.h"
#include "util_functions.h"


Node::Node(double **arr,  int dim,  int size,  int m_depth, Tree_Param * tParam) : arr_(arr), size_(size), dim_(dim), m_depth_(m_depth), t_param_(tParam) {
    // Node should probably also get the bounding box of the array.
    histogram_.setEpsilon(t_param_->epsilon_);
    histogram_.setNBuckets(t_param_->n_bucket_);
    int maxLength;
    if (t_param_->epsilon_ == 0)
        maxLength  = size; 
    else {
        maxLength = log(size)/log(1 + t_param_->epsilon_);
        if (maxLength > size) maxLength = size;
    }
    hist_mem_ = new HistogramMemory(maxLength, size_, tParam->n_bucket_);
    std::random_device rd;
    rng_.seed(rd()); 
}

//Samples an index i of an array with probability proportional to arr[i]
int Node::randomChoice(double *arr, int size) {
    if (size <= 0) {
        std::cerr<<"Negative size in randomChoice";
        return -1;
    }
    double sum = 0;
    for (int i = 0; i < size; i++) {
        if (arr[i] < 0) {
            std::cerr<<"Negative value in randomChoice";
            return -1; //todo: throw exception
        }
        sum += arr[i];
    }
    if (utils::isEqual(sum, 0)) {
        std::cerr<<"Sum zero in randomChoice";
        return -1;
    } 

    std::uniform_real_distribution<double> dist; // {0.0, 1.0};
    double x = dist(rng_);
    double curr = 0;
    for (int i = 0; i < size; i++) {
        curr += (arr[i]/sum);
        if (curr >= x)
            return i;
    }
    return size - 1;
}

void Node::findSplit() {
    Split splits[dim_];
    double err_red[dim_];
    for (int dim = 0 ;dim < dim_; dim ++) {
        histogram_.ComputeAppBuckets(arr_[dim], size_, hist_mem_, &splits[dim]);
        err_red[dim] = splits[dim].bucket_error_;
    }
    //sample to get the axis to split upon.
    int index = randomChoice(err_red, dim_);
    split_.bucket_error_ = splits[index].bucket_error_;
    split_.bucket_level_ = splits[index].bucket_level_;
    for (int k = 0; k < split_.bucket_level_; k++)
        split_.bucket_list_[k] = splits[index].bucket_list_[k];
    split_.dim_ = index;

    //todo: split the nodes and continue
}

Node::~Node() {
    delete hist_mem_;
}
