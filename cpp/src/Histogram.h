#ifndef C_ISOLATION_HISTOGRAM_H
#define C_ISOLATION_HISTOGRAM_H

#include "Split.h"
#include "HistogramMemory.h"


class Histogram {
public:
    void ComputeBuckets(const double * arr, const int length, double * sum, double * sqsum); //computes the optimal histogram using a dynamic program
    int **getBucketsList() const;
    double getBucketsErr(int index) const;
    Histogram(const int nBuckets, const double epsilon);
    //todo: put all parameters in a struct or class. Perhaps including the scratch space and out variables.
    /*
     * eps: approximation for bucket computation
     * sum: double[length] will hold partial sums
     * sqsum: double[length] will hold partial square sums
     * b_val: int[n_buckets][max_list_length] will hold internal b-values for histogram algorithm.
     * app_err_b: the error at the b-value
     * bucket_indices: bucket_indices[k][i] is the index pointing to the correct of list index of the k-1 bin.
     *
     * returns: the approximate error reduction.
     */
    //void ComputeAppBuckets(const double * arr, const int length,  double * sum, double * sqsum, int ** b_val, double ** app_err_b, int ** bucket_indices, int max_list_length, Split *split);

    int ComputeAppBuckets(double * arr, int length, HistogramMemory * his_mem, Split * split);
    double getAppBucketsErr(int index) const;
    virtual ~Histogram();
    int *getAppBucketsList(int index);
    const double * getBucketsErr() const;
    void getBestSplit();

private:
    int n_buckets_;
public:
    Histogram();

private:
    double epsilon_;
public:
    void setNBuckets(int nBuckets);
    void setEpsilon(double epsilon);

private:
    double  buckets_err_[MAX_BUCKETS]; //an array with the error accumulated per number of buckets used
    double  app_buckets_err_[MAX_BUCKETS]; //an array with the app error accumulated per number of buckets used
    int ** buckets_list_; //for each bucket number holds an array of the left boundaries of all buckets. Used by the DP only.
    //int ** app_buckets_list; //for each bucket number holds an array of the right boundaries of all buckets.
    int  app_buckets_list_[MAX_BUCKETS][MAX_BUCKETS];
    double IntervalError(double *sum, double *sqsum, int j, int i); // computes the error of interval [j...i]
    void ComputeLists(int ** indices, const int length); //computes the left boundaries of the bucket lists
    void ComputeAppLists(int ** indices, int ** b_val, int * list_length);
    void ComputeAppList(int ** indices, int ** b_val, int list_length, int bucket, int * out_list);
    bool EnoughMemory(HistogramMemory * hismem, int arrsize); //validates that there is enough memory to perform the computation
};


#endif //C_ISOLATION_HISTOGRAM_H
