#include <iostream>
#include <ostream>
#include <ctime>
#include <random>
#include "Histogram.h"
#include "util_functions.h"
#include "Node.h"

using namespace std;

int main() {
    
    int size = 300;
    double * arr = new double[size];
    for (int i = 0; i < size; i++)
        arr[i] = i;

    int max_length = size; // Total number of segments in a level, bounded by log_{1+eps}(size) // here segments may refer to the dp subproblems at each level
    int n_buckets = 10; // Number of intervals, k
    double epsilon = 0.01;
    Split split;
    clock_t begin;
    clock_t end;
    double total_time = 0;

    Histogram m_truth(n_buckets, epsilon);
    begin = clock();
    HistogramMemory hmem(max_length, size, n_buckets); // The memory that is required for computing a histogram
    m_truth.ComputeAppBuckets(arr, size, &hmem, &split);
    
    double * s = new double[size];
    double * sq = new double[size];
    m_truth.ComputeBuckets(*&arr, size, s, sq); // Computing exact buckets using Dynamic Programming
    
    
    int ** lists = m_truth.getBucketsList();
    const double * exact_err = m_truth.getBucketsErr();
    for (int i =0; i<n_buckets; i++)
        cout<<"exact error of bucket: "<< i + 1 << " is " << exact_err[i] <<endl;
    delete []s;
    delete []sq;

    end = clock();
    total_time += (double) (end - begin) / CLOCKS_PER_SEC;
    cout<<"TOTAL TIME: "<< total_time<< endl;
    for (int k = 0; k < n_buckets; k++)
        cout<<"app error of "<<k+1<<" buckets: "<<m_truth.getAppBucketsErr(k)<<endl;
    cout<<"Printing from split "<<split.bucket_level_ + 1<<" , "<<split.bucket_error_<<endl;
    for (int k = 0; k < split.bucket_level_ + 1; k++)
        cout<<split.bucket_list_[k]<<", ";
    cout<<endl;

    /*
    double **arrs = new double *[2];
    arrs[0]=arr;
    arrs[1] = arr;
    Tree_Param * m_param = new Tree_Param(epsilon, 5, n_buckets, 100); // 5: depth, n_buckets: k intervals, 100: nSamples
    Node m_node(arrs, 1, size, 5, m_param); // 1: dim, 5: depth
    m_node.findSplit();
    */
    return 0;
}

