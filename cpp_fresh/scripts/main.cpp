#include <iostream>
#include <ostream>
#include <ctime>
#include <random>
#include "Histogram.h"
#include "HistogramMemory.h"
#include "Split.h"
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
    std::vector<int> values;
    std::vector<int>counts;

    Histogram m_truth(values, counts, n_buckets, epsilon);
    begin = clock();
    HistogramMemory hmem(max_length, size, n_buckets); // The memory that is required for computing a histogram
    m_truth.compute_buckets(n_buckets);
    
    return 0;
}

