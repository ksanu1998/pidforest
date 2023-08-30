#include "HistogramMemory.h"

HistogramMemory::HistogramMemory(int maxLength, int arrLength, int nBuckets) : max_length_(maxLength), arr_length_(arrLength), n_buckets_(nBuckets) {
    sum_ = new double[arrLength];
    sqsum_ = new double[arrLength];
    total_count_ = new double[arrLength + 1]; 
    b_val_ = new int*[nBuckets];
    app_err_b_ = new double*[nBuckets];
    bucket_indices_ = new int*[nBuckets];
    for (int k = 0; k < nBuckets; k++){
        b_val_[k] = new int[maxLength];
        app_err_b_[k] = new double[maxLength];
        bucket_indices_[k] = new int[maxLength];
    }
}

HistogramMemory::~HistogramMemory() {
    delete[] sum_;
    delete[] sqsum_;
    delete[] total_count_;
    for (int k = 0; k < n_buckets_; k++ ){
        delete[] b_val_[k];
        delete[] app_err_b_[k];
        delete[] bucket_indices_[k];
    }
    delete[] b_val_;
    delete[] app_err_b_;
    delete[] bucket_indices_;
}