#ifndef C_ISOLATION_HISTOGRAMMEMORY_H
#define C_ISOLATION_HISTOGRAMMEMORY_H


class HistogramMemory {
public:
    double * sum_;
    double * sqsum_;
    double * total_count_;
    int ** b_val_;
    double ** app_err_b_;

    HistogramMemory(int maxLength, int arrLength, int nBuckets);

    virtual ~HistogramMemory();
    int ** bucket_indices_;
    const int max_length_; //the length of the lists which is a function of the array length and epsilon.
    const int n_buckets_;
    const int arr_length_;
};


#endif //C_ISOLATION_HISTOGRAMMEMORY_H
