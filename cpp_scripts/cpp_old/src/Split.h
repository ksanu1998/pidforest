#ifndef C_ISOLATION_SPLIT_H
#define C_ISOLATION_SPLIT_H
#define MAX_BUCKETS 20

struct Split {
    int bucket_level_;
    double bucket_error_;
    int bucket_list_[MAX_BUCKETS];
    int dim_;
};

#endif