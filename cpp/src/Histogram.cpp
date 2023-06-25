#include "Histogram.h"
#include "util_functions.h"
#include <ostream>
#include <string.h>
#include <iostream>
#include <math.h>

using namespace std;

Histogram::Histogram(const int nBuckets, const double epsilon) : n_buckets_(nBuckets), epsilon_(epsilon) {
    if (n_buckets_ > MAX_BUCKETS) {
        std::cerr<<"Number of buckets exceeds maximum of "<<MAX_BUCKETS<<endl;
        exit(1);
    }
}

double  Histogram::getBucketsErr(int index) const {
    return buckets_err_[index];
}

Histogram::~Histogram() {
  /*  for (int k = 0; k < n_buckets_; k++) {
        delete[] buckets_list_[k];
    }
    delete[] buckets_list_; */
}

/*** computes the optimal histogram using a dynamic program. To be used to verify correctness. Not the right interface.
 * 
    TODO: old code, needs recoding or removal
*/
void Histogram::ComputeBuckets(const double * arr, const int length, double * sum, double * sqsum) {
    //std::cout<<"In Compute Buckets"<<endl;
    memset(sum, 0, length * sizeof(double));
    memset(sqsum, 0, length * sizeof(double));
    sum[0] = arr[0];
    sqsum[0] = arr[0] * arr[0];
    buckets_list_ = new int *[n_buckets_];
    /*for (int k = 0; k < n_buckets_; k++) {
        buckets_list_[k] = new int[length]; //used for dynamic program. Should be recoded or removed
        buckets_list_[k][0] = 0;
    }*/
    // tmp_buckets[k][j] stores the error when using at most k+1 buckets for [0....j].
    double **tmp_buckets = new double*[n_buckets_];
    //holds the list of left bucket boundaries
    int **bucket_indices = new int*[n_buckets_];
    for (int k = 0; k < n_buckets_; k++) {
        tmp_buckets[k] = new double[length];
        bucket_indices[k] = new int[length];
        tmp_buckets[k][0] = 0;
        bucket_indices[k][0] = 0;
    }
    for (int i=1; i < length; i++) {
        sum[i] = sum[i - 1] + arr[i];
        sqsum[i] = sqsum[i - 1] + (arr[i] * arr[i]);
        tmp_buckets[0][i] = sqsum[i] - ((sum[i] * sum[i]) / (double)(i + 1)); //error for 1 bucket
        bucket_indices[0][i] = 0;
    }
    for (int i = 1; i < length; i++){
        for (int k = 1; k < n_buckets_; k++ ) {
            tmp_buckets[k][i] = tmp_buckets[k-1][i-1]; // the case the last item is a new bucket
            bucket_indices[k][i] = i;
            for (int j = 0; j < i - 1; j++){
                double tmp_error = tmp_buckets[k-1][j] + IntervalError(sum, sqsum, j + 1, i);
                if (tmp_error < tmp_buckets[k][i]) {
                    tmp_buckets[k][i] = tmp_error;
                    bucket_indices[k][i] = j + 1;
                }
            }
        }
    }
    for (int k = 0; k < n_buckets_; k++ ) {
        buckets_err_[k] = tmp_buckets[k][length - 1];
    }
    ComputeLists(bucket_indices, length);
    for (int k = 0; k < n_buckets_; k++) {
        delete[] tmp_buckets[k];
        delete[] bucket_indices[k];
    }
    delete[] tmp_buckets;
    delete[] bucket_indices;
}
/** Computes the actual split, used by DP only **/
void Histogram::ComputeLists(int **indices, const int length) {
    int b ,tmp,loc;
    for (int k = 0; k < n_buckets_; k++) {
        b = k;
        loc   = length - 1;
        while (b > 0) {
            // std::cout << k << " " << b << endl;
            tmp = indices[b][loc];
            // buckets_list_[k][b] = tmp; // this line was throwing segfault as buckets_list_ is not in scope of this function!
            // cout << "hist here" << endl;
            if (tmp > 0)
                tmp = tmp - 1;
            loc = indices[b - 1][tmp];
            b = b - 1;
        }
    }
}

// computes the error of a bucket  in indexes j,...,i
double Histogram::IntervalError(double *sum, double *sqsum, int j, int i) {
    double lb_s, lb_sq;
    if (i <= j) { return 0; }
    if (j == 0) { lb_s = 0; lb_sq = 0; }
        else {
            lb_s = sum[j-1]; 
            lb_sq = sqsum[j-1]; 
        }      
    double avg = (sum[i] - lb_s)/(i - j + 1);
    if (utils::isEqual(avg, 0))
        avg = 0;
    double ans = (sqsum[i] - lb_sq) - ((i - j + 1) * avg * avg);
    if (utils::isEqual(ans,0))
        ans = 0;
    if (ans < 0) {
        std::cout<<"ERROR in InteervalError "<<j<<i<<"..."<<sqsum[i]<<","<<lb_sq<<","<<sum[i]<<","<<lb_s<<","<<avg<<","<<ans<<endl;
    }
    return (sqsum[i] - lb_sq) - ((i - j + 1) * avg * avg);
}

int **Histogram::getBucketsList() const {
    return buckets_list_;
}

/**
 * Computing buckets approximately, following the paper "Approximation and Streaming Algorithms for Histogram Construction Problems" 
 * by Guha et. al.
 * 
 *  arr: the array with the data
 *  length: length of arr
 *  his_mem: the object that holds all the memory needed for histogram. Allocated and freed outside the class.
 *  split: The object that will hold the ouput of computations. Allocated and freed outside the class.
 * 
 **/
//todo: Optimize perf by putting sum and sqsum on the same array with even/odd placements
int Histogram::ComputeAppBuckets(double * arr, int length, HistogramMemory * his_mem, Split * split){
    if (!EnoughMemory(his_mem, length)){
        cout<<"ERROR: Not enought memory in HistogramMemory for computation";
        return 0;
        //todo Exception handling...
    }
    double * sum = his_mem->sum_;
    double * sq_sum = his_mem->sqsum_;
    //double * sq_sum = new double[his_mem->max_length_];
    int ** b_val = his_mem->b_val_; //buckets are [a,b]
    double ** app_err_b = his_mem->app_err_b_; // the approximate error at b
    int ** bucket_indices = his_mem->bucket_indices_;
    int max_length = his_mem->max_length_;
    int list_length[MAX_BUCKETS];
    double curr_err_a[MAX_BUCKETS];

/* Initializing variables */
    for (int k = 0; k < n_buckets_; k++) {
        list_length[k] = 0;
        b_val[k][0] = 0;
        app_err_b[k][0] = 0;
        bucket_indices[k] = new int[max_length];
        bucket_indices[k][0] = 0;
        curr_err_a[k] = 0;
        app_buckets_list_[k][0] = 0;
    }
    sum[0] = arr[0];
    sq_sum[0] = arr[0] * arr[0];
    //std::cout<<"sqsum[0] is: "<<sq_sum[0]<<endl;
    double tmp = 0;

/* Main loop on array items */
    for (int j = 1; j < length; j++) {
        //preparing the one bucket case. That's the recursion base.
        sum[j] = sum[j-1] + arr[j];
        sq_sum[j] = sq_sum[j - 1] + (arr[j] * arr[j]);
        tmp = sq_sum[j] - ((sum[j] * sum[j]) / (double)(j + 1)); //error for 1 bucket
        if (tmp < (1 + epsilon_) * curr_err_a[0]) { //segment ends in j
            b_val[0][list_length[0]] = j;
        } else { // open a new segment
            list_length[0]++;
            b_val[0][list_length[0]] = j;
            curr_err_a[0] = tmp;
            bucket_indices[0][list_length[0]] = 0;
        }
        app_err_b[0][list_length[0]] = tmp;

        // computing for buckets 2...
        for (int k = 1; k < n_buckets_; k++) {
            double tmp_err = sq_sum[j];
            //if (tmp_err < 0) std::cout<<"ERROR0"<<endl;
            int curr_index = 0;
            for (int i = 0; i <= list_length[k - 1]; i++) {
                int b  = b_val[k-1][i];
                if (tmp_err > app_err_b[k - 1][i] + IntervalError(sum, sq_sum, b + 1, j)) {
                    tmp_err = app_err_b[k - 1][i] + IntervalError(sum, sq_sum, b + 1, j);
                    curr_index = i;
                }
            }
            if ((k < n_buckets_ - 1 ) and (tmp_err > (1 + epsilon_) * curr_err_a[k])) { //open a new segment
                list_length[k]++;
                int l = list_length[k];
                b_val[k][l] = j;
                app_err_b[k][l] = tmp_err;
                bucket_indices[k][l] = curr_index;
                curr_err_a[k] = tmp_err;
            }
            else {
                int l = list_length[k];
                b_val[k][l] = j;
                app_err_b[k][l] = tmp_err;
                bucket_indices[k][l] = curr_index;
            }
        }
    }
    for (int k = 0; k < n_buckets_; k++){
        app_buckets_err_[k] = app_err_b[k][list_length[k]];
    }

   // Now find the best split:
   double var_red = 0;
   int best_bucket = 0;
  if (utils::isEqual(app_buckets_err_[0], 0)) {
      split->bucket_error_ = 0;
      split->bucket_level_ = 0;
  }
  else {
      for (int k = 1; k < n_buckets_; k++) {
          double err_red = (app_buckets_err_[0] - app_buckets_err_[k]) / app_buckets_err_[0];
          if (err_red > var_red) {
              var_red = err_red;
              best_bucket = k;
          }
      }
          split->bucket_error_ = var_red;
          split->bucket_level_ = best_bucket;
  }
  ComputeAppList(bucket_indices, b_val, list_length[best_bucket], best_bucket, split->bucket_list_);
  return 1;
}

//validates that there is enough memory to perform the computation
bool Histogram::EnoughMemory(HistogramMemory * hismem, int arrsize){
    bool sufficient_memory = true;
    if (arrsize < hismem->arr_length_)
        sufficient_memory = false;
    if (this->n_buckets_ <  hismem->n_buckets_)
        sufficient_memory  = false;
    double lb;
    if (utils::isEqual(this->epsilon_, 0)) lb = arrsize;
    else lb = log(arrsize)/log(1 + this->epsilon_);
    if (lb > arrsize) lb = arrsize; 
    if (hismem->max_length_ < int(lb))
        sufficient_memory = false;
    return sufficient_memory;
}

double Histogram::getAppBucketsErr(int index) const {
    return app_buckets_err_[index];
}

void Histogram::ComputeAppLists(int **indices, int ** b_val, int * list_length) {
    for (int k = 1; k < n_buckets_; k++) {
        int curr_level = k;
        int loc = list_length[curr_level];
        while (curr_level > 0) {
            app_buckets_list_[k][curr_level] = b_val[curr_level][loc];
            loc = indices[curr_level][loc];
            curr_level--;
        }
        app_buckets_list_[k][curr_level] = b_val[curr_level][loc];
    }
}

int* Histogram::getAppBucketsList(int index){
    return app_buckets_list_[index];
}

const double *Histogram::getBucketsErr() const {
    return buckets_err_;
}

void Histogram::ComputeAppList(int **indices, int **b_val, int list_length, int bucket, int *out_list) {
    int curr_level = bucket;
    int loc = list_length;
    while (curr_level > 0) {
        out_list[curr_level] = b_val[curr_level][loc];
        loc = indices[curr_level][loc];
        curr_level--;
    }
    out_list[curr_level] = b_val[curr_level][loc];
}

void Histogram::setNBuckets(int nBuckets) {
    n_buckets_ = nBuckets;
}

void Histogram::setEpsilon(double epsilon) {
    epsilon_ = epsilon;
}

Histogram::Histogram() {}
