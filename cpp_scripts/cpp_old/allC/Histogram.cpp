#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <float.h>
#include "Histogram.h"
#include "util.h"

Histogram::Histogram(int maxBuckets, int dim) : maxBuckets_(maxBuckets), dim_(dim) { }

/**
 * Gets vectors of gaps and vectors of counts, and computes the best split. 
 * fills the members: col_ the col number that is split. VarRed: the relative reduction in error.
 * Split: A vector with the indices of the left boundaries of the partitions in reverse order (so last item is 0)
 *  * ******/
void Histogram::computeSplit(const vector< vector<double> > & gaps, const vector< vector<int > > & counts, vector <int> & split, mt19937 & prg) {
    vector<double> varReds(dim_); // will hold the variance reduction obtained
    //vector< vector <int> > splits(dim_, vector<int>(maxBuckets_)); // will hold the actual split
    vector< vector <int> > splits(dim_);
    vector<double>newGaps(gaps[0].size());
    for (int col = 0; col < dim_; col++) {
        for (int j = 0; j < newGaps.size(); j++)
            newGaps[j] = gaps[col][j] / counts[col][j]; 
        //varReds[col] = colSplit(gaps[col], counts[col], splits[col]);
       // printVector(newGaps, "hist column of new gaps");
       // printVector(counts[col], "counts vector");
        if (newGaps.size() <= 1) {
            //cout<<"nothing to split in column: "<<col<<endl;;
            varReds[col] = 0;
            continue;
        }
        varReds[col] = colSplit(newGaps, counts[col], splits[col]);
        //printVector(splits[col], "Split");
    }
    double sumerr = 0;
    for (double er : varReds)
        sumerr += er;
    if (sumerr < 0.0001) // no point of splitting.
        {
            //cout<<"no error reduction possible, not splitting"<<endl;
            split.clear();
            return;
        }

    uniform_real_distribution<double> dis(0.0, sumerr);
    double sam = dis(prg);
    int i = 0;
    double tmpsum = varReds[i];
    while (sam > tmpsum) {
        i++;
        tmpsum += varReds[i];
    }
    //cout<<"sampled "<<i<<endl;
    col_ = i;
    errRed_ = varReds[i];
    split = splits[i]; /* TODO: Avoid this new allocation */
}

/**
 * computes the histogram with the simple opt DP algorithm.
 * returns: the error
 * returns: places the indices of the left boundries of the partition in 'split'
 * 
 * TODO: vsum, sqsum, totalCounds, hiseErrror, indices, should all be allocated once in 'computeSplit' and passed along to be reused.
 * This will save a factor of dim_ in the number of memory allocations. 
 * 
 * TODO: implement the approximation algorithm instead of the DP. 
 * ************/
double Histogram::colSplit(const vector<double> & gaps, const vector<int> & counts, vector<int> & split) {
    int size = gaps.size();
    int bucketNum;
    (maxBuckets_ < size) ? bucketNum = maxBuckets_ : bucketNum = size;
    vector<double> vsum(size);
    vector<double> sqsum(size);
    vector<int> totalCounts(size); // holds the total counts of elements
    vector< vector<double> > histError(size, vector<double>(bucketNum)); //histError[i][k] := optimal error for array 0...i using k+1 partitions
    vector< vector<int> > indices(size, vector<int>(bucketNum));
    double tmperr;

    vsum[0] = gaps[0] * counts[0]; // sum[i] holds the sum of 0...i
    sqsum[0] = gaps[0] * gaps[0] * counts[0]; // sqsum[i] holds the square sum of 0...i
    totalCounts[0] = counts[0];
  
    for (int i = 1; i < size; i++) {
        vsum[i] = vsum[i-1] + gaps[i] * counts[i];
        sqsum[i] = sqsum[i-1] + gaps[i] * gaps[i] * counts[i];
        totalCounts[i] = totalCounts[i-1] + counts[i];
        
    }
    for (int i = 0; i < size; i++) {
        histError[i][0] = intervalError(0,i,vsum, sqsum, totalCounts);
        for (int k = 1; k < bucketNum; k ++)
            histError[i][k] = DBL_MAX; /** TODO: change initialization of histError to DBL_MAX upon constructin **/
    }
    for (int i = 0; i < size; i++)
        for (int k = 1; k < bucketNum; k++)
            for (int j = 0;  j < i; j++) {
                tmperr = histError[j][k-1] + intervalError(j+1,i,vsum, sqsum, totalCounts);
                if (tmperr < histError[i][k]) {
                    indices[i][k] = j + 1;
                    histError[i][k] = tmperr; 
                }
            }
    double minError = histError[size - 1][0];
    double ErrorOneBucket = minError;
    int numPartitions = 1; // holds how many buckets are actually needed. Typicallly will be max_Buckets - 1.
    for (int i = 1; i < bucketNum; i++)
        if (minError > histError[size - 1][i]) {
            minError = histError[size - 1][i];
            numPartitions = i;
        }
    int currItem = size - 1;
    int currBucket = numPartitions;
    int tmpBucket;
    while (currItem > 0) {
        tmpBucket = indices[currItem][currBucket];
        split.push_back(tmpBucket);
        currItem = tmpBucket - 1;
        currBucket = currBucket - 1;
    }
    if (split.back()!= 0)
        split.push_back(0);
    return (ErrorOneBucket - minError) / ErrorOneBucket;
} 

double Histogram::intervalError(int i, int j, const vector <double> & asum, const vector <double> & sqsum, const vector<int> tCounts) {
    double tmpsum, tmpsq;
    int count;
    if (i==0) {
        tmpsum = asum[j];
        tmpsq = sqsum[j];
        count = tCounts[j];
    }
    else {
        tmpsum = asum[j] - asum[i-1];
        tmpsq = sqsum[j] - sqsum[i-1];
        count = tCounts[j] - tCounts[i-1]; 
    }
    return tmpsq - (tmpsum * tmpsum)/count;
}