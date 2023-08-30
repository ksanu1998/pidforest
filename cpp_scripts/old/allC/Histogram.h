/*************************
 * This class computes the actual histograms and splits
 * **************************/
#ifndef ALLC_HISTOGRAM_H
#define ALLC_HISTOGRAM_H
#include <vector>
#include <random>

using namespace std;

class Histogram {
    public:
    Histogram(int maxBuckets, int dim);
    //vector<int> split_; // holds the indices by which we split the column
    int col_; // holds the column id that needs to be split
    double errRed_; // holds the error of the split
    void computeSplit(const vector< vector<double> > & gaps, const vector< vector<int > > & counts, vector <int> & split, mt19937 & prg_);
    
    private:
    int maxBuckets_;
    int dim_;
    double intervalError(int i, int j, const vector<double> & asum, const vector <double> & sqsum, const vector<int> tCounts);
    // computes the split of a column and returns the variance reduction
    double colSplit(const vector <double> & gapsCol, const vector <int> & countsCol, vector<int> & splitCol);
 
};

#endif //ALLC_HISTOGRAM_H