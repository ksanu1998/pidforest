#ifndef ALLC_FOREST_H
#define ALLC_FOREST_H

#include<vector>
#include "ForestParam.h"
#include "Tree.h"
#include "DataSet.h"

using namespace std;

class Forest {
public: 
    Forest(const ForestParam & forestParam);
    ~Forest();
    void fit(const DataSet & dataSet);
    void predict(const DataSet & dataSet);
    void printForest();

private:
    ForestParam forestParam_; // global parameters of the forest
    Tree** trees_; // pointer to array of trees
    vector<double> scores_; // holds the scores the latest predicted data set.
    void indexSort(const vector<double> & scores, vector<int> & index);
};

#endif //ALLC_FOREST_H