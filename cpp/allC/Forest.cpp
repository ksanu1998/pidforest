#include "Forest.h"
#include <vector>
#include <iostream>

Forest::Forest(const ForestParam& forestParam) : forestParam_(forestParam) { 
    trees_ = new Tree*[forestParam_.n_trees_];
    for (int i=0; i < forestParam_.n_trees_ ; i++) {
        trees_[i] = new Tree(i, forestParam_);
    }
} 

Forest::~Forest() { 
    for (int i=0; i < forestParam_.n_trees_ ; i++) {
        delete trees_[i];
    }
    delete [] trees_;
}

void Forest::fit(const DataSet & dataSet) { 
    // will ask all trees to compute their trees. Here we will span multi threads.
    for (int i=0; i < forestParam_.n_trees_ ; i++) {
         trees_[i]->computeTree(& dataSet);
    }
}

void Forest::predict(const DataSet & dataSet) {
    //vector <double> scores(dataSet.rows_); 
    scores_.resize(dataSet.rows_);
    trees_[0]->predict(& dataSet,  scores_);
    cout<<"Scores\n";
    for (int i = 0; i < scores_.size(); i++)
        cout<<scores_[i]<<", ";
    cout<<endl;
    vector <int> s_anom(dataSet.rows_);
    std::iota(s_anom.begin(),s_anom.end(),0); //Initializing
   // auto compare = [](int i, int j) {return score}
    //std::sort(s_anom.begin(),s_anom.end());
    indexSort(scores_, s_anom);
    int anomaly = s_anom[0];
    cout<<"largest index is: "<<anomaly<<endl;
    cout<<"the point is  "<<dataSet.dataRows_[anomaly][0]<<","<<dataSet.dataRows_[anomaly][1]<<endl;
    
}

void Forest::printForest() {
    cout<<"*********  Printing Forest *************\n";
    for (int i=0; i < forestParam_.n_trees_; i++) {
        cout<<"<<<<<<< Printing Tree: "<<i<<" >>>>>>>>>>>>>>\n";
        trees_[i]->printTree();
    }
}

void Forest::indexSort(const vector<double> & scores, vector<int> & idx) { 
    stable_sort(idx.begin(), idx.end(), 
        [&](int i, int j) {return scores[i] < scores[j];});
}
