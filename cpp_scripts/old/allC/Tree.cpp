#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include "Tree.h"
#include "math.h"
#include "Histogram.h"
#include "util.h"

using namespace std;

Tree::Tree(int id, ForestParam & forestParam) :  forestParam_(forestParam), id_(id), root_(0)
{ 
    mPrg_.seed(time(0) + 2*id_);
    root_.setTree(this); 
} 

//Tree::~Tree() {
//}

/**
 * data: the entire data set.
 * Will compute a tree by first sampling a subset and then computing the tree
 * */
void Tree::computeTree(const DataSet * data) {
    dataSet_ = data;
    vector<int> indices;
    getSample(indices);
    root_.computeTree(indices, dataSet_->cubeStart_, dataSet_->cubeEnd_);
}

void Tree::predict(const DataSet * data, vector<double> & scores) {
    cout<<"In Predict\n";
    vector <int> indices(data->rows_);
    for (int i = 0; i < indices.size(); i++)
        indices[i] = i;
    root_.predict(data, scores, indices);
}


void Tree::getSample(vector <int> & indices) {
    int sampleSize = (dataSet_->rows_ >= forestParam_.n_samples_) ? forestParam_.n_samples_ : dataSet_->rows_;
    indices.reserve(sampleSize);
    if (sampleSize < dataSet_->rows_) {
        for (int i = 0; i < sampleSize; i++)
            indices.push_back(mPrg_() % dataSet_->rows_);
    }
    else {
        for (int i = 0; i < sampleSize; i++)
            indices.push_back(i);
    }
    sort(indices.begin(), indices.end());
    indices.erase( unique( indices.begin(), indices.end() ), indices.end() ); //indices are sorted and without duplicates
}

void Tree::printTree() {

    cout<<"**********************************"<<endl<<"Printing Tree"<<endl<<"**************************"<<endl;
    root_.printSubTree();
}

/**************************************************
 * Node is an internal class of Tree. The Tree is composed of nodes
 * ************************************************/

Tree::Node::Node(int Depth): depth_(Depth) { 
    mTree_ = NULL;
}

void Tree::Node::setTree(Tree * parent) { mTree_ = parent;}

/** filters the indices that are in the cube 
 * gets a vector of oldIndices and places those that are in the cube in newIndices
 */
void Tree::Node::filterIndices(const DataSet * data, vector<int> & newIndices, const vector<int> & oldIndices) {
    bool inside;
    double curr;
    for (int index : oldIndices) {
        inside = true;
        for (int i = 0; i < data->dim_; i++) {
            curr = data->dataRows_[index][i];
            if ((curr < cubeLB_[i]) or (curr > cubeUB_[i]))
                inside = false;
        }
        if (inside) {
            newIndices.push_back(index);
        }
    }
}

/** filters the indices that are in the cube 
 * gets a vector of oldIndices and places those that are in the cube in newIndices
 **/
void Tree::Node::filterIndices(vector<int> & oldIndices) {
    bool inside;
    double curr;
    for (int index : oldIndices) {
        inside = true;
        for (int i = 0; i < mTree_->dataSet_->dim_; i++) {
            curr = mTree_->dataSet_->dataRows_[index][i];
            if ((curr < cubeLB_[i]) or (curr > cubeUB_[i]))
                inside = false;
        }
        if (inside) {
            fitIndices_.push_back(index);
        }
    }
}

void Tree::Node::computeVolume() {
    volume_ = 0;
    double curr = 0;
    for (int i = 0; i < mTree_->dataSet_->dim_; i++) {
        curr = cubeUB_[i] - cubeLB_[i];
        if (curr > 0) {
            volume_ += log(curr);
        }
    }
}

/** 
 * Receives a vector of indices 
 * Creates the vectors of gaps and counts that would then be used to compute histograms.
 * This is a perf bottleneck and should be performant 
 * TODO: optimize performance by not computing uniqueCount and using counts as index over column 
 * **/
void Tree::Node::createGaps(const vector<int> & indices, vector< vector<double> > & gaps, vector< vector<int> > & counts) {
    const DataSet * data = mTree_->dataSet_;
    vector <double> column(indices.size());
    vector <double> uniqueColumn;
    
    int colSize;
    for (int col = 0; col < data->dim_; col++) {
        for (int i = 0; i < indices.size(); i++) { // get the data from the indices vector
            column[i] = data->data_[col][indices[i]];
        }
        sort(column.begin(), column.end()); // sort the column
        uniqueColumn = column;
        uniqueColumn.erase( unique( uniqueColumn.begin(), uniqueColumn.end() ), uniqueColumn.end() );
        if (uniqueColumn.size() <= 1) { // there is no entropy in the column
            gaps[col].push_back(0);
            counts[col].push_back(column.size());
            continue;  
        } 
        gaps[col].resize(uniqueColumn.size());
        //fill in gaps
        gaps[col][0] = (uniqueColumn[0] + uniqueColumn[1])/2 - cubeLB_[col];
        colSize = uniqueColumn.size();
        gaps[col][colSize - 1] = cubeUB_[col] - (uniqueColumn[colSize - 1] + uniqueColumn[colSize - 2])/2;
        for (int i = 1; i < colSize - 1; i++)
            gaps[col][i] = (uniqueColumn[i+1] -  uniqueColumn[i-1])/2;

        //compute counts:
        int i = 0;
        while(i < column.size()) {
            int currCount = 0;
            double currVal = column[i];
            while (column[i] == currVal) {
                currCount++;
                i++;
            }
            counts[col].push_back(currCount);
        }
    }
}

/******************************************************************* 
 * Initiates the computation of the tree. The cubes bounds are the ones of the new node. 
     * Need to compute the rest of the tree, will do the following steps:
     * 1. filter the indices and get the indices that are within the cube
     * 2. Compute an array of gaps per dimension. These are the values the histogram works on.
     * 3. Send this array of gaps to the Histogram class that would compute a split.
     * 4. Split the cube accordingly and create the children nodes.
 *****************************************************************/
void Tree::Node::computeTree(vector <int> & indices, const vector <double> & cubeLow, const vector <double> & cubeHigh) {
    //cout<<"In ComputeTree!"<<endl;   
    cubeLB_ = cubeLow;
    cubeUB_ = cubeHigh;
    //vector <int> myIndices; // the indices of this current node.
    //filterIndices(myIndices, indices); // filter the indices for the current node.
    filterIndices(indices);
    computeVolume();
    if (fitIndices_.size() <= 1) { //isolated a data point. No need to split further.
        density_ = 0 - volume_;
        return;
    }
    density_ = log(fitIndices_.size()) - volume_;
    if (depth_ >= mTree_->forestParam_.depth_) // reached the maximum depth allowed for a tree. No need to split further
        return;

    // Start splitting the indices across Children
    //printCube();
    vector < vector <double> > gaps( mTree_->dataSet_->dim_); // The vector of gaps already with the right size
    vector < vector <int> > counts( mTree_->dataSet_->dim_); // The vector of counts
    createGaps(fitIndices_, gaps, counts);

    // find split
   // cout<<"Finding next Split"<<endl;
    Histogram hist(mTree_->forestParam_.n_bucket_, mTree_->dataSet_->dim_);
    vector<int> split;
    hist.computeSplit(gaps, counts, split, mTree_->mPrg_);
    if (split.empty()) {
        return;
    }
    //cout<<"split size is "<<split.size()<<endl;
    
    double b1 = cubeUB_[hist.col_]; //upper bound of the cube
    double b2 = b1;
    int i = gaps[hist.col_].size()-1;
    int j = 0; //will iterate through the splits

    while (j < split.size()) {
        while (i >= split[j]) {
            b1 -= gaps[hist.col_][i];
            i--;
        }
        vector<double> newCubeU = cubeUB_;
        vector<double> newCubeL = cubeLB_;
       // cout<<"cube "<<j<<" is: "<<b1<<","<<b2<<endl;
        newCubeU[hist.col_] = b2;
        newCubeL[hist.col_] = b1;
       
        Node * new_newChild = new Node(depth_ + 1);
        new_newChild->setTree(mTree_);
        new_newChild->computeTree(fitIndices_, newCubeL, newCubeU);
        m_children_.push_back(new_newChild);
        

      /*  Node newChild = Node(depth_ + 1);
        newChild.setTree(mTree_);
        newChild.computeTree(fitIndices_, newCubeL, newCubeU);
        children_.push_back(newChild);*/
        b2 = b1;
        j++;
    }
    // find the subcubes and create the children with the subcube and 'myIndices', and the pointer to mTree.
    // childre.computeTree()
}

void Tree::Node::printCube() {
    cout<<"Lower bound of cube is: ";
    for (int j = 0; j < cubeLB_.size(); j++)
         cout<<cubeLB_[j]<<", ";
    cout<<"\nUpper bound of cube is: ";
    for (int j = 0; j < cubeUB_.size(); j++)
         cout<<cubeUB_[j]<<", ";
    cout<<"***"<<endl;
}

void Tree::Node::printSubTree() {
    cout<<"Printing Node"<<endl;
    cout<<"depth: "<<depth_<<endl;
    cout<<"volume: "<<volume_<<endl;
    cout<<"density: "<<density_<<endl;
    cout<<"Num of points: "<<fitIndices_.size()<<endl;
    printCube();
    cout<<"num of children is: "<<m_children_.size()<<endl;
    cout<<"***************************\n";
    if (m_children_.size() > 0) {
        for (auto child: m_children_)
            child->printSubTree();
    }
}

Tree::Node::~Node() {
    if (m_children_.empty()) 
        return;
    for (auto tmp : m_children_)
        delete(tmp);
}

void Tree::Node::predict(const DataSet * data, vector<double> & scores, const vector<int> & indices) {
    vector<int> m_indices;
    filterIndices(data, m_indices, indices);
    if (m_children_.size() > 0) {
        for (int i = 0; i < m_children_.size(); i++)
            m_children_[i]->predict(data, scores, m_indices);
    }
    else {
        for (int i = 0; i < m_indices.size(); i++)
            scores[m_indices[i]] = density_;
    }
}