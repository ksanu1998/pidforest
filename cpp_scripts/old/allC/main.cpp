#include <iostream>
#include <string>
#include <vector>
#include "ForestParam.h"
#include "Forest.h"
#include "DataSet.h"

using namespace std;

int main(int argc, char*argv[])
{
    //todo: proper usage
    string filename;
    if (argc == 1)
        filename = "data2.csv";
    else 
        filename = argv[1];
    //string filename = "bin_column.txt";
    DataSet data;
    data.readData(filename);

    int n_trees = 20;
    double epsilon = 0.1;
    int treeDepth = 4;
    int maxPartition = 3;
    int sampleSize = 200;
    ForestParam fParam(n_trees, epsilon, treeDepth, maxPartition, sampleSize);
    Forest mForest(fParam);

    mForest.fit(data);
    mForest.predict(data); 
    //mForest.printForest(); 
}