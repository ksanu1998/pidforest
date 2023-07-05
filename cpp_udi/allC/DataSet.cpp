#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "DataSet.h"

using namespace std;

DataSet::DataSet() { }

DataSet::~DataSet() { }

/**
 * TODO: test casting into double
 **/
void DataSet::readData( string filename) {
    cout<<"Reading file "<<filename<<endl;
    ifstream in(filename);
    try {
        if (in) {
            string line;
            while (getline(in, line)) {
                stringstream sep(line);
                string field;
                dataRows_.push_back(vector<double>()); // if casting into double fails will throw an exception
                while (getline(sep, field, ',')) {
                    dataRows_.back().push_back(stod(field));
                }
            }
        }
    }
    catch (const exception& e) { 
     cerr << e.what();
     exit(1); 
}
    rows_ = dataRows_.size();
    dim_ = dataRows_[0].size();

    if (!validData()) {
        cout<<"Error in data\n";
        //return;
        exit(-1); // todo: throw expection
    }
    else {
        cout<<"number of rows is "<<rows_<<endl;
        cout<<"size of row is "<<dim_<<endl;
    }
    // creating the columnar form.
    for (int i = 0; i < dim_; i++) {
        data_.push_back(vector<double>(rows_));
        for (int j = 0; j < rows_; j++)
            data_[i][j] = dataRows_[j][i]; 
    }
    // calculating the bounding box
    findBox();
}


//checks that data is not empty and that all rows are of same size
bool DataSet::validData() {
    if (dataRows_.empty())
        return false;
    int dim = dataRows_[0].size();
    for (vector<double> row : dataRows_) {
        if (row.size() != dim)
            return false;
    }
    return true;
}

/** caclualates the bounding box around the the data and marks which columns have no entropy
 */
void DataSet::findBox() {
    double min, max, min2, max2;
    bool entropy;
    for (int col = 0; col < dim_; col++) {
        min = max = data_[col][0];
        double curr;
        entropy = true;
        for (int row = 0; row < rows_; row++) {
            curr = data_[col][row];
            if (curr < min)
                min = curr;
            else if (curr > max)
                max = curr;
        }
        if (min == max) { // no entroy in the column
            entropy = false;
            cubeStart_.push_back(min - 0.5);
            cubeEnd_.push_back(max + 0.5);
            entropy_.push_back(entropy);
            continue;
        }
        min2 = max;
        max2 = min;
        for (int row = 0; row < rows_; row++) {
            curr = data_[col][row];
            if ((curr < min2) and (curr > min))
                min2 = curr;
            if ((curr > max2) and (curr < max))
                max2 = curr;
        }
        cubeStart_.push_back((3 * min - min2) / 2);
        cubeEnd_.push_back((3 * max - max2) / 2);
        entropy_.push_back(entropy);
    } // end for cols
}