#ifndef ALLC_DATASET_H
#define ALLC_DATASET_H

#include <vector>
#include <string>

using namespace std;

class DataSet {
public:
    DataSet();
    ~DataSet();

/* Reading an array from a csv file.
 * Creates both a row based and column based representation.
 **/
    // void readData(string filename);
    std::vector<std::vector<double>> readData(string filename);
    
    vector< vector< double > > data_; // holds the data in column format.
    vector< vector< double > > dataRows_; //holds the data in row format.
    vector <bool> entropy_; // boolean vector indicator if a column has no entropy and therefore should be dropped.
    vector<double> cubeStart_; // bounding box of data lower end
    vector<double> cubeEnd_; // bounding box of data upper end
    int dim_; // dimension of the data
    int rows_; // number of rows


private:
    bool validData();
    void findBox(); // calculates the binding box around the data
};

#endif //ALLC_DATASET_H