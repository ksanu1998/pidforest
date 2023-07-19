#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "Forest.h"
#include "Node.h"
#include "TimeSeries.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace std;

int main() {
    // Load your data and set up variables
    // string dataset = "nyc_taxi";
    int n_samples = 100;
    int max_depth = 10;
    int n_trees = 50;
    int max_samples = n_samples;
    int max_buckets = 3;
    double epsilon = 0.1;
    double sample_axis = 1.0;
    double threshold = 0.0;
    
    // Create the forest
    unordered_map<string, double> kwargs = {
        {"max_depth", max_depth},
        {"n_trees", n_trees},
        {"max_samples", max_samples},
        {"max_buckets", max_buckets},
        {"epsilon", epsilon},
        {"sample_axis", sample_axis},
        {"threshold", threshold}
    };
    Forest mForest(kwargs);
    
    string dataset = "nyc_taxi";  // Set your dataset name here

    // Read the CSV file
    string filename = "../data/numenta/" + dataset + ".csv";
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return 1;
    }

    vector<double> value;
    vector<int> label;
    string line;

    // Skip the first line (header)
    getline(file, line);

    // Read the file line by line
    while (getline(file, line)) {
        stringstream ss(line);
        string timestamp, val, anomaly_score, lab;
        
        if (getline(ss, timestamp, ',') && getline(ss, val, ',') &&
            getline(ss, anomaly_score, ',') && getline(ss, lab, ',')) {
            try {
                value.push_back(stod(val));
                label.push_back(stoi(lab));
            } catch (const exception& e) {
                cerr << "Error converting value or label: " << e.what() << endl;
                return 1;
            }
        }
    }

    // Shingle and transpose the data
    int shingleSize = 10;
    vector<vector<double>> X = TimeSeries::shingle(value, shingleSize);
    X = TimeSeries::transpose(X);

    // Print the shingled and transposed data
    cout << "X:" << endl;
    for (const auto& row : X) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }


    std::cout << "\n\n******" << dataset << " trial 1*******\n\n";
    
    std::cout << "\n******Our Algo*******\n";

    vector<vector<double>> pts;  // Provide your data here
    pts = X;
    std::cout << "\n******FIT*******\n";
    mForest.fit(pts);
    std::cout << "\n******FIT_DONE*******\n";
    
    // Make predictions
    vector<int> indices;
    unordered_map<int, vector<double>> outliers;
    unordered_map<int, vector<double>> scores;
    unordered_map<int, double> pst;
    vector<double> our_scores;
    double err = 0.1;
    double pct = 0.0;
    std::cout << "\n******PREDICT*******\n";
    
    tie(indices, outliers, scores, pst, our_scores) = mForest.predict(pts, err, pct);
    
    // Output the results
    cout << "Indices: ";
    for (int idx : indices) {
        cout << idx << " ";
    }
    cout << endl;
    
    // Print the outliers and their scores
    for (const auto& entry : outliers) {
        int idx = entry.first;
        const vector<double>& outlierValues = entry.second;
        cout << "Outlier at index " << idx << ": ";
        for (double val : outlierValues) {
            cout << val << " ";
        }
        cout << endl;
    }
    
    // Print the scores for each index
    for (const auto& entry : scores) {
        int idx = entry.first;
        const vector<double>& scoreValues = entry.second;
        cout << "Scores for index " << idx << ": ";
        for (double val : scoreValues) {
            cout << val << " ";
        }
        cout << endl;
    }
    
    // Print the PST values for each index
    for (const auto& entry : pst) {
        int idx = entry.first;
        double pstValue = entry.second;
        cout << "PST for index " << idx << ": " << pstValue << endl;
    }
    
    // Print the our_scores values
    cout << "Our Scores: ";
    for (double val : our_scores) {
        cout << val << " ";
    }
    cout << endl;
    
    return 0;
}
