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

double calculatePrecisionAll(const std::vector<double>& precision_our, const std::vector<double>& recall_our) {
    double max_precision_all = 0.0;

    for (size_t i = 0; i < precision_our.size(); ++i) {
        double precision = precision_our[i];
        double recall = recall_our[i];

        // Calculate precision_all for the current values of precision_our and recall_our
        double current_precision_all = 2.0 * precision * recall / (precision + recall);

        // Update max_precision_all if the current value is greater
        if (current_precision_all > max_precision_all) {
            max_precision_all = current_precision_all;
        }
    }

    return max_precision_all;
}

// Define a function to calculate precision, recall, and thresholds
void calculatePrecisionRecallThresholds(const std::vector<int>& y, const std::vector<double>& our_scores) {
    int t1 = y.size(); // Assuming t1 is the size of the y vector

    // Create arrays for precision, recall, and thresholds
    std::vector<double> precision_our;
    std::vector<double> recall_our;
    std::vector<double> thresholds_our;

    // Iterate through different threshold values
    for (double threshold = 0.0; threshold <= 1.0; threshold += 0.01) {
        // Calculate true positives, false positives, false negatives
        int true_positives = 0;
        int false_positives = 0;
        int false_negatives = 0;

        for (int i = 0; i < t1; ++i) {
            if (our_scores[i] <= -threshold && y[i] == 1) {
                true_positives++;
            } else if (our_scores[i] <= -threshold && y[i] == 0) {
                false_positives++;
            } else if (our_scores[i] > -threshold && y[i] == 1) {
                false_negatives++;
            }
        }

        // Calculate precision and recall
        double precision = static_cast<double>(true_positives) / (true_positives + false_positives);
        double recall = static_cast<double>(true_positives) / (true_positives + false_negatives);

        // Store precision, recall, and threshold values
        precision_our.push_back(precision);
        recall_our.push_back(recall);
        thresholds_our.push_back(threshold);
    }

    // Now, you have the precision, recall, and thresholds in the respective vectors
    // You can use them as needed
    // cout << "precision" << endl;
    // for (int i = 0; i < precision_our.size(); i++) {
    //     cout << precision_our[i] << " " << endl; 
    //  }
    //  cout << endl;
    double precision_all = calculatePrecisionAll(precision_our, recall_our);

    std::cout << "Precision All: " << precision_all << std::endl;

}


int main() {
    // Load your data and set up variables
    int n_trees = 50;
    int max_samples = 100;
    int max_depth = 10;
    int max_buckets = 3;
    double epsilon = 0.1;
    double sample_axis = 1.0;
    double threshold = 0.0;

    streambuf* orig_buf = cout.rdbuf();
    std::cout << "\n >> PIDForest Parameters\n";
    std::cout << "\n >> n_trees:" << n_trees << "\n";
    std::cout << "\n >> max_samples:" << max_samples << "\n";
    std::cout << "\n >> max_depth:" << max_depth << "\n";
    std::cout << "\n >> max_buckets:" << max_buckets << "\n";
    std::cout << "\n >> epsilon:" << epsilon << "\n";
    std::cout << "\n >> sample_axis:" << sample_axis << "\n";
    std::cout << "\n >> threshold:" << threshold << "\n";
    
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
    
    string dataset = "cpu_utilization_asg_misconfiguration";  // [nyc_taxi, ambient_temperature_system_failure, machine_temperature_system_failure, cpu_utilization_asg_misconfiguration] Set your dataset name here
    
    // Read the CSV file
    string filename = "../data/numenta/" + dataset + ".csv";
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return 1;
    }
    std::cout << "\n >> Using "<< dataset << " dataset, present at " << filename << "\n";

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
    // Shingle works correctly -- write tests to compare C++ and Python outputs of this function
    int shingleSize = 10;
    vector<vector<double>> X = TimeSeries::shingle(value, shingleSize);
    X = TimeSeries::transpose(X);

    // Print the shingled and transposed data
    std::cout << "\n >> There are " << X.size() << " shingles \n";
    
    /*
    int shingle_count = 0;
    for (const auto& row : X) {
        shingle_count++;
        std::cout << "\n >> Shingle #" << shingle_count <<"\n";
        std::cout << " >> "; 
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    */

    std::cout << "\n >> Running PIDForest algorithm on " << dataset << " dataset\n";

    vector<vector<double>> pts;  // Provide your data here
    pts = X;
    std::cout << "\n >> Running Forest::fit\n";
    mForest.fit(pts);
    std::cout << "\n >> Forest::fit DONE\n";
    
    // Make predictions
    vector<int> indices;
    unordered_map<int, vector<double>> outliers;
    unordered_map<int, vector<double>> scores;
    unordered_map<int, double> pst;
    vector<double> our_scores;
    double err = 0.1;
    double pct = 0.0;
    // cout.rdbuf(NULL);
    // cout.rdbuf(orig_buf);
    
    std::cout << "\n >> Running Forest::predict\n";
    
    tie(indices, outliers, scores, pst, our_scores) = mForest.predict(pts, err, pct);
    std::cout << "\n >> Forest::predict DONE\n";
    
    // Negate all values in the our_scores vector
    for (double& score : our_scores) {
        score = -score;
    }
    int t1 = X.size();
    std::vector<int> y_slice(label.begin(), label.begin() + t1);
    calculatePrecisionRecallThresholds(y_slice, our_scores);
    // Output the results
    /*
    std::cout << "\n >> Printing outputs Forest::predict\n";
    cout << "\n >> Indices: ";
    cout << "\n >> " << indices.size() << endl;

    for (int idx : indices) {
        cout << idx << "\n";
    }
    */
    /*
    
    // Print the outliers and their scores
    for (const auto& entry : outliers) {
        int idx = entry.first;
        const vector<double>& outlierValues = entry.second;
        cout << "\n >> Outlier at index " << idx << ": ";
        for (double val : outlierValues) {
            cout << val << " ";
        }
        cout << endl;
    }
    
    // Print the scores for each index
    for (const auto& entry : scores) {
        int idx = entry.first;
        const vector<double>& scoreValues = entry.second;
        cout << "\n >> Scores for index " << idx << ": ";
        for (double val : scoreValues) {
            cout << val << " ";
        }
        cout << endl;
    }
    
    // Print the PST values for each index
    for (const auto& entry : pst) {
        int idx = entry.first;
        double pstValue = entry.second;
        cout << "\n >> PST for index " << idx << ": " << pstValue << endl;
    }
    
    // Print the our_scores values
    cout << "\n >> Our Scores: ";
    for (double val : our_scores) {
        cout << val << " ";
    }
    cout << endl;
    */
    return 0;
}
