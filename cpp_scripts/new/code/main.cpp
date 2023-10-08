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
        double current_precision_all = 2.0 * precision * recall / (precision + recall);
        if (current_precision_all > max_precision_all) {
            max_precision_all = current_precision_all;
        }
    }

    return max_precision_all;
}

void calculatePrecisionRecallThresholds(const std::vector<int>& y, const std::vector<double>& our_scores) {
    int t1 = y.size();

    std::vector<double> precision_our;
    std::vector<double> recall_our;
    std::vector<double> thresholds_our;

    for (double threshold = 0.0; threshold <= 1.0; threshold += 0.01) {
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

        double precision = static_cast<double>(true_positives) / (true_positives + false_positives);
        double recall = static_cast<double>(true_positives) / (true_positives + false_negatives);

        precision_our.push_back(precision);
        recall_our.push_back(recall);
        thresholds_our.push_back(threshold);
        // std::cout << "True +ve: " << true_positives << ", False +ve: " << false_positives << ", False -ve: " << false_negatives << std::endl;
        // std::cout << "Threshold: " << threshold << ", Precision: " << precision << ", Recall: " << recall << std::endl;
    }

    double precision_all = calculatePrecisionAll(precision_our, recall_our);

    // std::cout << "Precision All: " << precision_all << std::endl;

}

void calculate_f1_score(const std::vector<int>& y, const std::vector<double>& our_scores, double threshold) {
    int t1 = y.size();
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

    double precision = static_cast<double>(true_positives) / (true_positives + false_positives);
    double recall = static_cast<double>(true_positives) / (true_positives + false_negatives);

    // std::cout << "True +ve: " << true_positives << ", False +ve: " << false_positives << ", False -ve: " << false_negatives << std::endl;
    // std::cout << "Threshold: " << threshold << ", Precision: " << precision << ", Recall: " << recall << std::endl;
    double f1_score = 2.0 * precision * recall / (precision + recall);
    // std::cout << "F1_Score: " << f1_score << std::endl;
}

int main() {
    std::cout.precision(10);
    int n_trees = 50; //50;
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

    getline(file, line);

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

    int shingleSize = 10;
    vector<vector<double>> X = TimeSeries::shingle(value, shingleSize);
    X = TimeSeries::transpose(X);

    std::cout << "\n >> There are " << X.size() << " shingles \n";

    std::cout << "\n >> Running PIDForest algorithm on " << dataset << " dataset\n";

    vector<vector<double>> pts; 
    pts = X;
    std::cout << "\n >> Running Forest::fit\n";
    mForest.fit(pts);
    std::cout << "\n >> Forest::fit DONE\n";
    
    vector<int> indices;
    unordered_map<int, vector<double>> outliers;
    unordered_map<int, vector<double>> scores;
    unordered_map<int, double> pst;
    vector<double> our_scores;
    double err = 0.1;
    double pct = 0.0;
    
    std::cout << "\n >> Running Forest::predict\n";
    
    tie(indices, outliers, scores, pst, our_scores) = mForest.predict(pts, err, pct);
    std::cout << "\n >> Forest::predict DONE\n";
    
    for (double& score : our_scores) {
        if (score < 0){
            score = -score;
        }
    }
    
    int t1 = X.size();
    std::vector<int> y_slice(label.begin(), label.begin() + t1);
    // for (int t=0; t < y_slice.size(); t++) {
    //     std::cout << y_slice[t] << std::endl;
    // }

    // Create a file stream
    std::ofstream outputFile("y_slice.txt");

    // Check if the file stream is open
    if (!outputFile.is_open()) {
        std::cerr << "Failed to open the file for writing." << std::endl;
        return 1; // Return an error code
    }

    // Store the elements from the beginning up to t1 in the file
    for (int i = 0; i < y_slice.size(); ++i) {
        outputFile << y_slice[i] << "\n"; // Write each element to a new line
    }

    // Close the file
    outputFile.close();

    // Create a file stream
    std::ofstream outputFile2("our_scores.txt");

    // Check if the file stream is open
    if (!outputFile2.is_open()) {
        std::cerr << "Failed to open the file for writing." << std::endl;
        return 1; // Return an error code
    }

    // Store the elements from the beginning up to t1 in the file
    for (int i = 0; i < our_scores.size(); ++i) {
        outputFile2 << our_scores[i] << "\n"; // Write each element to a new line
    }

    // Close the file
    outputFile2.close();
    calculatePrecisionRecallThresholds(y_slice, our_scores);
    calculate_f1_score(y_slice, our_scores, threshold);
    return 0;
}
