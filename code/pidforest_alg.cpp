#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <filesystem>
#include <chrono>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <Eigen/Dense>
#include <random>

#include "Forest.h"
#include "Node.h"
#include "PointSet.h"
#include "Cube.h"
#include "Histogram.h"

using namespace std;
using namespace Eigen;
namespace fs = std::filesystem;


#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>


vector < vector<string> > read_data(string fname)
{
	vector < vector<string> > content;
	vector <string> row;
	string line, word;
	fstream file (fname, ios::in);

	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
 
			stringstream str(line);
 
			while(getline(str, word, ','))
				row.push_back(word);
			content.push_back(row);
		}
	}

	// visualize
	/*
	for(int i=0; i < content.size(); i++) 
	{
		for(int j=0; j < content[i].size(); j++)
		{
			std::cout << ' ' << content[i][j] << ' ';
		}
		cout<<"\n";
	}
	*/
	return content;
}

vector <double> get_field(string name, vector < vector<string> > content)
{
	vector <double> field;
	map <string, int> map;
	map["timestamp"] = 0;
	map["value"] = 1;
	map["anomaly_score"] = 2;
	map["label"] = 3;
	
	for(int i=1; i < content.size(); i++) 
	{
		field.push_back(stod(content[i][map[name]]));
	}
	return field;
}

vector<vector<double>> shingle(vector<double> series, int dim) 
{
    int height = series.size() - dim + 1;
    vector<vector<double> > shingled(dim, vector<double>(height, 0));
    for (int i = 0; i < dim; i++) 
    {
        for (int j = 0; j < height; j++) 
        {
            shingled[i][j] = series[i+j];
        }
    }

    // visualize
	/*
    for (auto& row : shingled) {
        for (auto& val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
    */
    return shingled;
}

vector<vector<double>> transpose(vector<vector<double>> X)
{
    int rows = X.size();
    int cols = X[0].size();
    vector<vector<double>> X_transpose(cols, vector<double>(rows, 0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            X_transpose[j][i] = X[i][j];
        }
    }

    // visualize
    /*
    for (auto& row : X_transpose) 
    {
        for (auto& val : row) 
        {
            cout << val << " ";
        }
        cout << endl;
    }
    */
    return X_transpose;
}

vector<int> shape(vector<vector<double>> X) 
{
    int rows = X.size();
    int cols = X[0].size();
    vector<int> shape = {rows, cols};
    return shape;
}

double get_time()
{
	auto current_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = current_time.time_since_epoch();
    double timestamp = elapsed_seconds.count();
    return timestamp;
}


int main(int argc, char const *argv[])
{
	string datasets[] = {"nyc_taxi","ambient_temperature_system_failure","cpu_utilization_asg_misconfiguration","machine_temperature_system_failure"};
	int L = sizeof(datasets)/sizeof(datasets[0]);
	int trials = 1;
	int run_lof_svm = 1;
	for (int i = 0; i < L; ++i)
	{
		string fname = "../data/numenta/" + datasets[0] + ".csv";
		vector < vector<string> > content = read_data(fname);
		vector <double> arr = get_field("value", content);
		vector <double> y = get_field("anomaly_score", content);
		vector<vector<double> > X = shingle(arr, 10);
		vector<vector<double>> X_transpose = transpose(X);
		vector<int> shape_found = shape(X);
		int t1 = shape_found[0];
		std::string dir_name = "./experiment_results";
	    if (!fs::exists(dir_name)) {
	        fs::create_directory(dir_name);
	    }
	    string file_name = "./experiment_results/" + datasets[i] + ".txt";
	    std::ofstream file_object(file_name);
	    vector<vector<double> > time_all(trials, vector<double>(4, 0));
	    vector<vector<double> > precision_all(trials, vector<double>(4, 0));
	    vector<vector<double> > auc_all(trials, vector<double>(4, 0));

	    for (int j = 0; j < trials; j++)
	    {
	    	cout << "\n\n******" + datasets[i] + " trial "+to_string(j+1)+"*******\n\n";
	    	cout << "\n******Our Algo*******\n";
	    	double start = get_time();
	    	std::map<std::string, int> kwargs = {{"max_depth", 10}, {"n_trees", 50}, {"max_samples", n_samples}, 
                                     {"max_buckets", 3}, {"epsilon", 0.1}, {"sample_axis", 1}, {"threshold", 0}};
		}
	}
	
}