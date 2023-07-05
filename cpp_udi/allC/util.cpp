#include "util.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

void printVector(std::vector<double> & v, std::string s) {
    cout<<"length of "<<s<<": "<<v.size()<<endl;
    cout<<"Content of "<<s<<": ";
    for (int i = 0; i < v.size(); i++)
        cout<<v[i]<<", ";
    cout<<endl;
}

void printVector(std::vector<int> & v, std::string s) {
    cout<<"length of "<<s<<": "<<v.size()<<endl;
    cout<<"Content of "<<s<<": ";
    for (int i = 0; i < v.size(); i++)
        cout<<v[i]<<", ";
    cout<<endl;
}

void printVector(const std::vector<int> & v, std::string s) {
    cout<<"length of "<<s<<": "<<v.size()<<endl;
    cout<<"Content of "<<s<<": ";
    for (int i = 0; i < v.size(); i++)
        cout<<v[i]<<", ";
    cout<<endl;
}