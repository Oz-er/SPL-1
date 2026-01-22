#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>

using namespace std;

void load_dataset(string filename, vector<vector<double>> &X, vector<int> &Y) {

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    string line;
    getline(file, line); // skip header

    
    while(getline(file,line)){
        stringstream ss(line);
        string val;
        vector<double>row;


        while(getline(ss,val,',')){
            row.push_back(stod(val));
        }


        if(row.empty())continue;

        int label = (int)row.back();
        Y.push_back(label);
        row.pop_back();

        X.push_back(row);
    }


     if (!X.empty()) {
        cout << "Loaded " << filename << ": "
             << X.size() << " samples with "
             << X[0].size() << " features." << endl;
    }


}



int main(){

    vector<vector<double>>mat;
    vector<int>labels;

    string filename = "S.csv";
    load_dataset(filename,mat,labels);

    for(int i=0;i<mat.size();i++){
        for(int j=0;j<mat[0].size();j++){
            cout<<mat[i][j]<<" ";
        }
        cout << "| Label: " << labels[i]; 

        cout<<endl;
    }  
    
}