#include<iostream>
#include<vector>
#include<math.h>
#include<fstream>
#include<sstream>
#include<algorithm>
#include"custom_matrix.h"

// #define n 5
#define iter 100
#define dim 5


using namespace std;




//---------------------------------- ALGORITHM WRAPPERS ----------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------















  
int main(){

    cout<<"hello"<<endl;

    vector<vector<double>> source;
    vector<vector<double>> target;

    vector<int>source_labels;
    vector<int>target_labels;


    string filename = "S.csv";
    load_dataset(filename,source,source_labels);


    string filename2 ="Z.csv";
    load_dataset(filename2,target,target_labels);







}





