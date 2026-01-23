#include<bits/stdc++.h>
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


}


double get_distance(vector<double> &v1, vector<double> &v2){

    double sum = 0.0;
    int cols = v1.size();

    for(int i=0; i<cols;i++){
        double diff= v1[i]-v2[i];
        sum+=diff*diff;
    }

    return sqrt(sum);

}


bool compareDist(const pair<int,double>&a,const pair<int,double>&b){
    return a.second<b.second;
}


int knn_predict(vector<vector<double>> &train, vector<int>&train_labels, vector<double>&test_sample, int k){

    vector<pair<int,double>> distances;

    int train_rows = train.size();
    int cols = train[0].size();

    for(int i=0;i<train_rows;i++){
        double d= get_distance(train[i],test_sample);
        distances.push_back({i,d});
    }


    sort(distances.begin(),distances.end(),compareDist);


    int vote0 = 0;
    int vote1 = 0;
for (int i = 0; i < k; i++) {
        int neighbor_index = distances[i].first;
        int label = train_labels[neighbor_index];

        if (label == 0) vote0++;
        else vote1++;
    }

    // 4. Majority Rule
    if (vote1 > vote0) return 1;
    else return 0;

}



int main(){



    vector<vector<double>> source;
    vector<vector<double>> target;

    vector<int>source_labels;
    vector<int>target_labels;


    string filename = "S.csv";
    load_dataset(filename,source,source_labels);


    string filename2 = "Z.csv";
    load_dataset(filename2,target,target_labels);


    vector<vector<double>>Z;

    vector<vector<double>>train;
    vector<vector<double>>test;



    for(int i=0;i<source.size();i++){
        train.push_back(Z[i]);
    }

    for(int i=source.size();i<source.size()+target.size();i++){
        test.push_back(Z[i]);
    }

    int k_neighbors=1;



    int true_pos =0;
    int true_neg =0;
    int false_pos =0;
    int false_neg =0;


    for(int i=0;i<test.size();i++){
        int actual = target_labels[i];

        int predicted = knn_predict(train,source_labels,test[i],k_neighbors);

        if(predicted==1 && actual ==1)true_pos++;
        if(predicted==0 && actual ==0)true_neg++;
        if(predicted==1 && actual ==0)false_pos++;
        if(predicted==0 && actual ==1)false_neg++;
    }


    // 3. Calculate Metrics
    double accuracy = (double)(true_pos + true_neg) / test.size() * 100.0;
    
    // Prevent division by zero
    double precision = (true_pos + false_pos) > 0 ? (double)true_pos / (true_pos + false_pos) : 0.0;
    double recall = (true_pos + false_neg) > 0 ? (double)true_pos / (true_pos + false_neg) : 0.0;
    double f1 = (precision + recall) > 0 ? 2.0 * (precision * recall) / (precision + recall) : 0.0;

    cout << "------------------------------------------------" << endl;
    cout << "FINAL RESULTS (TCA + KNN)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "True Positives (Caught Bugs): " << true_pos << endl;
    cout << "False Negatives (Missed Bugs): " << false_neg << endl;
    cout << "------------------------------------------------" << endl;
    cout << "Accuracy:  " << accuracy << "%" << endl;
    cout << "Precision: " << precision << endl;
    cout << "Recall:    " << recall << endl;
    cout << "F1-Score:  " << f1 << endl;
    cout << "------------------------------------------------" << endl;





}