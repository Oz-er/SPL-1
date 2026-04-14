#include <iostream>
#include <vector>
#include "helper.h"
#include "algorithms.h"

using namespace std;

void run_coral(vector<vector<double>>& source, vector<int>& source_labels, 
               vector<vector<double>>& target, vector<int>& target_labels) {

    cout<<"\n--- Running CORAL ---" <<endl;

    vector<vector<double>> source_local = source;
    vector<vector<double>> target_local = target;

    z_score_normalize(source_local);
    z_score_normalize(target_local);



    auto covariance_source = get_covariance(source_local);
    auto covariance_target = get_covariance(target_local);
    //Regularize (Add Identity matrix so no eigenvalue is absolutely zero)
    int d = source[0].size();
    vector<vector<double>> I(d, vector<double>(d, 0.0));
    for(int i=0; i<d; i++) I[i][i] = 1.0; 
    covariance_source = matadd(covariance_source, I);
    covariance_target = matadd(covariance_target, I);





    //data whitening 
    auto cs_inverse_half = mat_power(covariance_source,-0.5);
    //data recoloring
    auto ct_half = mat_power(covariance_target, 0.5);
    //A = Cs^(-1/2) * Ct^(1/2)
    auto A = matmult(cs_inverse_half,ct_half);



    auto source_aligned = matmult(source_local,A);




    // --- FINAL EVALUATION ---
    int true_pos = 0, true_neg = 0, false_pos = 0, false_neg = 0;
    vector<double> prob_scores; 

    // Test the Target data against the newly Aligned Source data
    for(int i=0; i<target_local.size(); i++){
        int actual = target_labels[i];
        
        // We use K=3 for standard evaluation
        int predicted = knn_predict(source_aligned, source_labels, target_local[i], 3);
        double prob = get_knn_prob(source_aligned, source_labels, target_local[i], 3);
        
        prob_scores.push_back(prob);

        if(predicted==1 && actual==1) true_pos++;
        if(predicted==0 && actual==0) true_neg++;
        if(predicted==1 && actual==0) false_pos++;
        if(predicted==0 && actual==1) false_neg++;
    }

    double accuracy = (double)(true_pos + true_neg) / target_local.size() * 100.0;
    double precision = (true_pos + false_pos) > 0 ? (double)true_pos / (true_pos + false_pos) : 0.0;
    double recall = (true_pos + false_neg) > 0 ? (double)true_pos / (true_pos + false_neg) : 0.0;
    double f1 = (precision + recall) > 0 ? 2.0 * (precision * recall) / (precision + recall) : 0.0;
    double auc = calculate_auc(prob_scores, target_labels);

    cout << "------------------------------------------------" << endl;
    cout << "FINAL RESULTS (CORAL + KNN)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "Accuracy:  " << accuracy << "%" << endl;
    cout << "Precision: " << precision << endl;
    cout << "Recall:    " << recall << endl;
    cout << "F1-Score:  " << f1 << endl;
    cout << "AUC:       " << auc << endl;
    cout << "------------------------------------------------" << endl;


            










    }