#include <iostream>
#include <vector>
#include "helper.h"
#include "algorithms.h"

using namespace std;


void run_baseline(vector<vector<double>>& source, vector<int>& source_labels, 
                  vector<vector<double>>& target, vector<int>& target_labels) {
    
    cout << "\n--- Running Baseline (K-NN Only) ---" << endl;

    // 1. Stack and Normalize (Exactly like TCA for a fair comparison)
    auto stk = stacking(source, target); 
    z_score_normalize(stk);

    double mmd_before = calculate_mmd(stk, source.size(), target.size());
    cout << "\n------------------------------------------------" << endl;
    cout << "MMD Distance (Raw Data): " << mmd_before << endl;
    cout << "------------------------------------------------" << endl;

    // ==========================================
    // NOTICE: All TCA Matrix Math is DELETED here
    // ==========================================

    vector<vector<double>> train;
    vector<vector<double>> test;

    // 2. Split the normalized data directly (Using 'stk' instead of 'Z')
    for(int i = 0; i < source.size(); i++){
        train.push_back(stk[i]);
    }

    for(int i = source.size(); i < source.size() + target.size(); i++){
        test.push_back(stk[i]);
    }

    int k_neighbors = 1;

    int true_pos = 0;
    int true_neg = 0;
    int false_pos = 0;
    int false_neg = 0;

    vector<double> prob_scores; // Store probabilities for AUC  

    // Save predictions to a different file so it doesn't overwrite TCA
    ofstream outfile("baseline_predictions.csv");
    outfile << "Actual_Label,Predicted_Label,Probability_Score\n"; 

    for(int i = 0; i < test.size(); i++){
        int actual = target_labels[i];

        int predicted = knn_predict(train, source_labels, test[i], k_neighbors);
        double prob = get_knn_prob(train, source_labels, test[i], k_neighbors);
        
        prob_scores.push_back(prob);

        outfile << actual << "," << predicted << "," << prob << "\n";

        if(predicted == 1 && actual == 1) true_pos++;
        if(predicted == 0 && actual == 0) true_neg++;
        if(predicted == 1 && actual == 0) false_pos++;
        if(predicted == 0 && actual == 1) false_neg++;
    }

    outfile.close();

    double accuracy = (double)(true_pos + true_neg) / test.size() * 100.0;
    double precision = (true_pos + false_pos) > 0 ? (double)true_pos / (true_pos + false_pos) : 0.0;
    double recall = (true_pos + false_neg) > 0 ? (double)true_pos / (true_pos + false_neg) : 0.0;
    double f1 = (precision + recall) > 0 ? 2.0 * (precision * recall) / (precision + recall) : 0.0;
    double auc = calculate_auc(prob_scores, target_labels);

    cout << "FINAL RESULTS (BASELINE K-NN ONLY)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "Accuracy:  " << accuracy << "%" << endl;
    cout << "Precision: " << precision << endl;
    cout << "Recall:    " << recall << endl;
    cout << "F1-Score:  " << f1 << endl;
    cout << "AUC:       " << auc << endl; 
    cout << "------------------------------------------------" << endl;
    cout << ">>> Predicted labels successfully saved to: baseline_predictions.csv" << endl;
    cout << "------------------------------------------------" << endl;
}