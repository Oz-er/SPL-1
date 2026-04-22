#include <iostream>
#include <vector>
#include "helper.h"
#include "algorithms.h"

using namespace std;

void run_bda(vector<vector<double>>& source, vector<int>& source_labels, 
             vector<vector<double>>& target, vector<int>& target_labels) {
    
    cout << "\n--- Running BDA ---" << endl;
    double reg = 1; 


    auto stk=stacking(source,target); //20*10

    z_score_normalize(stk);

    double mmd_before = calculate_mmd(stk,source.size(),target.size());
    cout << "\n------------------------------------------------" << endl;
    cout << "MMD Distance (Original Data): " << mmd_before << endl;
    cout << "------------------------------------------------" << endl;


    auto K=kernel(stk); //(20*10).(10*20)=(20*20) 

    int m = K.size();

    auto H = Centering(K);
    auto tmp =matmult(H,K);
    auto kc = matmult(tmp,H); 

    // auto L = weighting_mat(K,source,target);
    // auto tmp2 = matmult(K,L);
    // auto kl = matmult(tmp2,K);
    //see line 800

    vector<vector<double>>muI(m,vector<double>(m));

    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
            if(i==j){
                muI[i][j] =reg;
            }
        }
    }
    
   


//------------------------------------------ BDA SETUP ---------------------------------------------------------


    double mu = 0.1;
    int iterations = 10;


    vector<int>pseudo_labels;

    for(int i=0;i<target.size();i++){
        int pred = knn_predict(source,source_labels,target[i],1);
        pseudo_labels.push_back(pred);
    }

    auto M0 = weighting_mat(K,source,target);

    vector<vector<double>> Z; 
    vector<vector<double>> train_data;
    vector<vector<double>> test_data;



//------------------------------------------BDA ITERATIONS--------------------------------------------------------




    for(int itr=0;itr<iterations;itr++){

        auto Mc_clean = conditional_weighting_mat(m,source.size(),target.size(),source_labels,pseudo_labels,0);
        auto Mc_buggy = conditional_weighting_mat(m,source.size(),target.size(),source_labels,pseudo_labels,1);



        auto Mc_combined = matadd(Mc_clean, Mc_buggy);
        auto part1 = mat_scalar_mult(M0, 1.0-mu);
        auto part2 = mat_scalar_mult(Mc_combined, mu);
        auto M_final = matadd(part1, part2);


        auto tmp2 = matmult(K,M_final);
        auto kl = matmult(tmp2,K);


        auto distance = matadd(kl, muI);
        auto structure = kc; // kc was calculated above the loop!
        auto W = matmult(mat_inverse(distance), structure);

        auto eigen_pairs = eigen_starter(W); 
        eigen_pairs = pair_sort(eigen_pairs);

        vector<vector<double>> W_matrix(m, vector<double>(dim));
        for(int i=0; i<m; i++){
            for(int j=0; j<dim; j++){
                W_matrix[i][j] = eigen_pairs[j].second[i];
            }
        }

        Z = matmult(K, W_matrix);

        train_data.clear();
        test_data.clear();
        for(int i=0; i<source.size(); i++) train_data.push_back(Z[i]);
        for(int i=source.size(); i<Z.size(); i++) test_data.push_back(Z[i]);



        int changes = 0;

        for(int i=0; i<test_data.size(); i++){
        int new_label = knn_predict(train_data, source_labels, test_data[i], 1);
        if(new_label != pseudo_labels[i]) changes++;
        pseudo_labels[i] = new_label;
        }


        cout << "Completed BDA Iteration: " << itr+1 << endl;    


        if(changes == 0) {
        cout << "Converged at iteration " << itr+1 << endl;
        break;
        }

    }


    
    double mmd_after = calculate_mmd(Z, source.size(), target.size());
    cout << "\n------------------------------------------------" << endl;
    cout << "MMD Distance (After BDA): " << mmd_after << endl;


//-----------------------------------EVALUATION--------------------------------------------------

    int k_neighbors= 1;
    
    int true_pos = 0, true_neg = 0, false_pos = 0, false_neg = 0;
    vector<double> prob_scores; 

     ofstream outfile("BDA_predictions.csv");
    outfile << "Actual_Label,Predicted_Label,Probability_Score\n"; 

    for(int i=0;i<test_data.size();i++){
        int actual = target_labels[i];

        int predicted = knn_predict(train_data,source_labels,test_data[i],k_neighbors);
        double prob = get_knn_prob(train_data, source_labels, test_data[i], k_neighbors);
        
        prob_scores.push_back(prob);

        outfile << actual << "," << predicted << "," << prob << "\n";

        if(predicted==1 && actual ==1)true_pos++;
        if(predicted==0 && actual ==0)true_neg++;
        if(predicted==1 && actual ==0)false_pos++;
        if(predicted==0 && actual ==1)false_neg++;
    }

    outfile.close();



    double accuracy = (double)(true_pos + true_neg) / test_data.size() * 100.0;

    double precision = (true_pos + false_pos) > 0 ? (double)true_pos / (true_pos + false_pos) : 0.0;
    
    double recall = (true_pos + false_neg) > 0 ? (double)true_pos / (true_pos + false_neg) : 0.0;
    
    double auc = calculate_auc(prob_scores, target_labels);

    double mmd_change_pct = 0.0;
    if (mmd_before != 0) {
        mmd_change_pct = ((mmd_before - mmd_after) / mmd_before) * 100.0;
    }

    cout << "FINAL RESULTS (BDA + KNN)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "MMD Change: " << mmd_change_pct << "%" << endl; 
    cout << "AUC:        " << auc << endl;
    cout << "Precision:  " << precision << endl;
    cout << "Recall:     " << recall << endl;
    cout << ">>> Predicted labels successfully saved to: BDA_predictions.csv" << endl;
    cout << "------------------------------------------------" << endl;

}