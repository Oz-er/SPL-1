#include <iostream>
#include <vector>
#include <fstream>
#include "helper.h"
#include "algorithms.h"

using namespace std;

void run_tca(vector<vector<double>>& source, vector<int>& source_labels, 
             vector<vector<double>>& target, vector<int>& target_labels) {
    
    cout << "\n--- Running TCA ---" << endl;
    double reg = 0.1; 


    auto stk=stacking(source,target); //20*10
    z_score_normalize(stk);


    double mmd_before = calculate_mmd(stk,source.size(),target.size());
    cout << "\n------------------------------------------------" << endl;
    cout << "MMD Distance (Original Data): " << mmd_before << endl;



//----------------------------------------------------------------------------
//------------ solves the generalized eigenvalue eqn--------------------------
//----------------------------------------------------------------------------

    auto K=kernel(stk); //(20*10).(10*20)=(20*20) 
    int m = K.size();
    auto H = Centering(K);
    auto tmp =matmult(H,K);
    auto kc = matmult(tmp,H); 
    auto L = weighting_mat(K,source,target);
    auto tmp2 = matmult(K,L);
    auto kl = matmult(tmp2,K);

    vector<vector<double>>muI(m,vector<double>(m));
    for(int i=0;i<m;i++){
    for(int j=0;j<m;j++){
        if(i==j){
        muI[i][j] =reg;
        }
    }
    }
    
    auto distance = matadd(kl,muI);
    auto structure = kc;
    auto W = matmult(mat_inverse(distance),structure);

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------






//----------------------------------------------------------------------------
//----------- capturing the top K eigenvectors -------------------------------
//----------------------------------------------------------------------------

    auto eigen_pairs = eigen_starter(W); 
    eigen_pairs = pair_sort(eigen_pairs);
    m = K.size(); 
    vector<vector<double>>W_matrix(m,vector<double>(dim));
    
    for(int i=0;i<m;i++){
        for(int j=0;j<dim;j++){
            W_matrix[i][j]=eigen_pairs[j].second[i];
        }
    }
    auto Z = matmult(K,W_matrix);  // Z = contains only 5 best features

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------







    double mmd_after = calculate_mmd(Z,source.size(),target.size());
    cout << "MMD Distance (Projected Data): " << mmd_after << endl;

    
    vector<vector<double>>train;
    vector<vector<double>>test;



//---------------------------------------------------------------------------
//---------------------- split the projected data ---------------------------
//---------------------- into train and test --------------------------------

    for(int i=0;i<source.size();i++){
    train.push_back(Z[i]);
    }

    for(int i=source.size();i<source.size()+target.size();i++){
    test.push_back(Z[i]);
    }

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



    int k_neighbors=1;

    int true_pos =0;
    int true_neg =0;
    int false_pos =0;
    int false_neg =0;

    vector<double> prob_scores;  

    ofstream outfile("TCA_predictions.csv");
    outfile << "Actual_Label,Predicted_Label,Probability_Score\n"; 

    for(int i=0;i<test.size();i++){
        int actual = target_labels[i];

        int predicted = knn_predict(train,source_labels,test[i],k_neighbors);
        double prob = get_knn_prob(train, source_labels, test[i], k_neighbors); //for AUC only
        
        prob_scores.push_back(prob);

        outfile << actual << "," << predicted << "," << prob << "\n";

        if(predicted==1 && actual ==1)true_pos++;
        if(predicted==0 && actual ==0)true_neg++;
        if(predicted==1 && actual ==0)false_pos++;
        if(predicted==0 && actual ==1)false_neg++;
    }




    outfile.close();

    double accuracy = (double)(true_pos + true_neg) / test.size() * 100.0;
    double precision = (true_pos + false_pos) > 0 ? (double)true_pos / (true_pos + false_pos) : 0.0;
    double recall = (true_pos + false_neg) > 0 ? (double)true_pos / (true_pos + false_neg) : 0.0;
    double f1 = (precision + recall) > 0 ? 2.0 * (precision * recall) / (precision + recall) : 0.0;
    double auc = calculate_auc(prob_scores, target_labels);

    double mmd_change_pct = 0.0;
    if (mmd_before != 0) {
        mmd_change_pct = ((mmd_before - mmd_after) / mmd_before) * 100.0;
    }

    cout << "FINAL RESULTS (TCA + KNN)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "MMD Change: " << mmd_change_pct << "%" << endl; 
    cout << "AUC:        " << auc << endl;
    cout << "Precision:  " << precision << endl;
    cout << "Recall:     " << recall << endl;
    cout << ">>> Predicted labels successfully saved to: TCA_predictions.csv" << endl;
    cout << "------------------------------------------------" << endl;

}