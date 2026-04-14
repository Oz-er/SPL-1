#include <iostream>
#include <vector>
#include "helper.h"
#include "algorithms.h"

using namespace std;

void run_tca(vector<vector<double>>& source, vector<int>& source_labels, 
             vector<vector<double>>& target, vector<int>& target_labels) {
    
    cout << "\n--- Running TCA ---" << endl;
    double reg = 0.1; 



    auto stk=stacking(source,target); //20*10

    z_score_normalize(stk);

    cout<<"before mmd"<<endl;

    double mmd_before = calculate_mmd(stk,source.size(),target.size());
    cout << "\n------------------------------------------------" << endl;
    cout << "MMD Distance (Original Data): " << mmd_before << endl;
    cout << "------------------------------------------------" << endl;

    cout<<"after mmd"<<endl;

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
    
   

    // cout<<"Stacked matrix: "<<endl;
    // for(int i=0;i<stk.size();i++){
    //     for(int j=0;j<stk[0].size();j++){
    //         cout<<stk[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    // cout<<endl;



    // cout<<"Kernel matrix: "<<endl;
    // for(int i=0;i<m;i++){
    //     for(int j=0;j<m;j++){
    //         cout<<K[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    // cout<<endl;



    // cout<<"Centering matrix: "<<endl;
    // for(int i=0;i<m;i++){
    //     for(int j=0;j<m;j++){
    //         cout<<kc[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    // cout<<endl;



    // cout<<"Weight matrix: "<<endl;
    // for(int i=0;i<m;i++){
    //     for(int j=0;j<m;j++){
    //         cout<<kl[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }

    // cout<<endl;


    auto distance = matadd(kl,muI);
    auto structure = kc;


    auto W = matmult(mat_inverse(distance),structure);

    // cout<<"1 done"<<endl;

    auto eigen_pairs = eigen_starter(W); 
    eigen_pairs = pair_sort(eigen_pairs);


    // cout<<"2 done"<<endl;


    // for(int i=0;i<eigen_pairs.size();i++){
    //     cout<<"value: "<<eigen_pairs[i].first<<" "<<endl;

    //     for(int j=0;j<eigen_pairs[i].second.size();j++){
    //         cout<<"vector: "<<eigen_pairs[i].second[j]<<" ";
    //     }
    //     cout<<endl;
    // }


    // vector<vector<double>>top_eigenvectors;

    // for(int i=0;i<dim;i++){
    //     top_eigenvectors.push_back(eigen_pairs[i].second);
    // }


    m = K.size(); // Total samples (Source + Target)
    vector<vector<double>>W_matrix(m,vector<double>(dim));
    
    for(int i=0;i<m;i++){
        for(int j=0;j<dim;j++){
            W_matrix[i][j]=eigen_pairs[j].second[i];
        }
    }
    cout<<"3 done"<<endl;

    auto Z = matmult(K,W_matrix);

    // cout<<"4 done"<<endl;



    double mmd_after = calculate_mmd(Z,source.size(),target.size());
    cout << "\n------------------------------------------------" << endl;
    cout << "MMD Distance (Projected Data): " << mmd_after << endl;
    cout << "------------------------------------------------" << endl;  



    
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


    vector<double> prob_scores; // Store probabilities for AUC  

    for(int i=0;i<test.size();i++){
        int actual = target_labels[i];

        int predicted = knn_predict(train,source_labels,test[i],k_neighbors);
    // 2. Get Probability Score (for AUC) 
        double prob = get_knn_prob(train, source_labels, test[i], k_neighbors);
        prob_scores.push_back(prob);
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

    double auc = calculate_auc(prob_scores, target_labels);

    cout << "------------------------------------------------" << endl;
    cout << "FINAL RESULTS (TCA + KNN)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "Accuracy:  " << accuracy << "%" << endl;
    cout << "Precision: " << precision << endl;
    cout << "Recall:    " << recall << endl;
    cout << "F1-Score:  " << f1 << endl;
    cout << "AUC:       " << auc << endl; // Print AUC
    cout << "------------------------------------------------" << endl;
    // cout<<"5 done"<<endl;



}