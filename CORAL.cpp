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



double reg = 0.1;



double dotproduct(vector<double> &v1,vector<double> &v2){
    double ans=0;

    int m = v1.size();

    for(int i=0;i<m;i++){
        ans+= v1[i]*v2[i];
    }

    return ans;
}




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



//normalize each column independently
void z_score_normalize(vector<vector<double>> &data){
    if(data.empty())return;
    int rows = data.size();
    int cols = data[0].size();


    for(int i=0;i<cols;i++){
        double sum =0.0;
        double sq_sum = 0.0;

        for(int j=0;j<rows;j++){
            sum += data[j][i];
        }
        double mean = sum/rows;


        for(int j=0;j<rows;j++){
            sq_sum += pow(data[j][i]-mean,2);
        }
        double sd = sqrt(sq_sum/rows);
        
        if(sd<1e-9) sd=1.0;


        for (int j=0;j<rows;j++) {
            data[j][i] = (data[j][i] - mean) / sd;
        }
    }
}




double eclnorm(vector<double> &v){
    
    double ans=0;

    for(int i=0;i<v.size();i++){
        ans += (v[i]*v[i]);
    }

    ans = sqrt(ans);

    return ans;

}


vector<double> normalize(vector<double> &v){
    double norm =eclnorm(v);

    vector<double>v2(v.size());


    // const double epsilon = 1e-10;
    // if(norm < epsilon){
    //     for(int i=0; i<n; i++) v2[i] = 0;
    //     if(n > 0) v2[0] = 1.0;  // Make it a unit vector
    //     return v2;
    // }

    for(int i=0;i<v.size();i++){
        v2[i] = v[i]/norm;
    }

    return v2;
}


vector<double> subtract(vector<double> &first, vector<double> &second){

    int n= first.size();

    vector<double>result(n);

    for(int i=0;i<n;i++){
        result[i] = first[i] - second[i];
    }

    return result;
}


vector<double> proj (vector<double> &ground,vector<double> &v){
    double scalar = dotproduct(ground,v);

    int n = ground.size();
    vector<double>result(n);

    for(int i=0;i<n;i++){
        result[i] = ground[i]*scalar;
    }

    return result;
}


void qr_decmopose ( vector<vector<double>>&a ,vector<vector<double>>&q ,vector<vector<double>>&r ){
    q=a;


    for(int i=0;i<r.size();i++){
        for(int j=0;j<r[0].size();j++){
            r[i][j]=0;
        }
    }

    // //check for near 0 values
    // const double epsilon = 1e-10;

    //iterating through the column

    int n = a.size();

    for(int i=0;i<n;i++){ 

        //extracting the column
        vector<double>q_i(n);
        q_i = getcolumn(q,i);


        //diagonal of r will be the norms
        r[i][i] = eclnorm(q_i);

        // if(r[i][i] < epsilon){
        //     continue;
        // }

        //setting the unit vector
        q_i = normalize(q_i);
        setcolumn(q,q_i,i);

        for(int j=i+1;j<n;j++){

            
            vector<double>q_j(n);
            q_j = getcolumn(q,j);


            //setting the off diagonal of r
            double shadow = dotproduct(q_i,q_j);
            r[i][j]=shadow;


            //subtracting the projection vector 

            vector<double>projection(n);
            projection = proj(q_i,q_j);
            q_j = subtract(q_j,projection);


            setcolumn(q,q_j,j);
        }
    }
}


void printmat (vector<vector<double>> &v){

    for(int i=0;i<v.size();i++){
        for(int j=0;j<v[0].size();j++){
            cout<<v[i][j]<<" ";
        }

        cout<<endl;
    }
}





vector<pair<double,vector<double>>> geteigens (vector<vector<double>> &a,vector<vector<double>> &q,vector<vector<double>> &r){
   
    int n = a.size();
    vector<pair<double,vector<double>>>eigens(n);

    // vector<double> eigens(n);

    vector<vector<double>>identity(n,vector<double>(n));

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j){
                identity[i][j] =1.00;
            }
        }
    }
    

    for(int i=0;i<iter;i++){
        qr_decmopose(a,q,r);
        a=matmult(r,q);

        identity = matmult(identity,q);
        

    }

    for(int i=0;i<n;i++){
        eigens[i].first=a[i][i];
    }

    for(int i=0;i<n;i++){
        eigens[i].second = getcolumn(identity,i);
    }

    return eigens;
}



auto eigen_starter(vector<vector<double>> &a){

    int s = a.size();

    vector<vector<double>> q(s,vector<double>(s));
    vector<vector<double>> r(s,vector<double>(s));

    return geteigens(a,q,r);

}

void printeigens (vector<double> &v){
    
    for(int i=0;i<v.size();i++){
        cout<<v[i]<<endl;
    }
    
}


vector<pair<double,vector<double>>> pair_sort(vector<pair<double,vector<double>>> &eigens){
    
    bool swapped;

    int n=eigens.size();
    
    for (int i=0;i<n-1;i++) {
        swapped = false;
        for (int j=0;j<n-i-1;j++) {
            if (eigens[j].first < eigens[j + 1].first) {
                swap(eigens[j], eigens[j + 1]);
                swapped = true;
            }
        }
    
        if (!swapped)
            break;
    }
   return eigens;
}





//-------------------------calculating covariance matrix ---------------------



vector<vector<double>> get_covariance(vector<vector<double>> &m){

    int rows = m.size();
    int cols = m[0].size();




    vector<double>means(cols,0.0);
    for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
    means[j]+=m[i][j];
    }
    }
    for(int j=0;j<cols;j++){
    means[j]/=rows;
    }



    //subtracting the means of each feature 
    //column from every single data
    //to center the data
    vector<vector<double>> m_centered(rows,vector<double>(cols));
    for(int i=0; i<rows; i++){
    for(int j=0; j<cols; j++){
    m_centered[i][j]=m[i][j]-means[j];
    }
    }



    // Covariance = (X_c^T * X_c) / (n - 1)
    auto m_centered_transpose = transpose(m_centered);
    auto cov = matmult(m_centered_transpose,m_centered);
    for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
    cov[i][j]=cov[i][j]/(rows-1);
    }
    }

    return cov;


}













vector<vector<double>> mat_power  (vector<vector<double>> &m,  double power){

    int rows = m.size();
    int cols = m[0].size();


    auto eigs = eigen_starter(m);
    vector<vector<double>> eig_vectors;
    vector<double> eig_values;
    
    for(int  i=0;i<eigs.size();i++){
    eig_values.push_back(eigs[i].first);
    eig_vectors.push_back(eigs[i].second);
    }


    vector<vector<double>> idt;

    for(int  i=0;i<eigs.size();i++){
        
        
    }

    //note down the formula before proceeding 






}









double calculate_mmd(vector<vector<double>>&data, int source_rows, int target_rows){
    

    int dimensions = data[0].size();

    vector<double> mean_s(dimensions,0.0);
    vector<double> mean_t(dimensions,0.0);

    for(int i=0;i<source_rows;i++){
        for(int j=0;j<dimensions;j++){
            mean_s[j]+=data[i][j];
        }
    }

    for(int j=0;j<dimensions;j++){
        mean_s[j]/=source_rows;
    }


    for(int i=source_rows;i<source_rows+target_rows;i++){
        for(int j=0;j<dimensions;j++){
            mean_t[j]+=data[i][j];
        }
    }

    for(int j=0;j<dimensions;j++){
        mean_t[j]/=target_rows;
    }


    double dist_sq=0.0;

    for(int j=0; j<dimensions;j++){
        dist_sq += pow(mean_s[j]-mean_t[j],2);
    }

    return sqrt(dist_sq);


}




//---------------------------------- training a classifier ----------------------------------




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

double get_knn_prob(vector<vector<double>> &train, vector<int>&train_labels, vector<double>&test_sample, int k){

    vector<pair<int,double>>distances;

    int train_rows = train.size();

    for (int i = 0; i < train_rows; i++) {
        double d = get_distance(train[i],test_sample);
        distances.push_back({i,d});
    }

    sort(distances.begin(),distances.end(),compareDist);

    int vote1 = 0;

    for(int i=0;i<k;i++){
        int neighbor_index = distances[i].first;

        if(train_labels[neighbor_index]==1){
            vote1++;
        }
    }

    return (double)vote1/k;

}



//--------------------------------calculating arewa under the curve--------------------------------------




double calculate_auc(vector<double>&probs,vector<int>&actual_labels){

    vector<double>pos_scores;
    vector<double>neg_scores;

    for(int i = 0; i < probs.size(); i++) {
        if(actual_labels[i] == 1) pos_scores.push_back(probs[i]);
        else neg_scores.push_back(probs[i]);
    }
    
    if(pos_scores.empty() || neg_scores.empty()) return 0.0;

    double valid_pairs = 0.0;
    double total_pairs = 0.0;

    // Compare every Positive vs Negative sample
    for(double p_score : pos_scores) {
        for(double n_score : neg_scores) {
            total_pairs++;
            if(p_score > n_score) valid_pairs += 1.0;       // Correct ranking (Bug score > Clean score)
            else if(p_score == n_score) valid_pairs += 0.5; // Tie
        }
    }

    return valid_pairs / total_pairs;
}


  
int main(){

    cout<<"hello"<<endl;

    vector<vector<double>> source;
    vector<vector<double>> target;

    vector<int>source_labels;
    vector<int>target_labels;


    string filename = "datasets/A.csv";
    load_dataset(filename,source,source_labels);


    string filename2 ="datasets/Z.csv";
    load_dataset(filename2,target,target_labels);


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
    
   


//---------BDA SETUP-------------


    double mu = 0.5;
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

        // d. Project the data
        Z = matmult(K, W_matrix);

        // e. Split Z back into Train and Test
        train_data.clear();
        test_data.clear();
        for(int i=0; i<source.size(); i++) train_data.push_back(Z[i]);
        for(int i=source.size(); i<Z.size(); i++) test_data.push_back(Z[i]);

        // f. Update Pseudo-Labels for the NEXT iteration
        for(int i=0; i<test_data.size(); i++){
            pseudo_labels[i] = knn_predict(train_data, source_labels, test_data[i], 1);
        }
        
        cout << "Completed BDA Iteration: " << itr+1 << endl;
    }


    
    double mmd_after = calculate_mmd(Z, source.size(), target.size());
    cout << "\n------------------------------------------------" << endl;
    cout << "MMD Distance (After BDA): " << mmd_after << endl;
    cout << "------------------------------------------------" << endl;



    // --- FINAL EVALUATION ---
    int true_pos = 0, true_neg = 0, false_pos = 0, false_neg = 0;
    vector<double> prob_scores; 

    // Test using the FINAL latents space (Z) after 10 iterations
    for(int i=0; i<test_data.size(); i++){
        int actual = target_labels[i];
        
        // Using k=3 for final evaluation 
        int predicted = knn_predict(train_data, source_labels, test_data[i], 3);
        double prob = get_knn_prob(train_data, source_labels, test_data[i], 3);
        
        prob_scores.push_back(prob);

        if(predicted==1 && actual==1) true_pos++;
        if(predicted==0 && actual==0) true_neg++;
        if(predicted==1 && actual==0) false_pos++;
        if(predicted==0 && actual==1) false_neg++;
    }

    double accuracy = (double)(true_pos + true_neg) / test_data.size() * 100.0;
    double precision = (true_pos + false_pos) > 0 ? (double)true_pos / (true_pos + false_pos) : 0.0;
    double recall = (true_pos + false_neg) > 0 ? (double)true_pos / (true_pos + false_neg) : 0.0;
    double f1 = (precision + recall) > 0 ? 2.0 * (precision * recall) / (precision + recall) : 0.0;
    double auc = calculate_auc(prob_scores, target_labels);

    cout << "------------------------------------------------" << endl;
    cout << "FINAL RESULTS (BDA + KNN)" << endl;
    cout << "------------------------------------------------" << endl;
    cout << "Accuracy:  " << accuracy << "%" << endl;
    cout << "Precision: " << precision << endl;
    cout << "Recall:    " << recall << endl;
    cout << "F1-Score:  " << f1 << endl;
    cout << "AUC:       " << auc << endl;
    cout << "------------------------------------------------" << endl;

    return 0; // End of main()


}





