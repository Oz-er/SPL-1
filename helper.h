#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "custom_matrix.h"

#define iter 20
#define dim 5

using namespace std;



inline double reg = 0.1;



inline double dotproduct(vector<double> &v1,vector<double> &v2){
    double ans=0;

    int m = v1.size();

    for(int i=0;i<m;i++){
        ans+= v1[i]*v2[i];
    }

    return ans;
}




inline void load_dataset(string filename, vector<vector<double>> &X, vector<int> &Y) {

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
inline void z_score_normalize(vector<vector<double>> &data){
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




inline double eclnorm(vector<double> &v){
    
    double ans=0;

    for(int i=0;i<v.size();i++){
        ans += (v[i]*v[i]);
    }

    ans = sqrt(ans);

    return ans;

}


inline vector<double> normalize(vector<double> &v){
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


inline vector<double> subtract(vector<double> &first, vector<double> &second){

    int n= first.size();

    vector<double>result(n);

    for(int i=0;i<n;i++){
        result[i] = first[i] - second[i];
    }

    return result;
}


inline vector<double> proj (vector<double> &ground,vector<double> &v){
    double scalar = dotproduct(ground,v);

    int n = ground.size();
    vector<double>result(n);

    for(int i=0;i<n;i++){
        result[i] = ground[i]*scalar;
    }

    return result;
}


inline void qr_decmopose ( vector<vector<double>>&a ,vector<vector<double>>&q ,vector<vector<double>>&r ){
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


inline void printmat (vector<vector<double>> &v){

    for(int i=0;i<v.size();i++){
        for(int j=0;j<v[0].size();j++){
            cout<<v[i][j]<<" ";
        }

        cout<<endl;
    }
}





inline vector<pair<double,vector<double>>> geteigens (vector<vector<double>> &a,vector<vector<double>> &q,vector<vector<double>> &r){
   
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



inline auto eigen_starter(vector<vector<double>> &a){

    int s = a.size();

    vector<vector<double>> q(s,vector<double>(s));
    vector<vector<double>> r(s,vector<double>(s));

    return geteigens(a,q,r);

}

inline void printeigens (vector<double> &v){
    
    for(int i=0;i<v.size();i++){
        cout<<v[i]<<endl;
    }
    
}


inline vector<pair<double,vector<double>>> pair_sort(vector<pair<double,vector<double>>> &eigens){
    
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


//------------------------------------------------------- ingredients -----------------------------------------------------



inline vector<vector<double>> stacking (vector<vector<double>>&a, vector<vector<double>> &b){

   vector<vector<double>> mat;
   
    for(int i=0;i<a.size();i++){
    mat.push_back(a[i]);
    }

    for(int i=0;i<b.size();i++){
    mat.push_back(b[i]);
    }


    return mat;
}




inline vector<vector<double>> kernel(vector<vector<double>>&stk){

    int samples = stk.size();
    vector<vector<double>>k(samples,vector<double>(samples));

    for(int i=0;i<samples;i++){
        for(int j=i;j<samples;j++){
            double val= dotproduct(stk[i],stk[j]);
            
            k[i][j]=val;
            k[j][i]=val;
        
        }
    }


    return k;
}





inline vector<vector<double>> Centering(vector<vector<double>>&kernel){

    int samples = kernel.size();
    vector<vector<double>>c(samples,vector<double>(samples));

    for(int i=0;i<samples;i++){
        for(int j=0;j<samples;j++){
            

            if(i==j){
                c[i][j] = (1.00 - 1.00/samples);
            }

            else{
                c[i][j]= - 1.00/samples;
            }
        
        }
    }

    return c;
}









inline vector<vector<double>> weighting_mat(vector<vector<double>>&kernel ,vector<vector<double>>&source, vector<vector<double>>&target){

    int samples = kernel.size();
    vector<vector<double>>l(samples,vector<double>(samples));

    int ns = source.size();
    int nt = target.size();

    int total = ns + nt;


    for(int i=0;i<samples;i++){
        for(int j=0;j<samples;j++){

            if(i<ns && j<ns){
                l[i][j] = 1.00/(ns*ns);
            }

            else if(i>=ns && j>=ns){
                l[i][j] = 1.00/(nt*nt);
            }

            else{
                l[i][j] = -1.00/(ns*nt);
            }
        
        }
    }

    return l;
}






//---------------------------------- BDA SPECIFIC  ----------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------







inline vector<vector<double>> conditional_weighting_mat(int total_samples,int ns, int nt , 
                                                vector<int>&source_labels,
                                                vector<int>&target_pseudo_labels,
                                                int class_label){


    vector<vector<double>> Mc(total_samples,vector<double>(total_samples,0.0));
    

    int source_samples_inthis_class =0;
    int target_samples_inthis_class =0;


    for(int i=0;i<ns;i++){
        if(source_labels[i]==class_label){
            source_samples_inthis_class++;
        }
    }

    for(int i=0;i<nt;i++){
        if(target_pseudo_labels[i]==class_label){
            target_samples_inthis_class++;
        }
    }


    

    if(source_samples_inthis_class==0 || target_samples_inthis_class==0){
        return Mc;
    }




    for(int i=0;i<total_samples;i++){
        for(int j=0;j<total_samples;j++){
            
            bool i_is_in_class = false;
            if(i<ns){
                if(source_labels[i]==class_label){
                    i_is_in_class = true;
                }
            }
            else{
                if(target_pseudo_labels[i-ns]==class_label){
                    i_is_in_class=true;
                }
            }



            bool j_is_in_class = false;
            if(j<ns){
                if(source_labels[j]==class_label){
                    j_is_in_class = true;
                }
            }
            else{
                if(target_pseudo_labels[j-ns]==class_label){
                    j_is_in_class=true;
                }
            }







            if(i_is_in_class && j_is_in_class){
                if(i<ns && j<ns){
                    Mc[i][j] = 1.00/(source_samples_inthis_class*source_samples_inthis_class);
                }

                else if(i>=ns && j>=ns){
                    Mc[i][j] = 1.00/(target_samples_inthis_class*target_samples_inthis_class);
                }

                else{
                    Mc[i][j] = -1.00/(source_samples_inthis_class*target_samples_inthis_class);
                }   
            }

            


        }
    }


    return Mc;





    }






inline double calculate_mmd(vector<vector<double>>&data, int source_rows, int target_rows){
    

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




//----------------------------------  CORAL SPECIFIC ----------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------





inline vector<vector<double>> get_covariance(vector<vector<double>> &m){

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
    auto cov = matmult(m_centered_transpose,m_centered);//cols*cols
    for(int i=0;i<cols;i++){
    for(int j=0;j<cols;j++){
    cov[i][j]=cov[i][j]/(rows-1);
    }
    }

    return cov;


}








//mat mat_power = U*D*U^T
//U = eigenvctors stored in column
//D = power*identity matrix(diagonal matrix with the eigenvalue 
//raised to the power as the only value)
//U^T = Transpose of U
//only will work for square matrices
inline vector<vector<double>> mat_power  (vector<vector<double>> &m,  double power){

    int rows = m.size();
    int cols = m[0].size();


    auto eigs = eigen_starter(m);
    vector<vector<double>> U(rows,vector<double>(cols));
    vector<double> eig_values(rows);
    
    for(int  i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
    U[i][j]=eigs[j].second[i];
    }
    }


    for(int i=0;i<rows;i++){
    eig_values[i]=eigs[i].first;
    }



    auto U_T = transpose(U);



    vector<vector<double>> D(rows,vector<double>(cols,0.0));
    for(int i=0;i<rows;i++){
    for(int j=0;j<cols;j++){
    if(i==j){
    if(eig_values[i]>0.0){
    D[i][j]=pow((eig_values[i]),power);
    }
    }
    }
    }



    vector<vector<double>> res(rows,vector<double>(cols));
    res = matmult(U,D);
    res = matmult(res,U_T);

    return res;

}














//---------------------------------- training a classifier ----------------------------------




inline double get_distance(vector<double> &v1, vector<double> &v2){

    double sum = 0.0;
    int cols = v1.size();

    for(int i=0; i<cols;i++){
        double diff= v1[i]-v2[i];
        sum+=diff*diff;
    }

    return sqrt(sum);

}


inline bool compareDist(const pair<int,double>&a,const pair<int,double>&b){
    return a.second<b.second;
}


inline int knn_predict(vector<vector<double>> &train, vector<int>&train_labels, vector<double>&test_sample, int k){

    vector<pair<int,double>> distances;

    int train_rows = train.size();
    int cols = train[0].size();


//-------------calculating the distance between-----------------
//-------------this test row and every train row----------------
//------------then sort for the shortest ones-------------------
    for(int i=0;i<train_rows;i++){
    double d= get_distance(train[i],test_sample);
    distances.push_back({i,d});
    }
    sort(distances.begin(),distances.end(),compareDist);
//---------------------------------------------------------------
//---------------------------------------------------------------
//---------------------------------------------------------------

    int vote0 = 0;
    int vote1 = 0;


//------------------ count the labels of the --------------------
//----------------- nearest neighbors found above------------------

    for (int i = 0; i < k; i++) {
    int neighbor_index = distances[i].first;
    int label = train_labels[neighbor_index];
    if (label == 0) vote0++;
    else vote1++;
    }
//----------------------------------------------------------------
    // 4. Majority Rule
    if (vote1 > vote0) return 1;
    else return 0;

}



inline double get_knn_prob(vector<vector<double>> &train, vector<int>&train_labels, vector<double>&test_sample, int k){

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


inline double calculate_auc(vector<double>&probs,vector<int>&actual_labels){

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



#endif