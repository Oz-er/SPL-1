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


//------------------------------------------------------- ingredients -----------------------------------------------------



vector<vector<double>> stacking (vector<vector<double>>&a, vector<vector<double>> &b){

   vector<vector<double>> mat;
   
    for(int i=0;i<a.size();i++){
    mat.push_back(a[i]);
    }

    for(int i=0;i<b.size();i++){
    mat.push_back(b[i]);
    }


    return mat;
}




vector<vector<double>> kernel(vector<vector<double>>&stk){

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





vector<vector<double>> Centering(vector<vector<double>>&kernel){

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









vector<vector<double>> weighting_mat(vector<vector<double>>&kernel ,vector<vector<double>>&source, vector<vector<double>>&target){

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








int main(){

    cout<<"hello"<<endl;

    vector<vector<double>> source;
    vector<vector<double>> target;

    vector<int>source_labels;
    vector<int>target_labels;


    string filename = "A.csv";
    load_dataset(filename,source,source_labels);


    string filename2 = "S.csv";
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

    cout<<"1 done"<<endl;

    auto eigen_pairs = eigen_starter(W); 
    eigen_pairs = pair_sort(eigen_pairs);


    cout<<"2 done"<<endl;


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



    // double mmd_after = calculate_mmd(Z,source.size(),target.size());
    // cout << "\n------------------------------------------------" << endl;
    // cout << "MMD Distance (Projected Data): " << mmd_after << endl;
    // cout << "------------------------------------------------" << endl;  



    
    vector<vector<double>>train;
    vector<vector<double>>test;



    for(int i=0;i<source.size();i++){
        train.push_back(Z[i]);
    }

    for(int i=source.size();i<source.size()+target.size();i++){
        test.push_back(Z[i]);
    }

    int k_neighbors=3;



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


    cout<<"5 done"<<endl;








}





