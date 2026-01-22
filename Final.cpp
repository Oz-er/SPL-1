#include<iostream>
#include<vector>
#include<math.h>
#include<fstream>
#include<sstream>
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







int main(){

    
    vector<vector<double>> source;
    vector<vector<double>> target;

    vector<int>source_labels;
    vector<int>target_labels;


    string filename = "S.csv";
    load_dataset(filename,source,source_labels);


    string filename2 = "Z.csv";
    load_dataset(filename2,target,target_labels);

    z_score_normalize(source);
    z_score_normalize(target);

    auto stk=stacking(source,target); //20*10
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

    auto eigen_pairs = eigen_starter(W); 
    eigen_pairs = pair_sort(eigen_pairs);


    // for(int i=0;i<eigen_pairs.size();i++){
    //     cout<<"value: "<<eigen_pairs[i].first<<" "<<endl;

    //     for(int j=0;j<eigen_pairs[i].second.size();j++){
    //         cout<<"vector: "<<eigen_pairs[i].second[j]<<" ";
    //     }
    //     cout<<endl;
    // }


    vector<vector<double>>top_eigenvectors;

    for(int i=0;i<dim;i++){
        top_eigenvectors.push_back(eigen_pairs[i].second);
    }


    // int m = K.size(); // Total samples (Source + Target)
    // vector<vector<double>>W_matrix(m,vector<double>(dim));
    
    // for(int i=0;i<m;i++){
    //     for(int j=0;j<dim;j++){
    //         W_matrix[i][j]=eigen_pairs[j].second[i];
    //     }
    // }

    // auto Z = matmult(K,W_matrix);






    



}





