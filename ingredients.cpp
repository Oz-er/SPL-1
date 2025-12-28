#include<iostream>
#include<vector>


using namespace std;

double dotproduct(vector<double> &v1,vector<double> &v2){
    double ans=0;

    int n = v1.size();

    for(int i=0;i<n;i++){
        ans+= v1[i]*v2[i];
    }

    return ans;
}




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
        for(int j=i;j<samples;j++){
            

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






vector<vector<double>> MMD(vector<vector<double>>&kernel ,vector<vector<double>>source, vector<vector<double>>target){

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





