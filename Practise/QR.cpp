#include<iostream>
#include<vector>
#include<math.h>
#include "custom_matrix.h"

using namespace std;

#define n 150
#define iter 1000


double dotproduct(vector<double> &v1,vector<double> &v2){
    double ans=0;

    for(int i=0;i<n;i++){
        ans+= v1[i]*v2[i];
    }

    return ans;
}



double eclnorm(vector<double> &v){
    
    double ans=0;

    for(int i=0;i<n;i++){
        ans += (v[i]*v[i]);
    }

    ans = sqrt(ans);

    return ans;

}


vector<double> normalize(vector<double> &v){
    double norm =eclnorm(v);

    vector<double>v2(n);


    // const double epsilon = 1e-10;
    // if(norm < epsilon){
    //     for(int i=0; i<n; i++) v2[i] = 0;
    //     if(n > 0) v2[0] = 1.0;  // Make it a unit vector
    //     return v2;
    // }

    for(int i=0;i<n;i++){
        v2[i] = v[i]/norm;
    }

    return v2;
}


vector<double> subtract(vector<double> &first, vector<double> &second){

    vector<double>result(n);

    for(int i=0;i<n;i++){
        result[i] = first[i] - second[i];
    }

    return result;
}


vector<double> proj (vector<double> &ground,vector<double> &v){
    double scalar = dotproduct(ground,v);

    vector<double>result(n);

    for(int i=0;i<n;i++){
        result[i] = ground[i]*scalar;
    }

    return result;
}


void qr_decmopose ( vector<vector<double>>&a ,vector<vector<double>>&q ,vector<vector<double>>&r ){
    q=a;


    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            r[i][j]=0;
        }
    }

    // //check for near 0 values
    // const double epsilon = 1e-10;

    //iterating through the column

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
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<v[i][j]<<" ";
        }

        cout<<endl;
    }
}


vector<double> geteigenvalues (vector<vector<double>> &a,vector<vector<double>> &q,vector<vector<double>> &r){
   
    vector<double> eigens(n);

    for(int i=0;i<iter;i++){
        qr_decmopose(a,q,r);
        a=matmult(r,q);
    }

    for(int i=0;i<n;i++){
        eigens[i]=a[i][i];
    }

    return eigens;
}


void printeigens (vector<double> &v){
    
    for(int i=0;i<n;i++){
        cout<<v[i]<<endl;
    }
    
}



int main(){

    cout<<"enter the matrix :"<<endl;


    vector<vector<double>> a(n,vector<double>(n));

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cin>>a[i][j];
        }
    }


    vector<vector<double>> q(n,vector<double>(n));
    vector<vector<double>> r(n,vector<double>(n));
   


    vector<double> eigens(n);

    eigens = geteigenvalues(a,q,r);
    printeigens(eigens);


    

}




