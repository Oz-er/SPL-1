#include<bits/stdc++.h>
using namespace std;
#define N 5



void getCofactor(vector<vector<int>>&mat,vector<vector<int>>&temp,int p,int q, int size){
   
   
   int row=p;
   int col=q;
   int k=0;
   int l;

   
    for(int i=0;i<size;i++){     
        if(i==row)continue; 
        l=0;
        for(int j=0;j<size;j++){
            if(j==col)continue;
            temp[k][l]=mat[i][j];
            l++;
        }
        k++;
    }


}

int determinantOfMatrix(vector<vector<int>> &mat,int size){

    int d=0;

    if(size==1){
        return mat[0][0];
    }

    vector<vector<int>>temp(size,vector<int>(size));
    int sign =1;

    for(int col=0;col<size;col++){
        getCofactor(mat,temp,0,col,size);

        d+= sign*mat[0][col]*determinantOfMatrix(temp,size-1);

        sign=-sign;
    }
    return d;  

}



int main(){

    int n;
    cin>>n;

    vector<vector<int>>v(n,vector<int>(n));

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cin>>v[i][j];
        }
    }

    // Function call
    cout << "Determinant of the matrix is : " << 
             determinantOfMatrix(v, n);
    return 0;
}

