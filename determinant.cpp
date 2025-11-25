#include<bits/stdc++.h>
using namespace std;
#define n 3



int det(vector<vector<int>> &v){

    int ans;

    ans = (v[0][0]*v[1][1])-(v[0][1]*v[1][0]);
    return ans;
    
}


int main(){

    vector<vector<int>>v(n,vector<int>(n));

    for(int i=0;i<n;i++){
        for(int j=0;j<3;j++){
            int a;
            cin>>a;

            v[i][j]=a;
        }
    }


    int row=0;
    int col=0;
    vector<vector<int>> v2(2,vector<int>(2));


    int k=0,l=0;

    int ans=0;
    bool odd=true;

    int num;

    while(col!=n){

        int num = v[row][col];

    

        
    for(int i=0;i<n;i++){     
        if(i==row)continue; 
        l=0;
        for(int j=0;j<n;j++){
            if(j==col)continue;
            v2[k][l]=v[i][j];
            l++;
        }
        k++;
    }



    if(odd){
        ans+= num*det(v2);
    }
    else{
        ans-= num*det(v2);
    }

    

    col++;
    k=0;
    l=0;
    odd=!odd;

    }


    cout<<ans<<endl;






}





    //     for(int i=0;i<2;i++){
    //     for(int j=0;j<2;j++){
    //         cout<<v2[i][j]<<" ";
    //     }                                 
    //     cout<<endl;

    // }
