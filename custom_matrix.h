#ifndef CUSTOM_MATRIX_H
#define CUSTOM_MATRIX_H

#include<vector>
#include<utility>


std::pair<int,int> get_rowcol(std::vector<std::vector<double>> &v){

    if(v.empty()){
        return {0,0};
    }

    return {v.size(),v[0].size()};
}



std::vector<double> getcolumn (std::vector<std::vector<double>> &m , int col){
    
    std::pair<int,int>a;

    a = get_rowcol(m);
    int n = a.first;    
    std::vector<double> v(n);

    for(int i=0;i<n;i++){
        v[i]=m[i][col];
    }

    return v;
}



//for square matrices only
std::vector<std::vector<double>> matmult(std::vector<std::vector<double>> &a,std::vector<std::vector<double>> &b){

    int n = a.size();
    std::vector<std::vector<double>> result(n,std::vector<double>(n));

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            double sum =0;
            for(int k=0; k<n ; k++){
                sum += a[i][k]*b[k][j];
            }
            result[i][j]=sum;
        }
    }

    return result;
}





std::vector<double> getrow (std::vector<std::vector<double>> &m , int row){
    
    int n = m.size();
    std::vector<double> v(n);

    for(int i=0;i<n;i++){
        v[i]=m[row][i];
    }

    return v;
}






void setcolumn (std:: vector<std:: vector<double>> &m , std::vector<double> &v, int col){
    
    int n = m.size();
    for(int i=0;i<n;i++){
        m[i][col] = v[i];
    }

}







std::vector<std::vector<double>> mat_inverse(std::vector<std::vector<double>>mat){

    //creating the identity

    int size = mat.size();
    std::vector<std::vector<double>>idm(size,std::vector<double>(size));


    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            if(i==j){
                idm[i][j]=1.0;
            }
            else{
                idm[i][j]=0.0;
            }
        }
    }


    for(int i=0;i<size;i++){

        //pivoting starts
        int pivotRow = i;
        double maxVal = std::abs(mat[i][i]);

        for(int k=i+1; k<size; k++){
            if(std::abs(mat[k][i]) > maxVal){
                maxVal = std::abs(mat[k][i]);
                pivotRow = k;
            }
        }

        if(pivotRow != i){
            std::swap(mat[i], mat[pivotRow]);
            std::swap(idm[i], idm[pivotRow]);
        }
        
        if(std::abs(mat[i][i]) < 1e-9){
            std::cerr << "Warning: Matrix is singular or close to it!" << std::endl;
        }

        double pivot = mat[i][i];


        //pivoting ends

        //divide all with the diagonal
        for(int j=0;j<size;j++){
            mat[i][j] /= pivot;
            idm[i][j] /= pivot;
        }

        //eliminate otherr rows

        for(int k=0;k<size;k++){
            if(k!=i){
                double factor = mat[k][i];

                for(int j=0;j<size;j++){
                    mat[k][j] -= factor*mat[i][j];
                    idm[k][j] -= factor*idm[i][j];
                }
            }
        }
    }


    return idm ;

}







#endif