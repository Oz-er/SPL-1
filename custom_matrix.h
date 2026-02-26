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



std::vector<std::vector<double>> mat_scalar_mult(std::vector<std::vector<double>> &mat, double scalar){
    int rows = mat.size();
    int cols = mat[0].size();
    std::vector<std::vector<double>> result(rows,std::vector<double>(cols));

    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            result[i][j] = mat[i][j] * scalar;
        }
    }
    return result;
}



std::vector<std::vector<double>> matmult(const std::vector<std::vector<double>> &a,const std::vector<std::vector<double>> &b){

    if(a.empty() || b.empty())return{};

    int rows_a = a.size();
    int cols_a = a[0].size();
    int cols_b = b[0].size();

    // int n = a.size();
    std::vector<std::vector<double>> result(rows_a,std::vector<double>(cols_b,0.0));

    for(int i=0;i<rows_a;i++){
        for(int j=0;j<cols_b;j++){
            double sum =0;
            for(int k=0; k<cols_a;k++){
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


std::vector<std::vector<double>>matadd(std::vector<std::vector<double>> &a,std::vector<std::vector<double>> &b){
    int rows=a.size();
    int cols=a[0].size();

    std::vector<std::vector<double>> result(rows,std::vector<double>(cols));

    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            result[i][j]=a[i][j] + b[i][j];
        }
    }

    return result;
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
            std::cerr << "Warning: Matrix is singular or close to it" << std::endl;
        }

        double pivot = mat[i][i];


        //pivoting ends

        //divide all with the pivot(dividing by a large pivot keeps the error shrunk)
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