//
//  LinAlg.cpp
//
//

// Taken from my machine learning library @ github.com/novak-99/MLPP

#include "LinAlg.hpp"
#include <iostream>
#include <map>
#include <cmath>

namespace NeoQuant{

    std::vector<std::vector<std::complex<double>>> LinAlg::addition(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B){
        std::vector<std::vector<std::complex<double>>> C;
        C.resize(A.size());
        for(int i = 0; i < C.size(); i++){
            C[i].resize(A[0].size());
        }
        
        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[0].size(); j++){
                C[i][j] = A[i][j] + B[i][j];
            }
        }
        return C;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::subtraction(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B){
        std::vector<std::vector<std::complex<double>>> C;
        C.resize(A.size());
        for(int i = 0; i < C.size(); i++){
            C[i].resize(A[0].size());
        }

        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[0].size(); j++){
                C[i][j] = A[i][j] - B[i][j];
            }
        }
        return C;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::matmult(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B){
        std::vector<std::vector<std::complex<double>>> C;
        C.resize(A.size());
        for(int i = 0; i < C.size(); i++){
            C[i].resize(B[0].size());
        }
        
        for(int i = 0; i < A.size(); i++){ 
            for(int k = 0; k < B.size(); k++){ 
                for(int j = 0; j < B[0].size(); j++){ 
                    C[i][j] += A[i][k] * B[k][j]; 
                } 
            } 
        } 
        return C;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::hadamardProduct(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B){
        std::vector<std::vector<std::complex<double>>> C;
        C.resize(A.size());
        for(int i = 0; i < C.size(); i++){
            C[i].resize(A[0].size());
        }
        
        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[0].size(); j++){
                C[i][j] = A[i][j] * B[i][j];
            }
        }
        return C;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::kroneckerProduct(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B){
        std::vector<std::vector<std::complex<double>>> C;

        // [1,1,1,1]   [1,2,3,4,5]
        // [1,1,1,1]   [1,2,3,4,5]    
        //             [1,2,3,4,5]

        // [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5]
        // [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5]
        // [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5]
        // [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5]
        // [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5]
        // [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5] [1,2,3,4,5]

        // Resulting matrix: A.size() * B.size()
        //                   A[0].size() * B[0].size()

        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < B.size(); j++){
                std::vector<std::vector<std::complex<double>>> row;
                for(int k = 0; k < A[0].size(); k++){
                    row.push_back(scalarMultiply(A[i][k], B[j]));
                } 
                C.push_back(LinAlg::flatten(row));
            }
        }
        return C;    
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::transpose(std::vector<std::vector<std::complex<double>>> A){
        std::vector<std::vector<std::complex<double>>> AT;
        AT.resize(A[0].size());
        for(int i = 0; i < AT.size(); i++){
            AT[i].resize(A.size());
        }
        
        for(int i = 0; i < A[0].size(); i++){
            for(int j = 0; j < A.size(); j++){
                AT[i][j] = A[j][i];
            }
        }
        return AT;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::scalarMultiply(std::complex<double> scalar, std::vector<std::vector<std::complex<double>>> A){
        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[i].size(); j++){
                A[i][j] *= scalar;
            }
        }
        return A;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::scalarAdd(std::complex<double> scalar, std::vector<std::vector<std::complex<double>>> A){
        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[i].size(); j++){
                A[i][j] += scalar;
            }
        }
        return A;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::matrixPower(std::vector<std::vector<std::complex<double>>> A, int n){
        std::vector<std::vector<std::complex<double>>> B = identity(A.size());
        if(n == 0){
            return identity(A.size());
        }
        // don't need this branch...
        // else if(n < 0){
        //     A = inverse(A);
        // }
        for(int i = 0; i < std::abs(n); i++){
            B = matmult(B, A);
        }
        return B;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::abs(std::vector<std::vector<std::complex<double>>> A){
        std::vector<std::vector<std::complex<double>>> B;
        B.resize(A.size());
        for(int i = 0; i < B.size(); i++){
            B[i].resize(A[0].size());
        }
        for(int i = 0; i < B.size(); i++){
            for(int j = 0; j < B[i].size(); j++){
                B[i][j] = std::abs(A[i][j]);
            }
        }
        return B;
    }

    std::complex<double> LinAlg::trace(std::vector<std::vector<std::complex<double>>> A){
        std::complex<double> trace = 0;
        for(int i = 0; i < A.size(); i++){
            trace += A[i][i];
        }
        return trace;
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::identity(int d){
        std::vector<std::vector<std::complex<double>>> identityMat; 
        identityMat.resize(d);
        for(int i = 0; i < identityMat.size(); i++){
            identityMat[i].resize(d);
        }
        for(int i = 0; i < identityMat.size(); i++){
            for(int j = 0; j < identityMat.size(); j++){
                if(i == j){
                    identityMat[i][j] = 1;
                }
                else { identityMat[i][j] = 0; }
            }
        }
        return identityMat;
    }

    std::vector<std::complex<double>> LinAlg::flatten(std::vector<std::vector<std::complex<double>>> A){
        std::vector<std::complex<double>> a; 
        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[i].size(); j++){
                a.push_back(A[i][j]);
            }
        }
        return a;
    }

    void LinAlg::printMatrix(std::vector<std::vector<std::complex<double>>> A){
        for(int i = 0; i < A.size(); i++){
            for(int j = 0; j < A[i].size(); j++){
                std::cout << A[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    std::vector<std::vector<std::complex<double>>> LinAlg::outerProduct(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b){
        std::vector<std::vector<std::complex<double>>> C;
        C.resize(a.size());
        for(int i = 0; i < C.size(); i++){
            C[i] = scalarMultiply(a[i], b);
        }
        return C;
    }

    std::vector<std::complex<double>> LinAlg::scalarMultiply(std::complex<double> scalar, std::vector<std::complex<double>> a){
        for(int i = 0; i < a.size(); i++){
            a[i] *= scalar;
        }
        return a;
    }

    std::vector<std::complex<double>> LinAlg::scalarAdd(std::complex<double> scalar, std::vector<std::complex<double>> a){
        for(int i = 0; i < a.size(); i++){
            a[i] += scalar;
        }
        return a;
    }

    std::vector<std::complex<double>> LinAlg::addition(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b){
        std::vector<std::complex<double>> c;
        c.resize(a.size());
        for(int i = 0; i < a.size(); i++){
            c[i] = a[i] + b[i];
        }
        return c;
    }

    std::vector<std::complex<double>> LinAlg::subtraction(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b){
        std::vector<std::complex<double>> c;
        c.resize(a.size());
        for(int i = 0; i < a.size(); i++){
            c[i] = a[i] - b[i];
        }
        return c;
    }

    std::complex<double> LinAlg::dot(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b){
        std::complex<double> c = 0;
        for(int i = 0; i < a.size(); i++){
            c += a[i] * b[i];
        }
        return c;
    }

    std::vector<std::complex<double>> LinAlg::kroneckerProduct(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b){
        std::vector<std::complex<double>> c;
        c.resize(a.size() * b.size()); 

        for(int i = 0; i < a.size(); i++){
            std::vector batch = LinAlg::scalarMultiply(a[i], b);
            for(int j = 0; j < batch.size(); j++){
                c[i * b.size() + j] = batch[j];
            }
        }
        return c;
    }

    std::vector<std::complex<double>> LinAlg::abs(std::vector<std::complex<double>> a){
        std::vector<std::complex<double>> b; 
        b.resize(a.size());
        for(int i = 0; i < b.size(); i++){
            b[i] = std::abs(a[i]);
        }
        return b;
    }

    void LinAlg::printVector(std::vector<std::complex<double>> a){
        for(int i = 0; i < a.size(); i++){
            std::cout << a[i] << " ";
        }
        std::cout << std::endl;
    }

    std::vector<std::complex<double>> LinAlg::matVecMult(std::vector<std::vector<std::complex<double>>> A, std::vector<std::complex<double>> b){
        std::vector<std::complex<double>> c;
        c.resize(A.size());
            
        for(int i = 0; i < A.size(); i++){
            for(int k = 0; k < b.size(); k++){
                c[i] += A[i][k] * b[k];
            }
        }
        return c;
    }
}