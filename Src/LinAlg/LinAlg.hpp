//
//  LinAlg.hpp
//
//

// Taken from my machine learning library @ github.com/novak-99/MLPP

#ifndef LINALG_HPP
#define LINALG_HPP

#include <vector>
#include <complex>

namespace NeoQuant{
    class LinAlg{
        public:
        
        // MATRIX FUNCTIONS

        static std::vector<std::vector<std::complex<double>>> addition(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B);

        static std::vector<std::vector<std::complex<double>>> subtraction(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B);
        
        static std::vector<std::vector<std::complex<double>>> matmult(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B);
        
        static std::vector<std::vector<std::complex<double>>> hadamardProduct(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B);

        static std::vector<std::vector<std::complex<double>>> kroneckerProduct(std::vector<std::vector<std::complex<double>>> A, std::vector<std::vector<std::complex<double>>> B);
        
        static std::vector<std::vector<std::complex<double>>> transpose(std::vector<std::vector<std::complex<double>>> A);
        
        static std::vector<std::vector<std::complex<double>>> scalarMultiply(std::complex<double> scalar, std::vector<std::vector<std::complex<double>>> A);

        static std::vector<std::vector<std::complex<double>>> scalarAdd(std::complex<double> scalar, std::vector<std::vector<std::complex<double>>> A);

        static std::vector<std::vector<std::complex<double>>> matrixPower(std::vector<std::vector<std::complex<double>>> A, int n);

        static std::vector<std::vector<std::complex<double>>> abs(std::vector<std::vector<std::complex<double>>> A);

        static std::complex<double> trace(std::vector<std::vector<std::complex<double>>> A); 

        static std::vector<std::vector<std::complex<double>>> identity(int d);

        static std::vector<std::complex<double>> flatten(std::vector<std::vector<std::complex<double>>> A);
        
        static void printMatrix(std::vector<std::vector<std::complex<double>>> A);
        
        // VECTOR FUNCTIONS

        static std::vector<std::vector<std::complex<double>>> outerProduct(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b); // This multiplies a, bT 

        static std::vector<std::complex<double>> elementWiseDivision(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);
        
        static std::vector<std::complex<double>> scalarMultiply(std::complex<double> scalar, std::vector<std::complex<double>> a);

        static std::vector<std::complex<double>> scalarAdd(std::complex<double> scalar, std::vector<std::complex<double>> a);
        
        static std::vector<std::complex<double>> addition(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);
        
        static std::vector<std::complex<double>> subtraction(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);
        
        static std::complex<double> dot(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);

        static std::vector<std::complex<double>> kroneckerProduct(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);

        static std::vector<std::complex<double>> abs(std::vector<std::complex<double>> a);
        
        static void printVector(std::vector<std::complex<double>> a);

        static std::vector<std::complex<double>> matVecMult(std::vector<std::vector<std::complex<double>>> A, std::vector<std::complex<double>> b);
    
        private:
    };

}

#endif /* LINALG_HPP */