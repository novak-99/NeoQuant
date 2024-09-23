//
//  Gates.hpp
//
//

#ifndef GATES_HPP
#define GATES_HPP

#include <vector>
#include <utility>
#include <complex>

namespace NeoQuant{
    class Gates {
        public:
            Gates();
            std::vector<std::complex<double>> entangle(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);
            std::vector<std::complex<double>> bellMeasurement(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b);

            // Uniform means all qubits 
            // Specified = to specific qubits. 
            std::vector<std::complex<double>> applyUniformGate(std::vector<std::complex<double>> state, std::vector<std::vector<std::complex<double>>> gate);
            std::vector<std::complex<double>> applySpecifiedGate(std::vector<std::complex<double>> state, std::vector<int> qubits, std::vector<std::vector<std::complex<double>>> gate); 

            // for 2x2 gates. 
            std::vector<std::complex<double>> applyUniform4DGate(std::vector<std::complex<double>> state, std::vector<std::vector<std::complex<double>>> gate);
            std::vector<std::complex<double>> applySpecified4DGate(std::vector<std::complex<double>> state, std::vector<std::pair<int, int>> qubits, std::vector<std::vector<std::complex<double>>> gate); 

            std::vector<std::complex<double>> applyX(std::vector<std::complex<double>> state); 
            std::vector<std::complex<double>> applyX(std::vector<std::complex<double>> state, int qubit); 
            std::vector<std::complex<double>> applyX(std::vector<std::complex<double>> state, std::vector<int> qubits); 

            std::vector<std::complex<double>> applyY(std::vector<std::complex<double>> state); 
            std::vector<std::complex<double>> applyY(std::vector<std::complex<double>> state, int qubit); 
            std::vector<std::complex<double>> applyY(std::vector<std::complex<double>> state, std::vector<int> qubits); 

            std::vector<std::complex<double>> applyZ(std::vector<std::complex<double>> state); 
            std::vector<std::complex<double>> applyZ(std::vector<std::complex<double>> state, int qubit); 
            std::vector<std::complex<double>> applyZ(std::vector<std::complex<double>> state, std::vector<int> qubits); 

            std::vector<std::complex<double>> applyHadamard(std::vector<std::complex<double>> state); 
            std::vector<std::complex<double>> applyHadamard(std::vector<std::complex<double>> state, int qubit); 
            std::vector<std::complex<double>> applyHadamard(std::vector<std::complex<double>> state, std::vector<int> qubits); 

            std::vector<std::complex<double>> applyCNOT(std::vector<std::complex<double>> state); 
            std::vector<std::complex<double>> applyCNOT(std::vector<std::complex<double>> state, std::pair<int, int> qubits); 
            std::vector<std::complex<double>> applyCNOT(std::vector<std::complex<double>> state, std::vector<std::pair<int, int>> qubits); 


            std::vector<std::vector<std::complex<double>>> getX();
            std::vector<std::vector<std::complex<double>>> getY();
            std::vector<std::vector<std::complex<double>>> getZ();

            std::vector<std::vector<std::complex<double>>> getHadamard();
            std::vector<std::vector<std::complex<double>>> getCNOT();


            std::vector<std::complex<double>> phaseOracle(std::vector<std::complex<double>> state, std::vector<int> functionOutputs);
        private:
            std::vector<std::vector<std::complex<double>>> X; 
            std::vector<std::vector<std::complex<double>>> Y;
            std::vector<std::vector<std::complex<double>>> Z;

            std::vector<std::vector<std::complex<double>>> hadamard; 
            std::vector<std::vector<std::complex<double>>> CNOT; 
    };

}

#endif /* GATES_HPP */