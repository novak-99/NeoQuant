//
//  Gates.cpp
//
//

#include <iostream>
#include <cmath>
#include "Gates/Gates.hpp"
#include "LinAlg/LinAlg.hpp"

 using namespace std::complex_literals;

namespace NeoQuant {

    Gates::Gates()
    : 
    X({
        {0, 1}, 
        {1, 0}
    }), 

    Y({
        {0, -1i},
        {1i, 0}
    }),

    Z({
        {1, 0},
        {0, -1}
    }),

    hadamard({
        {1/std::sqrt(2), 1/std::sqrt(2)},
        {1/std::sqrt(2), -1/std::sqrt(2)}
    }),

    CNOT({
        {1, 0, 0, 0},
        {0, 1, 0, 0}, 
        {0, 0, 0, 1},
        {0, 0, 1, 0}
    })

    {

    }


    std::vector<std::complex<double>> Gates::entangle(std::vector<std::complex<double>> a, std::vector<std::complex<double>> b){
        return applyCNOT(LinAlg::kroneckerProduct(applyHadamard(a), b)); 
    }

    std::vector<std::complex<double>> Gates::applyUniformGate(std::vector<std::complex<double>> state, std::vector<std::vector<std::complex<double>>> gate){
        int numQubits = std::log2(state.size());

        std::vector<std::vector<std::complex<double>>> identity = LinAlg::identity(2);

        std::vector<std::vector<std::complex<double>>> extendedGate = {{1}}; 

        for(int i = 1; i <= numQubits; i++){
            extendedGate = LinAlg::kroneckerProduct(extendedGate, gate); 
        }

        return LinAlg::matVecMult(extendedGate, state);
    }

    std::vector<std::complex<double>> Gates::applySpecifiedGate(std::vector<std::complex<double>> state, std::vector<int> qubits, std::vector<std::vector<std::complex<double>>> gate){
        int numQubits = std::log2(state.size());

        std::vector<std::vector<std::complex<double>>> identity = LinAlg::identity(2);

        std::vector<std::vector<std::complex<double>>> extendedGate = {{1}}; 

        for(int i = 1; i <= numQubits; i++){
            if(std::find(qubits.begin(), qubits.end(), i) != qubits.end()){
                extendedGate = LinAlg::kroneckerProduct(extendedGate, gate); 
            }
            else{
                extendedGate = LinAlg::kroneckerProduct(extendedGate, identity); 
            }
        }

        return LinAlg::matVecMult(extendedGate, state);
    }

    // Apply to all consecutive gates. 
    std::vector<std::complex<double>> Gates::applyUniform4DGate(std::vector<std::complex<double>> state, std::vector<std::vector<std::complex<double>>> gate){
        int numQubits = std::log2(state.size());

        std::vector<std::vector<std::complex<double>>> identity = LinAlg::identity(2);

        std::vector<std::vector<std::complex<double>>> extendedGate = {{1}}; 

        for(int i = 1; i < numQubits; i+=2){
            extendedGate = LinAlg::kroneckerProduct(extendedGate, gate); 
        }

        // Well, if odd num of qubits, the last qubit will need to have an identity applied 
        if(numQubits % 2 != 0){
            extendedGate = LinAlg::kroneckerProduct(extendedGate, identity);
        }

        return LinAlg::matVecMult(extendedGate, state);
    }
    
    std::vector<std::complex<double>> Gates::applySpecified4DGate(std::vector<std::complex<double>> state, std::vector<std::pair<int, int>> qubits, std::vector<std::vector<std::complex<double>>> gate){
        int numQubits = std::log2(state.size());

        std::vector<std::vector<std::complex<double>>> identity = LinAlg::identity(2);

        std::vector<std::vector<std::complex<double>>> extendedGate = {{1}}; 

        for(int i = 1; i <= numQubits; i++){
            if(std::find(qubits.begin(), qubits.end(), std::pair<int,int>({i, i+1})) != qubits.end()){
                extendedGate = LinAlg::kroneckerProduct(extendedGate, gate);
                i++;  
            }
            else{
                extendedGate = LinAlg::kroneckerProduct(extendedGate, identity); 
            }
        }
        return LinAlg::matVecMult(extendedGate, state);
    }

    std::vector<std::complex<double>> Gates::applyX(std::vector<std::complex<double>> state){
        return applyUniformGate(state, X);
    }

    std::vector<std::complex<double>> Gates::applyX(std::vector<std::complex<double>> state, int qubit){
        return applySpecifiedGate(state, {qubit}, X);
    }

    std::vector<std::complex<double>> Gates::applyX(std::vector<std::complex<double>> state, std::vector<int> qubits){
        return applySpecifiedGate(state, qubits, X);
    }

    std::vector<std::complex<double>> Gates::applyY(std::vector<std::complex<double>> state){
        return applyUniformGate(state, Y);
    }

    std::vector<std::complex<double>> Gates::applyY(std::vector<std::complex<double>> state, int qubit){
        return applySpecifiedGate(state, {qubit}, Y);
    }

    std::vector<std::complex<double>> Gates::applyY(std::vector<std::complex<double>> state, std::vector<int> qubits){
        return applySpecifiedGate(state, qubits, Y);
    }

    std::vector<std::complex<double>> Gates::applyZ(std::vector<std::complex<double>> state){
        return applyUniformGate(state, Z);
    }

    std::vector<std::complex<double>> Gates::applyZ(std::vector<std::complex<double>> state, int qubit){
        return applySpecifiedGate(state, {qubit}, Z);
    }

    std::vector<std::complex<double>> Gates::applyZ(std::vector<std::complex<double>> state, std::vector<int> qubits){
        return applySpecifiedGate(state, qubits, Z);
    }

    std::vector<std::complex<double>> Gates::applyHadamard(std::vector<std::complex<double>> state){
        return applyUniformGate(state, hadamard);
    }

    std::vector<std::complex<double>> Gates::applyHadamard(std::vector<std::complex<double>> state, int qubit){
        return applySpecifiedGate(state, {qubit}, hadamard);
    }

    std::vector<std::complex<double>> Gates::applyHadamard(std::vector<std::complex<double>> state, std::vector<int> qubits){
        return applySpecifiedGate(state, qubits, hadamard);
    }

    std::vector<std::complex<double>> Gates::applyCNOT(std::vector<std::complex<double>> state){
        return applyUniform4DGate(state, CNOT); 
    }

    std::vector<std::complex<double>> Gates::applyCNOT(std::vector<std::complex<double>> state, std::pair<int, int> qubits){
        return applySpecified4DGate(state, {qubits}, CNOT);
    }
    
    std::vector<std::complex<double>> Gates::applyCNOT(std::vector<std::complex<double>> state, std::vector<std::pair<int, int>> qubits){
        return applySpecified4DGate(state, qubits, CNOT);
    }
    
    std::vector<std::vector<std::complex<double>>> Gates::getX(){
        return X; 
    }

    std::vector<std::vector<std::complex<double>>> Gates::getY(){
        return Y; 
    }

    std::vector<std::vector<std::complex<double>>> Gates::getZ(){
        return Z; 
    }

    std::vector<std::vector<std::complex<double>>> Gates::getHadamard(){
        return hadamard;
    }

    std::vector<std::vector<std::complex<double>>> Gates::getCNOT(){
        return CNOT; 
    }

    std::vector<std::complex<double>> Gates::phaseOracle(std::vector<std::complex<double>> state, std::vector<int> functionOutputs){
        int size = state.size(); // 2^N size vector. 
        std::vector<std::vector<std::complex<double>>> phaseOracleGate = LinAlg::identity(size);

        for(int i = 0; i < functionOutputs.size(); i++){
            phaseOracleGate[i][i] = std::pow(phaseOracleGate[i][i], functionOutputs[i]);  // 0 =>  no change, so 1, 1 => switch to -1. 
        }

        return LinAlg::matVecMult(phaseOracleGate, state);

    }


}