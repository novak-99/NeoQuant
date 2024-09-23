//
//  States.cpp
//
//

#include <iostream>
#include <cmath>
#include "States/States.hpp"
#include "LinAlg/LinAlg.hpp"

using namespace std::complex_literals;

namespace NeoQuant {

    States::States()
    : 
    zeroState({
        1,
        0
    }), 

    oneState({
        0,
        1
    }),

    phiPlus({
        1/std::sqrt(2),
        0,
        0,
        1/std::sqrt(2),
    }),

    phiMinus({
        1/std::sqrt(2),
        0,
        0,
        -1/std::sqrt(2),
    }),

    psiPlus({
        0,
        1/std::sqrt(2),
        1/std::sqrt(2),
        0,
    }),

    psiMinus({
        0,
        1/std::sqrt(2),
        -1/std::sqrt(2),
        0,
    })


    {

    }
    
    std::vector<std::complex<double>> States::createState(double theta, double phi){
        std::complex<double> i(0.0,1.0); 
        return LinAlg::addition(LinAlg::scalarMultiply(std::cos(theta/2), zeroState), LinAlg::scalarMultiply(std::exp(phi * i) * std::sin(theta/2), oneState));
    }

    std::vector<std::complex<double>> States::createMultiQubitState(std::vector<int> qubits){

        std::vector<std::complex<double>> currentState = !qubits[0] ? zeroState : oneState; 

        
        for(int i = 1; i < qubits.size(); i++){
            if(qubits[i] == 0){
                std::vector<std::complex<double>> newState = LinAlg::kroneckerProduct(currentState, zeroState);
                currentState = newState;
            }
            else{ // qubits_i == 1. 
                std::vector<std::complex<double>> newState = LinAlg::kroneckerProduct(currentState, oneState);
                currentState = newState;
            }
        }
        return currentState;
    }

    std::vector<std::complex<double>> States::getZeroState(){
        return zeroState;
    }

    std::vector<std::complex<double>> States::getOneState(){
        return oneState;
    }

    std::vector<std::complex<double>> States::getPhiPlus(){
        return phiPlus; 
    }

    std::vector<std::complex<double>> States::getPhiMinus(){
        return phiMinus;
    }

    std::vector<std::complex<double>> States::getPsiPlus(){
        return psiPlus; 
    }

    std::vector<std::complex<double>> States::getPsiMinus(){
        return psiMinus;
    }

    double States::zeroProb(std::vector<std::complex<double>> state){
        return std::real(std::pow(LinAlg::dot(zeroState, state), 2));
    }

    double States::oneProb(std::vector<std::complex<double>> state){
        return std::real(std::pow(LinAlg::dot(oneState, state), 2)); 
    }

    double States::stateProb(std::vector<std::complex<double>> state, std::vector<std::complex<double>> qubitResult){
        return std::real(std::pow(LinAlg::dot(state, qubitResult), 2)); 
    }

    std::vector<std::complex<double>> States::createUniformZeroState(int n){

        std::vector<std::complex<double>> state = zeroState;

        for(int i = 1; i < n; i++){
            state = LinAlg::kroneckerProduct(state, zeroState);
        }

        return state; 
    }

    std::vector<std::complex<double>> States::createUniformOneState(int n){

        std::vector<std::complex<double>> state = oneState;

        for(int i = 1; i < n; i++){
            state = LinAlg::kroneckerProduct(state, oneState);
        }

        return state; 
    }
}