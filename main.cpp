#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "Circuits/Circuits.hpp"
#include "States/States.hpp"
#include "LinAlg/LinAlg.hpp"
#include "Gates/Gates.hpp"
#include "Measure/Measure.hpp"
using namespace NeoQuant; 
std::vector<std::string> answer;

int main(){
    
    // States states; 
    // LinAlg::printVector(LinAlg::kroneckerProduct(states.getOneState(), LinAlg::kroneckerProduct(states.getOneState(), states.getOneState())));

    States states; 
    Circuits circuits; 
    std::vector<std::complex<double>> teleportationQubit = LinAlg::addition(LinAlg::scalarMultiply(std::sqrt(0.7), states.getZeroState()), LinAlg::scalarMultiply(std::sqrt(0.3), states.getOneState()));


    const int N = 1000; 
    double k = 0; 
    for(int i = 0; i < N; i++){
        auto [a, b, c] = circuits.quantumTeleportation(states.getZeroState(), states.getZeroState(), teleportationQubit);

        if(c == 1) k++; 
    }

    std::cout << double(k) / double(N) << "\n";

    // std::vector<std::complex<double>> exampleState = {1,2,3,4,5,6,7,8, 1,2,3,4,5,6,7,8}; 

    // Gates gates; 
    // LinAlg::printVector(gates.applyCNOT(exampleState, {{1,2}, {3,4}} ));


    // auto [a, b, c] = circuits.quantumTeleportation(states.getZeroState(), states.getZeroState(), teleportationQubit);
    // std::cout << a << b << c << "\n";

    // return 0; 




    // Circuits circuits; 
    // std::cout << circuits.deutschJosza({0,0,1,0,1,1}) << "\n";
}