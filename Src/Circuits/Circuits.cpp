//
//  Circuits.hpp
//
//

#include "Circuits/Circuits.hpp"
#include "LinAlg/LinAlg.hpp"
#include "Gates/Gates.hpp"
#include "States/States.hpp"
#include "Measure/Measure.hpp"
#include <cmath>
#include <cstdlib>    

#include <random>

#include <iostream>

namespace NeoQuant{

    std::tuple<int, int, int> Circuits::quantumTeleportation(std::vector<std::complex<double>> Alice, std::vector<std::complex<double>> Bob, std::vector<std::complex<double>> state){
        Gates gates; 
        States states;
        // Entangle Alice & Bob 

        std::vector<std::complex<double>> entangledState = gates.entangle(Alice, Bob); 
        //std::vector<std::complex<double>> entangledState = states.getPhiPlus();

        // Include teleportation qubit in total state. 
        std::vector<std::complex<double>> netState = LinAlg::kroneckerProduct(state, entangledState);

        // Bell measurement of Alice + teleportation qubit

        netState = gates.applyCNOT(netState, {1, 2}); 
        netState = gates.applyHadamard(netState, 1);

    

        std::pair<int, int> measuredQubits = measureAliceAndTeleport(netState);

        int i = measuredQubits.first; 
        int j = measuredQubits.second;

        if(j == 1) netState = gates.applyX(netState, 3);
        if(i == 1) netState = gates.applyZ(netState, 3);


        return {i, j, measureBob(netState, i, j)}; 

    }

    bool Circuits::deutschJosza(std::vector<int> functionOutputs){
        Gates gates; 
        States states; 

        std::vector<std::complex<double>> state = states.createMultiQubitState(functionOutputs);



        state = gates.applyHadamard(state); // uniform hadamard. 

        state = gates.phaseOracle(state, functionOutputs);

        state = gates.applyHadamard(state);
        

        // Final NUMERICAL, not quantum, state. 
        Measure measure(state);

        std::vector<int> finalBitString = measure.measureState();

        // You could use a for loop, but this is nicer. 
        // https://stackoverflow.com/questions/20132485/check-if-entire-vector-is-zero
        bool allZero = std::all_of(finalBitString.begin(), finalBitString.end(), [](int i) { return i==0; }); 
        bool allOne = std::all_of(finalBitString.begin(), finalBitString.end(), [](int i) { return i==1; }); 

        return true ? (allZero || allOne) : false; 

        // Remark: 
        // Return TRUE = CONSTANT FUNCTION. [ all zeroes, all ones ]
        // Return FALSE = BALANCED FUNCTION. [ mixed ]
    }

    std::pair<int, int> Circuits::measureAliceAndTeleport(std::vector<std::complex<double>> netState){
        States states; 
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> qubitResultSpace; 
        // Make indexing a lot easier. 
        qubitResultSpace.resize(2); 
        for(int i = 0; i < 2; i++){
            qubitResultSpace[i].resize(2);
            for(int j = 0; j < 2; j++){
                qubitResultSpace[i][j].resize(2);
                for(int k = 0; k < 2; k++){
                    // Remark: phi1 [x] phi2 [x] phi3 yields an 8-D vector. 
                    qubitResultSpace[i][j][k].resize(8); 
                }
            }
        }

        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    qubitResultSpace[i][j][k] = states.createMultiQubitState({i, j, k});
                }
            }
        }

        // Create probabilities. 
        std::vector<std::vector<double>> qubitResultProbs; 
        qubitResultProbs.resize(2); 
        for(int i = 0; i < 2; i++){
            qubitResultProbs[i].resize(2);
        }

        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    qubitResultProbs[i][j] += states.stateProb(qubitResultSpace[i][j][k], netState); // squared dot yields probs. 
                }
            }
        }

        int e = RNG({qubitResultProbs[0][0], qubitResultProbs[0][1], qubitResultProbs[1][0],  qubitResultProbs[1][1]});

        if(e == 0){
            return {0, 0};
        }
        else if(e == 1){
            return {0, 1};
        }
        else if(e == 2){
            return {1, 0};
        }
        else{
            return {1, 1};
        }


    }

    int Circuits::measureBob(std::vector<std::complex<double>> netState, int i, int j){
        States states; 
        std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> qubitResultSpace; 
        // Make indexing a lot easier. 
        qubitResultSpace.resize(2); 
        for(int i = 0; i < 2; i++){
            qubitResultSpace[i].resize(2);
            for(int j = 0; j < 2; j++){
                qubitResultSpace[i][j].resize(2);
                for(int k = 0; k < 2; k++){
                    // Remark: phi1 [x] phi2 [x] phi3 yields an 8-D vector. 
                    qubitResultSpace[i][j][k].resize(8); 
                }
            }
        }

        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    qubitResultSpace[i][j][k] = states.createMultiQubitState({i, j, k});
                }
            }
        }

        // Create probabilities. 
        std::vector<std::vector<std::vector<double>>> qubitResultProbs; 
        qubitResultProbs.resize(2); 
        for(int i = 0; i < 2; i++){
            qubitResultProbs[i].resize(2);
            for(int j = 0; j < 2; j++){
                qubitResultProbs[i][j].resize(2);
            }
        }

        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    qubitResultProbs[i][j][k] += states.stateProb(qubitResultSpace[i][j][k], netState); // squared dot yields probs. 
                }
            }
        }

        double totalProb = qubitResultProbs[i][j][0] + qubitResultProbs[i][j][1];

        double probZero = qubitResultProbs[i][j][0] / totalProb;
        double probOne = qubitResultProbs[i][j][1] / totalProb;

        int e = RNG({ probZero, probOne });

        if(e == 0){
            return 0;
        }
        else{
            return 1; 
        }
    }

    // Redo later in "measuring"
    int Circuits::RNG(const std::vector<double> probs){
        // Generate a random number between 0 and 1
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double randomNum = dis(gen);

        // Accumulate the probabilities
        double cumulativeProbability = 0.0;
        for (int i = 0; i < probs.size(); ++i) {
            cumulativeProbability += probs[i];
            if (randomNum <= cumulativeProbability) {
                return i; // Return the index of the chosen event
            }
        }

        // This line should never be reached, but added to avoid compiler warnings
        return -1;
    }

}