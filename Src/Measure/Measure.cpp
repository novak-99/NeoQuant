//
//  Measure.cpp
//
//

#include <algorithm>
#include <random>
#include <string>
#include <iostream>

#include "Measure/Measure.hpp"
#include "States/States.hpp"

namespace NeoQuant {

    Measure::Measure(std::vector<std::complex<double>> state)
    : state(state), n(std::log2(state.size())), allStates(allPossibleQubitStates(n)), qubitResultSpace(createQubitResultSpace(n)), qubitResultProbs(createQubitResultProbs(state, n))
    {

    }

    // Easy peasy lemon squeezy 
    std::vector<int> Measure::measureState(){
        std::vector<double> probs; 
        for(auto const& [state, prob] : qubitResultProbs){
            std::cout << "probs :" << prob << "\n";
            probs.push_back(prob); 
        }
        std::cout << RNG(probs) << "\n";

        return allStates[RNG(probs)]; 
    }

    std::map<std::vector<int>, std::vector<std::complex<double>>> Measure::createQubitResultSpace(int n){
        States states; 
        //std::vector<std::vector<int>> allStates = allPossibleQubitStates(n);
        std::map<std::vector<int>, std::vector<std::complex<double>>> qubitResultSpace; 

        for(int i = 0; i < allStates.size(); i++){
            qubitResultSpace[allStates[i]] = states.createMultiQubitState(allStates[i]);
        }

        return qubitResultSpace; 
    }

    std::map<std::vector<int>, double> Measure::createQubitResultProbs(std::vector<std::complex<double>> state, int n){
        States states; 
        //std::vector<std::vector<int>> allStates = allPossibleQubitStates(n);
        std::map<std::vector<int>, double> qubitResultProbs; 

        for(int i = 0; i < allStates.size(); i++){
            qubitResultProbs[allStates[i]] = states.stateProb(qubitResultSpace[allStates[i]], state); // squared dot yields probs. 
        }

        return qubitResultProbs; 
    }



    std::vector<std::vector<int>> Measure::allPossibleQubitStates(int n){
        std::vector<std::vector<int>> allStates; 
        allPossibleBitStrings(allStates, "", n);
        return allStates; 
    }

    void Measure::allPossibleBitStrings(std::vector<std::vector<int>>& allStates, std::string s, int digitsLeft){
        if(digitsLeft == 0){ // the length of string is n
            std::vector<int> currentPerm; 
            for(int i = 0; i < s.length(); i++){
                currentPerm.push_back(int(s[i] - '0'));
            }
            allStates.push_back(currentPerm);
        }
        else {
            allPossibleBitStrings(allStates, s + "0", digitsLeft - 1);
            allPossibleBitStrings(allStates, s + "1", digitsLeft - 1);
        }
    }

    int Measure::RNG(const std::vector<double> probs){
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