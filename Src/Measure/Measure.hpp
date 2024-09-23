//
//  Gates.hpp
//
//

#ifndef MEASURE_HPP
#define MEASURE_HPP

#include <vector>
#include <map> 
#include <complex>

namespace NeoQuant{
    class Measure {
        public:
            Measure(std::vector<std::complex<double>> state); 
            // n-d vector with measured values. 
            std::vector<int> measureState(); 



        private:

            std::map<std::vector<int>, std::vector<std::complex<double>>> createQubitResultSpace(int n); 
            std::map<std::vector<int>, double> createQubitResultProbs(std::vector<std::complex<double>> state, int n); 

            int RNG(const std::vector<double> probs);



            std::vector<std::complex<double>> state; 
            int n; 
            std::vector<std::vector<int>> allStates; 

            std::map<std::vector<int>, std::vector<std::complex<double>>> qubitResultSpace; 
            std::map<std::vector<int>, double> qubitResultProbs; 

            // So this is numerical:
            // 0 0 0 
            // 0 1 0
            // etc. 
            std::vector<std::vector<int>> allPossibleQubitStates(int n); 
            void allPossibleBitStrings(std::vector<std::vector<int>>& allPossibleQubitStates, std::string s, int digitsLeft);
            

    };

}

#endif /* GATE_HPP */