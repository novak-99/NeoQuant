//
//  States.hpp
//
//

#ifndef STATES_HPP
#define STATES_HPP

#include <vector>
#include <complex>

namespace NeoQuant{
    class States {
        public:
            States();
            std::vector<std::complex<double>> createState(double theta, double phi);
            std::vector<std::complex<double>> createMultiQubitState(std::vector<int> qubits);

            std::vector<std::complex<double>> getZeroState();
            std::vector<std::complex<double>> getOneState();

            std::vector<std::complex<double>> getPhiPlus();
            std::vector<std::complex<double>> getPhiMinus();
            std::vector<std::complex<double>> getPsiPlus();
            std::vector<std::complex<double>> getPsiMinus();

            double zeroProb(std::vector<std::complex<double>> state);
            double oneProb(std::vector<std::complex<double>> state);
            double stateProb(std::vector<std::complex<double>> state, std::vector<std::complex<double>> qubitResult);

            std::vector<std::complex<double>> createUniformZeroState(int n); 
            std::vector<std::complex<double>> createUniformOneState(int n); 

        private:
            std::vector<std::complex<double>> zeroState;
            std::vector<std::complex<double>> oneState;


            std::vector<std::complex<double>> phiPlus;
            std::vector<std::complex<double>> phiMinus;

            std::vector<std::complex<double>> psiPlus;
            std::vector<std::complex<double>> psiMinus;
    };

}

#endif /* STATES_HPP */