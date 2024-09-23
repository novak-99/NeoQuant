//
//  Circuits.hpp
//
//

#ifndef CIRCUITS_HPP
#define CIRCUITS_HPP

#include <vector>
#include <complex>

#include <tuple>
#include <utility>

namespace NeoQuant{
    class Circuits {
        public:
            std::tuple<int, int, int> quantumTeleportation(std::vector<std::complex<double>> Alice, std::vector<std::complex<double>> Bob, std::vector<std::complex<double>> state);
            bool deutschJosza(std::vector<int> functionOutputs);
        private: 
            std::pair<int, int> measureAliceAndTeleport(std::vector<std::complex<double>> netState);
            int measureBob(std::vector<std::complex<double>> netState, int i, int j);

            // We will move this function to a "measuring" class at a later point. 
            int RNG(const std::vector<double> probs);

    };

}

#endif /* CIRCUITS_HPP */