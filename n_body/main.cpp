#include <string>
#include <utility>
#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include "integrate.hpp"
#include "n_body_solver.hpp"
#include "parameters_simulation.hpp"


int main(){

    // Generate vector of StellarObjects to simulate time evolution
    // and set up all the other parameters, where T is the simulated time 
    // interval and N the number of time-steps we will take in the simulation.
    auto objects = solar_system();
    double const G = 2.95912208286e-4;
    int const N = 40000;
    double const T = 60000;

    // Running the calculation and timing it
    auto start = std::chrono::high_resolution_clock::now();
    auto res = n_body_solver(objects, T, N, G);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "The relevant calculation took:\t" << diff.count() << " s" <<"\n";

    // Writing the data into a file for later plotting
    std::string path="data.txt";
    std::ofstream out(path);

    out << objects.size() << "\n";
    for(auto stellar_obj : objects){
        out << stellar_obj.get_name() << ",";
    }
    out << "\n";
    auto start2 = std::chrono::high_resolution_clock::now();
    out << res << "\n";
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2 = end2 - start2;
    std::cout << "Time to write to output-file:\t" << diff2.count() << " s" <<"\n";

}
