#include "parameters_simulation.hpp"
#include "naive_sol.hpp"
#include "ode_int_sol.hpp"

using planet_container_type = std::vector<StellarObject>;

int main(){
    // Generate vector of StellarObjects to simulate time evolution
    // and set up all the other parameters, where T is the simulated time 
    // interval and N the number of time-steps we will take in the simulation.
    planet_container_type objects = two_body();
    double const G = 2.95912208286e-4;
    int const N = 40000;
    double const T = 200000;
    
    // Do simulation using naive integrators
    perform_simulation(objects, T, N, G, 1);
    perform_simulation(objects, T, N, G, 2);
    perform_simulation(objects, T, N, G, 3);

    // Do simulation using ODE int
    perform_ode_int_simulation(objects, T, N);

}
