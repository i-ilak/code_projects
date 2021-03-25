#include "parameters_simulation.hpp"
#include "naive_sol.hpp"
#include "ode_int_sol.hpp"

int main(){
    // Generate vector of StellarObjects to simulate time evolution
    // and set up all the other parameters, where T is the simulated time 
    // interval and N the number of time-steps we will take in the simulation.
    auto objects = solar_system();
    double const G = 2.95912208286e-4;
    int const N = 40000;
    double const T = 200000;
    
    // Do simulation using naive integrators (Velocity Verlet)
    perform_simulation(objects, T, N, G, 3);

    // Do simulation using ODE int
    perform_ode_int_simulation(objects, T, N, G);

}
