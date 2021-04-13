#include "include/parameters_simulation.hpp"
#include "include/perform_simulation.hpp"
#include <omp.h>

using planet_container_type = std::vector<StellarObject>;

int main(){
    // We want to run some of the code in parallel, so we set the number of available processors.
    int nProcessors = omp_get_max_threads();

    std::cout<<"Number of available processors:\t" << nProcessors<<std::endl;
    omp_set_num_threads(nProcessors);           

    // Generate vector of StellarObjects to simulate time evolution
    // and set up all the other parameters, where T is the simulated time 
    // interval and N the number of time-steps we will take in the simulation.
    planet_container_type objects = solar_system();
    double const G = 2.95912208286e-4;
    int const N = 100000;
    double const T = 200000;
    
    // Time simulation
    auto start = std::chrono::high_resolution_clock::now();
    // Do simulation in parallel for all integrators
    #pragma omp parallel for
    for (int i =0; i<=3;i++){
        if(i==3){
            perform_ode_int_simulation(objects, T, N);
        }
        else{
            perform_simulation(objects, T, N, G, i+1);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    std::cout << "Entire calculation took:\t" << diff.count() << " s" <<"\n";

}
