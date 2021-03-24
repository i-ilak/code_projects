#include "fishing_pop.hpp"
#include <iostream>
#include <chrono>
#include <fstream>

int main(){
    std::string const path="data.txt";
    std::ofstream out(path);

    Penna::Genome::set_mutation_rate(2);
    Penna::Animal::set_threshold(2);
    Penna::Animal::set_reproduction_age(8);
    Penna::Animal::set_prob_to_get_pregnant(1);

    Penna::FishingPopulation fish(300, 10000, 0, 0);

    auto start = std::chrono::high_resolution_clock::now();
    for(std::size_t k = 0; k < 5000; ++k ){
        fish.step();
        out << k << "\t" << fish.size() << "\n";
        if(k == 500){
            fish.change_fishing(0.17, 8);
        }
        if(k == 3500){
            fish.change_fishing(0.22, 0);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Runtime: " << diff.count() << " s" << "\n";
    out.close();

    return 0;
}