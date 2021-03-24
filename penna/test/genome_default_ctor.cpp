#include "genome.hpp"
#include <stdexcept>
#include <iostream>


void check_default_ctor(){

	Penna::Genome test = Penna::Genome();
	Penna::Genome test2;

	for(std::size_t k=0; k< Penna::Genome::genome_size; ++k){
		if(test.get_genome()[k]!=false or test2.get_genome()[k]!=false){
			throw std::runtime_error("Genome::Genome (Ctor) not working!");
		}
	}
	return;
}

int main(){
	Penna::Genome::set_mutation_rate(2);
	try{
		check_default_ctor();
	}
	catch(std::exception & e){
		std::cout << e.what() << "\n";
	}
}