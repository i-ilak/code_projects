#include "genome.hpp"
#include <stdexcept>
#include <iostream>


void check_copy_ctor(){

	Penna::Genome parent;
	Penna::Genome copy_parent(parent);

	for(std::size_t k = 0; k<Penna::Genome::genome_size; ++k){
		if(parent.get_genome()[k] != copy_parent.get_genome()[k]){
			std::runtime_error("Copy-Ctor is not working!");
		}
	}

	return;
}

int main(){
	Penna::Genome::set_mutation_rate(2);
	try{
		check_copy_ctor();
	}
	catch(std::exception & e){
		std::cout << e.what() << "\n";
	}
}