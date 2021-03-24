#include "genome.hpp"
#include <cassert>
#include <vector>
#include <algorithm>

namespace Penna{

age_t Genome::mutation_rate_;

void Genome::set_mutation_rate(age_t num){
	assert(num <= Genome::genome_size);
	mutation_rate_ = num;
}

age_t Genome::get_mutation_rate(){
	return mutation_rate_;
}

std::bitset<Genome::genome_size> Genome::get_genome() const{
	return genes_;
}

age_t Genome::count_bad(age_t num) const {
	assert(num <= Genome::genome_size);
	age_t cnt = 0;
	for(std::size_t k = 0; k < num; ++k){
		if(genes_[k] == true){
			cnt +=1 ;
		}
	}
	return cnt;
}

void Genome::mutate(){
	std::vector<age_t> indices;
	for(std::size_t k = 0; k < Genome::genome_size; ++k){
		indices.push_back(k);
	}

	std::random_shuffle(indices.begin(), indices.end());

	for(std::size_t k=0; k < mutation_rate_; ++k){
		genes_[indices[k]] = not genes_[indices[k]];
	}

	return;
}

}; // Penna namespace