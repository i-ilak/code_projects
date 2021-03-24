#include "population.hpp"
#include <cassert>
#include <functional>
#include <algorithm>

// Returns true if aniaml should die
struct DeathPrediacte{
	double n_over_nmax;
	explicit DeathPrediacte(double const & num) : n_over_nmax(num) {};
	bool operator()(Penna::Animal const & animal){
        return n_over_nmax > 1. or drand48() < n_over_nmax;
	}	
};

namespace Penna {

Population::Population(age_t const & how_big, age_t const & nmax){
	    for(std::size_t k = 0; k < how_big; ++k){
		    pop_.push_back(Animal());
	    }
	    max_size_=nmax;
}

age_t Population::size() const{
	return pop_.size();
}

void Population::step() {
	std::for_each( pop_.begin(),
				   pop_.end(),
				   std::mem_fn( &Animal::advance_by_one_year ));

	pop_.remove_if( std::mem_fn( &Animal::is_dead ) );
	pop_.remove_if( DeathPrediacte( pop_.size()/double(max_size_) ) );


	cont_t parents;
	std::copy_if( pop_.begin(),
				  pop_.end(),
				  std::back_inserter(parents),
				  std::mem_fn(&Animal::is_pregnant));

	std::transform( parents.begin(),
					parents.end(),
					std::back_inserter(pop_),
					std::mem_fn(&Animal::give_birth));
}

};
