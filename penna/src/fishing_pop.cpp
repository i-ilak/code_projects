#include "fishing_pop.hpp"

// Returns true if animal gets fished.
struct FishingPredicate{
	using age_t = unsigned int;
	double f_rate;
	age_t f_age;
	FishingPredicate( double const & rate, 
					  age_t const & age ) : f_rate(rate), f_age(age) {};
	bool operator()(Penna::Animal const & animal){
		return (animal.age() > f_age) and (drand48() < f_rate);
	}
};

namespace Penna {

void FishingPopulation::change_fishing(	double const & f_rate, 
										age_t const & f_age){
	fishing_rate_ = f_rate;
	fishing_age_  = f_age;
}

void FishingPopulation::step(){
	Population::step();
	pop_.remove_if( FishingPredicate(fishing_rate_, fishing_age_) );
}

};