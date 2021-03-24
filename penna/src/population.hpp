#ifndef POPULATION_HPP
#define POPULATION_HPP

#include <list>
#include "animal.hpp"
#include "penna_vector.hpp"

namespace Penna {

using cont_t = std::list<Penna::Animal>;
//using cont_t = penna_vector<Penna::Animal>;

class Population {
public:
	Population(age_t const &, age_t const &);
		age_t size() const;
		void step();

	protected:
		cont_t pop_;
	private:
		age_t max_size_;
};

}; // Penna namespace

#endif // POPULATION_HPP
