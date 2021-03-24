#ifndef FISHINGPOPULATION_HPP
#define FISHINGPOPULATION_HPP

#include "population.hpp"

namespace Penna {

class FishingPopulation : public Population {
    public:
        FishingPopulation( age_t  const & how_many,
                           age_t  const & nmax,
                           double const & f_rate,
                           age_t  const & f_age) 
                                : Population(how_many, nmax),
                                  fishing_rate_(f_rate),
                                  fishing_age_(f_age) {};
        void change_fishing(double const &, age_t const &);
        void step();

    private:
        double fishing_rate_;
        age_t fishing_age_;
};

}; // Penna namespace

#endif // POPULATION_HPP