#ifndef __ANIMAL_HPP__
#define __ANIMAL_HPP__

#include "genome.hpp"

namespace Penna {

class Animal {
    
    public:
        Animal() : pregnant_(0), age_(0), genome_(Genome()) {}

        Animal(Genome const & genes) : pregnant_(0), 
                                     age_(0), 
                                     genome_(genes) {}

        bool is_pregnant() const;
        age_t age() const;
        bool is_dead() const;

        static void set_threshold(age_t);
        static age_t get_threshold();

        static void set_reproduction_age(age_t);
        static age_t get_reproduction_age();

        static void set_prob_to_get_pregnant(age_t);
        static double get_prob_to_get_pregnant();

        void advance_by_one_year();

        Animal give_birth() const ;

    private:
        static age_t reporoduction_age_;
        static age_t bad_threshold_;
        static double prob_to_get_pregnant_;

        bool pregnant_;
        age_t age_;
        Genome genome_; // can be const if we change cont_t to std::list
};

}; // Penna namespace

#endif // __ANIMAL_HPP__
