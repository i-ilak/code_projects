#include "animal.hpp"
#include <cassert>
#include <random>

namespace Penna {


age_t Animal::reporoduction_age_;
age_t Animal::bad_threshold_;
double Animal::prob_to_get_pregnant_;

bool Animal::is_pregnant() const {
    return pregnant_;
}

age_t Animal::age() const {
    return age_;
}

void Animal::set_threshold(age_t num){
	bad_threshold_ = num;
}

age_t Animal::get_threshold(){
	return bad_threshold_;
}

void Animal::set_reproduction_age(age_t num){
    reporoduction_age_ = num;
}

age_t Animal::get_reproduction_age(){
    return reporoduction_age_;
}

void Animal::set_prob_to_get_pregnant(age_t num) {
    prob_to_get_pregnant_ = num;
}

double Animal::get_prob_to_get_pregnant(){
    return prob_to_get_pregnant_;
}

bool Animal::is_dead() const{
    return (age_ > Genome::genome_size) 
            or (genome_.count_bad(age_) > bad_threshold_);
}

void Animal::advance_by_one_year() {
    assert(! is_dead());
    age_+=1;
    if(age_ > reporoduction_age_ and !pregnant_){
        if(drand48() < prob_to_get_pregnant_){
            pregnant_ = true;
        }
    }
    return;
}


Animal Animal::give_birth() const {
            assert( pregnant_ );
            Genome baby_genes = genome_;
            baby_genes.mutate();
            Animal baby = Animal(baby_genes);
            return baby;
        }

}; // Penna namespace
