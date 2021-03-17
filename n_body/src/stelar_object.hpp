#ifndef STELLAR_OBJECT_HPP
#define STELLAR_OBJECT_HPP

#include "integrate.hpp"

// Class representing a stellar object
/*
 * The masses are taken relative to the Sun, 
 * distances are in A.U. (astronomical  units),
 * time in Earth days and 
 * the gravitational constant is G = 2.95912208286e-4.
 */
class StellarObject{
    public:
        explicit StellarObject( vector_t pos,vector_t speed,
                                mass_t const& mass, name_t const& name) 
                                : mass_(mass), name_(name),
                                  position_(std::move(pos)),
                                  velocity_(std::move(speed)) {};
        vector_t get_position() const {return position_;}
        vector_t get_velocity() const {return velocity_;}
        double get_mass() const {return mass_;}
    private:
        vector_t position_;
        vector_t velocity_;
        name_t name_;
        mass_t mass_;
};

#endif //STELLAR_OBJECT_HPP