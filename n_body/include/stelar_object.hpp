#ifndef STELLAR_OBJECT_HPP
#define STELLAR_OBJECT_HPP

#include <string>
#include <eigen3/Eigen/Dense>

using mass_t = double;
using name_t = std::string const;
using vector_t = Eigen::Vector3d;

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
        name_t get_name() const {return name_;}
    private:
        vector_t position_;
        vector_t velocity_;
        name_t name_;
        mass_t mass_;
};

#endif //STELLAR_OBJECT_HPP