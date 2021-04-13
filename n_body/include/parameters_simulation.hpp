#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include "stelar_object.hpp"

using mass_t = double;
using vector_t = Eigen::Vector3d;
using planet_container_type = std::vector<StellarObject>;

/*
 * We simulate a part of the Solar system. The initial data is given down 
 * below and corresponds to September 5th, 1994 at 00:00h.
 * The masses are taken relative to the Sun, distances are in A.U. (astronomical  units),
 * time in Earth days and the gravitational constant is G = 2.95912208286e-4.
 */
planet_container_type solar_system(){
    mass_t mSun = 1.00000597682;
    vector_t qSun = vector_t(0,0,0);
    vector_t vSun = vector_t(0,0,0);

    mass_t mj = 0.00095486104043;
    vector_t qj = vector_t(-3.5023653, -3.8169847, -1.5507963);
    vector_t vj = vector_t(0.00565429, -0.00412490, -0.00190589);

    mass_t ms = 0.000285583733151;
    vector_t qs = vector_t(9.0755314, -3.0458353, -1.6483708);
    vector_t vs = vector_t(0.00168318, 0.00483525, 0.00192462);

    mass_t mu = 0.0000437273164546;
    vector_t qu = vector_t (8.3101420, -16.2901086, -7.2521278);
    vector_t vu = vector_t (0.00354178, 0.00137102, 0.00055029);

    mass_t mn = 0.0000517759138449;
    vector_t qn = vector_t (11.4707666, -25.7294829, -10.8169456);
    vector_t vn = vector_t (0.00288930, 0.00114527, 0.00039677);

    mass_t mp = 1.0 / ( 1.3e8 );
    vector_t qp = vector_t (-15.5387357, -25.2225594, -3.1902382);
    vector_t vp = vector_t (0.00276725, -0.00170702, -0.00136504);

    StellarObject Sun(qSun  , vSun,mSun, "Sun");
    StellarObject Jupiter(qj, vj,    mj, "Jupiter");
    StellarObject Saturn( qs, vs,    ms, "Saturn");
    StellarObject Uranus( qu, vu,    mu, "Uranus");
    StellarObject Neptun( qn, vn,    mn, "Neptun");
    StellarObject Pluto(  qp, vp,    mp, "Pluto");

    planet_container_type planets = {Sun, Jupiter, Saturn, Uranus, Neptun, Pluto};
    return planets;
}

#endif 