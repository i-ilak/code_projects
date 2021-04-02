#ifndef N_BODY_HPP
#define N_BODY_HPP

#include <vector>
#include "stelar_object.hpp"
#include "integrate.hpp"

using mass_t = double;
using mass_container_t = std::vector<mass_t>;
using vector_t = Eigen::Vector3d;
using phase_t = Eigen::VectorXd;    // type used for elements of "phase space"
using planet_container_t = std::vector<StellarObject>;


/*
 * Note that the following solution is based entirely on the approach
 * of pushing each ODE into the form dz=f(t,z), where z is an element
 * of phase-space, i.e. has the form:
 * z = [r_(1,x), r_(1,y), r_(1,z), r_(2,x), .... r_(N,z), 
 *                                  v_(1,x), v_(1,y), ..., v_(N, z)].
 * Here, r_1 is the position vector of object 1 and v_1 its velocity
 * vector, etc.
 */

/*
 * We need to calculate the norm of a vector at some point, i.e.
 *      ¦v¦=Sqrt(v_1^2+v_2^2+...),
 * Since the input is given in terms of Eigen::VectorXd, we need to
 * write a separate function for this. 
 */
// The outer functor we want to apply later.
double square(double x);
// The actual function, using the previously defined functor. 
double norm(vector_t const& v);


/*
 * Function that calculates the force of (N-1) objects onto one object.
 * PRE:
 *  - qk:           position of the object for which we want to calculate the force
 *  - masses:       masses of all the objects, where masses[k] is the mass of of the
 *                  object at qk. (see last parameter for k)
 *  - positions:    positions of all the objects. Note that positions[k]==qk.
 *  - N:            number of objects
 *  - k:            which object from masses are we considering
 * 
 * POST:
 *  - returns the force acting on object k. Note that the force is a vector, 
 *    i.e. direction and amplitude.
 */
vector_t force_on_object_k(vector_t const & qk, mass_container_t const& masses, 
                           phase_t const & positions, short const & N, short const & k);

/*
 * The rhs of the n-body problem involves only one complicated part,
 * i.e. the sum over all objects, calculating the forces that they
 * exert on each other. The following function does that part.
 * 
 * PRE:
 *  - z0:       initial conditions
 *  - masses:   masses of the involved bodies
 *  - G:        Gravitational constant
 * 
 * Note that G appears here because one needs to fix the unit system.
 * 
 * POST:
 *  - Returns ddz
 */
phase_t nbody_prod(phase_t const & z0, mass_container_t const& masses, double const & G);


/*
 * Function that puts everything together to create rhs of the n-body problem.
 */
phase_t nbody_rhs(phase_t const & z, mass_container_t const& masses, double const G);

struct nbody_rhs_helper{
    const mass_container_t masses_;
    const double grav_const_;
    nbody_rhs_helper(mass_container_t const & masses, double const G) : masses_(masses), grav_const_(G) {}
    phase_t operator()(double const & t, phase_t const & z){
        return nbody_rhs(z, masses_, grav_const_); 
    }
};

/*
 * PRE:
 *  - planets:      vector with stelar_obj instances, which we want to simulate.
 *  - T:            Time span we want to simulate
 *  - bins:         number of steps we want to take in simulation. The higher this
 *                  number the more accurate the simulation, but also increased
 *                  simulation time.
 *  - gravitational_constant:
 *                  Self-explanatory. We included this, together with the T because
 *                  we do not know in which unit system the user works. (usually
 *                  determined by the way how he chooses to represent the masses of 
 *                  planets.)
 *  - method:       is either 2 or 3. 
 *                   - 1: stands for Explicit Euler (default)
 *                   - 2: stands for Explicit midpoint
 *                   - 3: stands for Velocity Verlet 
 * POST:
 *  - Eigen::MatrixXd, where row k contains the phase space coordinates at time k*T/N,
 *    of all the planets. Elements of 3*u to 3*u+6 are the coordinates of u-th object. 
 */
Eigen::MatrixXd n_body_solver(planet_container_t const & planets, 
                              double const & T, int const & bins,
                              double const & gravitational_constant,
                              int const & method);

#endif //N_BODY_HPP