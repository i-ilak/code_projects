#include "n_body_solver.hpp"
#include <omp.h>

using mass_t = double;
using mass_container_t = std::vector<mass_t>;
using vector_t = Eigen::Vector3d;
using phase_t = Eigen::VectorXd;    // type used for elements of "phase space"
using planet_container_type = std::vector<StellarObject>;


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
double square(double x){  
    return std::pow(x,2);
}
// The actual function, using the previously defined functor. 
double norm(vector_t const& v){
    return std::sqrt(v.unaryExpr(&square).sum());
}


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
                           phase_t const & positions, short const & N, short const & k){
    vector_t force;
    force.fill(0);
    #pragma omp parallel for
    for(std::size_t i=0; i < N; ++i){
        if(i!=k){
            vector_t diff = positions.segment(3*i,3)-qk;
            double normq = std::pow(norm(diff),3);
            force += masses[i]/normq * diff;
        }
    }
    return force;
};

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
phase_t nbody_prod(phase_t const & z0,
                   mass_container_t const& masses, 
                   double const & G){
    short const N = masses.size();
    phase_t positions(3*N);                     //space for positions of the particles
    positions << z0.head(3*N);                  //fill it with the relevant elements from z
    phase_t ddz(3*N);                           //space for the rhs of the ODE
    ddz.fill(0);
    #pragma omp parallel for                                
    for(std::size_t k = 0; k<N; ++k){
        ddz.segment(3*k,3) = G * force_on_object_k(positions.segment(3*k,3), masses, positions, N, k);
    }
    return ddz;
}


/*
 * Function that puts everything together to create rhs of the n-body problem.
 */
phase_t nbody_rhs(phase_t const & z, mass_container_t const& masses,
                  double const G){
    int const N = masses.size();
    phase_t rhs(6*N);                               // space for the rhs
    rhs.head(3*N) = z.tail(3*N);                    // fill velocities in first 3*N elements of rhs
    rhs.tail(3*N) = nbody_prod(z, masses, G);       // fill last 3*N elements with nbody_prod
    return rhs;
}

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
Eigen::MatrixXd n_body_solver(planet_container_type const & planets, 
                              double const & T, int const & bins,
                              double const & gravitational_constant,
                              int const & method){
    mass_container_t masses;
    int const N = planets.size();
    phase_t z0(6*N);
    #pragma omp parallel for 
    for(std::size_t k=0; k < N; ++k) {
        z0.segment(3*k,3) = planets[k].get_position();
    }
    #pragma omp parallel for 
    for(std::size_t k=N; k < 2*N; ++k) {
        z0.segment(3*k,3) = planets[k%N].get_velocity();
    }
    for(auto & planet : planets){
        masses.push_back(planet.get_mass());
    }
    nbody_rhs_helper rhs = nbody_rhs_helper(masses, gravitational_constant);
    switch (method)
    {
    default:
        return explicit_euler(rhs, z0, T, bins);
        break;
    case 2:
        return explicit_midpoint(rhs, z0, T, bins);
        break;
    case 3:
        return velocity_verlet(rhs, z0, T, bins);
        break;
    }
    
}
