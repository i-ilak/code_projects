#include "integrate.hpp"

/*
 For both methods we divide the computation into steps. Each step gets perform N times. 
 The following function implement this step. Later, these will be called and summed up.
 */

/*
  PRE:
    - rhs_t rhs:        f(t,z) for a given ODE.
    - phase_t const z0: initial condition in phase space coordinates. 
    - double t0:        time at which we start the step
    - double dt:        time step, usually defined as time interval of the simulation
                        divided by number of total steps taken in the simulation.
  POST:
    - phase_t z+dz:     phase space type coordinates
 */
phase_t explicit_euler_step(rhs_t rhs, phase_t const & z0, double t0, double dt){
    return z0 + dt * rhs(t0,z0);
}
phase_t implicit_euler_step(rhs_t rhs, phase_t const & z0, double t0, double dt) {
    auto F = [&](phase_t const& z){ return z0 - z + dt * rhs(t0, z);};
    phase_t initial_guess = z0 + dt * rhs(t0, z0);
    return phase_t(0);// scipy.optimize.fsolve(F, initial_guess);
}


/*
 Explicit Euler method.     For more details see: Wikipedia https://en.wikipedia.org/wiki/Euler_method
 Explicit Midpoint method.  For more detials see: Wikipedia https://en.wikipedia.org/wiki/Midpoint_method
 PRE:
    rhs_t f:            characterizes the ODE. 
    phase_t const z0:   initial conditions of the problem in form of phase-space
                        coordinates.
    double const T:     time interval you want to simulate, i.e. has to be a
                        positive real number.
    unsigned const N:   number of individual steps you want to take in simulation.
                        The bigger N is the more precise the result will be and  
                        the longer the computation will take. 
 POST:
    Eigen::MmatrixXd result:    row k contains the phase space coordinates at time k*T/N.
 */
Eigen::MatrixXd explicit_euler(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    Eigen::MatrixXd z(N+1, z0.size());
    z.fill(0);
    double const t0 = 0;
    double const dt = T/N;
    z.row(0) = z0;
    for(std::size_t i = 0; i<N; ++i) {
        z.row(i + 1) = explicit_euler_step(rhs, z.row(i), t0 + i * dt, dt);
    }
    phase_t t(N+1);
    for(std::size_t k = 0; k<t.size(); ++k){
        t(k) = (k+1) * dt;
    }
    z.conservativeResize(Eigen::NoChange, z.cols()+1);
    z.col(z.cols()-1) = t;
    return z;
}
Eigen::MatrixXd explicit_midpoint(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    Eigen::MatrixXd z(N+1, z0.size());
    z.fill(0);
    double const t0 = 0;
    double const dt =T/N;
    z.row(0) = z0;
    phase_t t(N+1);
    for(std::size_t k = 0; k<t.size(); ++k){
        t(k) = (k+1) * dt;
    }
    for(std::size_t i = 0; i<N; ++i) {
        phase_t tmp = z.row(i);
        tmp += (.5 * dt) * rhs(t(i), z.row(i));
        phase_t tmp2 = dt * rhs(t(i) + .5*dt, tmp);
        z.row(i + 1) = z.row(i);
        z.row(i + 1) += tmp2;
    }

    z.conservativeResize(Eigen::NoChange, z.cols()+1);
    z.col(z.cols()-1) = t;
    return z;
}
