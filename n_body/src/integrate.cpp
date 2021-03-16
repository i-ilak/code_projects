#include "integrate.hpp"

Eigen::MatrixXd integrate(step_t step, rhs_t rhs, phase_t const & z0,
                          double const & time_interval, unsigned const & bins){
    Eigen::MatrixXd z(bins+1, z0.size());
    z.fill(0);
    double const t0 = 0;
    double const dt = time_interval/bins;
    z.row(0) = z0;
    for(std::size_t i = 0; i<bins; ++i) {
        z.row(i + 1) = step(rhs, z.row(i), t0 + i * dt, dt);
    }
    phase_t t(bins+1);
    for(std::size_t k = 0; k<t.size(); ++k){
        t(k) = (k+1) * dt;
    }
    z.conservativeResize(Eigen::NoChange, z.cols()+1);
    z.col(z.cols()-1) = t;
    return z;
}

phase_t explicit_euler_step(rhs_t rhs, phase_t const & z0, double t0, double dt){
    return z0 + dt * rhs(t0,z0);
}
phase_t implicit_euler_step(rhs_t rhs, phase_t const & z0, double t0, double dt) {
    auto F = [&](phase_t const& z){ return z0 - z + dt * rhs(t0, z);};
    phase_t initial_guess = z0 + dt * rhs(t0, z0);
    return phase_t(0);// scipy.optimize.fsolve(F, initial_guess);
}



Eigen::MatrixXd explicit_euler(rhs_t rhs, phase_t const &z0, double const & T, unsigned const & N){
    return integrate(explicit_euler_step,rhs,z0,T,N);
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
