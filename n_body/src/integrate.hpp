#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <string>
#include <eigen3/Eigen/Dense>
#include <utility>
#include <iostream>
#include <cmath>
#include <fstream>

using mass_t = double const;
using name_t = std::string const;
using vector_t = Eigen::Vector3d;
using phase_t = Eigen::VectorXd;
using rhs_t = phase_t (*)(double const &, phase_t const &);
using step_t = phase_t (*)(rhs_t , phase_t const &, double, double);

phase_t explicit_euler_step(rhs_t, phase_t const &, double, double);
Eigen::MatrixXd integrate(step_t, rhs_t, phase_t const &, double const &, unsigned const &);
Eigen::MatrixXd explicit_euler(rhs_t, phase_t const &, double const &, unsigned const &);

Eigen::MatrixXd explicit_midpoint(rhs_t, phase_t const &, double const &, unsigned const &);

#endif // INTEGRATE_HPP