#include <string>
#include <eigen3/Eigen/Dense>
#include <utility>
#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include "integrate.hpp"
#include <vector>


double square(double x){  // the functor we want to apply
    return std::pow(x,2);
}

double norm(vector_t const& v){
    return std::sqrt(v.unaryExpr(&square).sum());
}

// Class representing a stellar object
class StellarObject{
    public:
        explicit StellarObject(  vector_t pos = vector_t(0,0,0),       //Position of the sun relative to itself
                                 vector_t speed = vector_t(0,0,0),     //Velocity of the sun relative to itself
                                 mass_t const& mass = 1.00000597682,         //Mass of the Sun [kg]
                                 name_t const& name = "Sun") : mass_(mass),
                                                               name_(name),
                                                               position_(std::move(pos)),
                                                               velocity_(std::move(speed)) {};
        vector_t get_position() const {
            return position_;
        }
        vector_t get_velocity() const {
            return velocity_;
        }
        double get_mass() const {
            return mass_;
        }
    private:
        vector_t position_;
        vector_t velocity_;
        name_t name_;
        mass_t mass_;
};


phase_t nbody_prod(double const & t, phase_t const & z,
                   std::vector<double> const& masses, double const & G){
    int const N = masses.size();
    phase_t q(3*N);                     //space for positions of the particles
    q << z.head(3*N);                   //fill it with the relevant elements from z
    phase_t ddx(3*N);                   //space for the rhs of the ODE
    ddx.fill(0);                       //make sure that set to zero everywhere
    for(std::size_t k = 0; k<N; ++k){
        vector_t tmp;                       //tmp object to store ddx
        tmp.fill(0);
        vector_t qk = q.segment(3*k,3); //get position of k-th object
        for(std::size_t i=0; i < N; ++i){
            if(i!=k){                       //calculate acceleration of k-th object
                vector_t qi = q.segment(3*i,3);
                vector_t diff = qi-qk;
                double normq = std::pow(norm(diff),3);
                tmp += masses[i]/normq * diff;
            }
        }
        ddx.segment(3*k,3) = G * tmp;
    }
    return ddx;
}

// z = [r_(1,x), r_(1,y), r_(1,z), r_(2,x), .... r_(N,z), v_(1,x), v_(1,y), ..., v_(N, z)]
phase_t nbody_rhs(double const & t, phase_t const & z, std::vector<double> const& masses,
                  double const G){
    int const N = masses.size();
    phase_t rhs(6*N);                     //space for the rhs
    rhs.head(3*N) = z.tail(3*N);        // fill velocities in first 3*N elements of rhs
    rhs.tail(3*N) = nbody_prod(t, z, masses, G); // fill last 3*N elements with nbody_prod
    return rhs;
}

static double G;
static std::vector<double> masses;

Eigen::MatrixXd n_body_solver(std::vector<StellarObject> const & planets){
    int const N = planets.size();
    phase_t z0(6*N);
    for(std::size_t k=0; k < N; ++k) {
        z0.segment(3*k,3) = planets[k].get_position();
    }
    for(std::size_t k=N; k < 2*N; ++k) {
        z0.segment(3*k,3) = planets[k%N].get_velocity();
    }
    for(auto & planet : planets){
        masses.push_back(planet.get_mass());
    }
    G = 2.95912208286e-4;
    auto reduced_rhs = [](double const & t, phase_t const & z0) {return nbody_rhs(t, z0, masses, G);};
    return explicit_midpoint(reduced_rhs, z0, 60000, 40000);
}

/*
class n_body_solver{
    public:
        n_body_solver(std::vector<StellarObject> const & planets, 
                      double const & G=2.95912208286e-4) : N_(planets.size()) {
            G_= G;
            for(std::size_t k=0; k < N_; ++k) {
                z0_.segment(3*k,3) = planets[k].get_position();
            }
            for(std::size_t k=N_; k < 2*N_; ++k) {
                z0_.segment(3*k,3) = planets[k%N_].get_velocity();
            }
            for(auto & planet : planets){
                masses_.push_back(planet.get_mass());
            }
        }

        phase_t nbody_prod(double const & t, phase_t const & z,
                           std::vector<double> const& masses, double const & G){
            int const N = masses.size();
            phase_t q(3*N);                     //space for positions of the particles
            q << z.head(3*N);                   //fill it with the relevant elements from z
            phase_t ddx(3*N);                   //space for the rhs of the ODE
            ddx.fill(0);                       //make sure that set to zero everywhere
            for(std::size_t k = 0; k<N; ++k){
                vector_t tmp;                       //tmp object to store ddx
                tmp.fill(0);
                vector_t qk = q.segment(3*k,3); //get position of k-th object
                for(std::size_t i=0; i < N; ++i){
                    if(i!=k){                       //calculate acceleration of k-th object
                        vector_t qi = q.segment(3*i,3);
                        vector_t diff = qi-qk;
                        double normq = std::pow(norm(diff),3);
                        tmp += masses[i]/normq * diff;
                    }
                }
                ddx.segment(3*k,3) = G * tmp;
            }
            return ddx;
        }

        // z = [r_(1,x), r_(1,y), r_(1,z), r_(2,x), .... r_(N,z), v_(1,x), v_(1,y), ..., v_(N, z)]
        phase_t nbody_rhs(double const & t, phase_t const & z, std::vector<double> const& masses,
                          double const G){
            int const N = masses.size();
            phase_t rhs(6*N);                     //space for the rhs
            rhs.head(3*N) = z.tail(3*N);        // fill velocities in first 3*N elements of rhs
            rhs.tail(3*N) = nbody_prod(t, z, masses, G); // fill last 3*N elements with nbody_prod
            return rhs;
        }

        void solve(double const & time_intervall, unsigned int const & bins){
            auto reduced_rhs = [](double const & t, phase_t const & z0) {return nbody_rhs(t, z0, masses_, G_);};
            result_ = explicit_midpoint(reduced_rhs, z0_, time_intervall, bins);
        }

        Eigen::MatrixXd get_result() const {
            return result_;
        }

    private:
        phase_t z0_;
        unsigned short const N_;
        static double G_;
        static std::vector<double> masses_;
        Eigen::MatrixXd result_;
};
*/

int main(){
    /*
     * Simple 2-body problem to see if the code above works as expected.
     * Make sure to set G = 1.
    vector_t r0 = vector_t(2,0,0);
    vector_t v0 = vector_t(0, std::sqrt(500/2), 0);

    StellarObject Sun;
    StellarObject Earth(r0, v0, 1, "Earth");


    std::vector<StellarObject> planets = {Sun, Earth};
     */
    /*
     * Special case of the 3-body problem as an extended test.
     * Make sure to set G = 1.
    vector_t q1 = vector_t(0.97000436, -0.24308753, 0);
    vector_t v1 = vector_t(0.46620368, 0.43236573, 0);
    vector_t q2 = vector_t(-0.97000436, 0.24308753, 0);
    vector_t v2 = vector_t(0.46620368, 0.43236573, 0);
    vector_t q3 = vector_t(0, 0, 0);
    vector_t v3 = vector_t(-0.93240737, -0.86473146, 0);

    StellarObject Earth1(q1, v1, 1, "Earth");
    StellarObject Earth2(q2, v2, 1, "Earth");
    StellarObject Earth3(q3, v3, 1, "Earth");

    std::vector<StellarObject> planets = {Earth1, Earth2, Earth3};
    */

    /*
     * Simple simulation of some planets in the solar system.
     * Make sure to set G = 2.95912208286e-4,
     *                  msun = 1.00000597682.
     */
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

    mass_t mp = 7.692307692307693e-09;
    vector_t qp = vector_t (-15.5387357, -25.2225594, -3.1902382);
    vector_t vp = vector_t (0.00276725, -0.00170702, -0.00136504);

    StellarObject Sun;
    StellarObject Jupiter(qj,vj, mj, "Jupiter");
    StellarObject Saturn(qs,vs, ms, "Saturn");
    StellarObject Neptun(qn,vn, mn, "Neptun");
    StellarObject Pluto(qp,vp, mp, "Pluto");

    std::vector<StellarObject> planets = {Sun, Jupiter, Saturn, Neptun, Pluto};

    // Running the calculation and timing it
    auto start = std::chrono::high_resolution_clock::now();
    auto res = n_body_solver(planets);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "The relevant calculation took:\t" << diff.count() << " s" <<"\n";

    // Writing the data into a file for later plotting
    std::string path="data.txt";
    std::ofstream out(path);

    out << planets.size() << "\n";
    auto start2 = std::chrono::high_resolution_clock::now();
    out << res << "\n";
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2 = end2 - start2;
    std::cout << "Time to write to output-file:\t" << diff2.count() << " s" <<"\n";



    /*
    // Running the calculation and timing it
    
    n_body_solver test(planets);
    auto start = std::chrono::high_resolution_clock::now();
    test.solve(60000,40000);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "The relevant calculation took:\t" << diff.count() << " s" <<"\n";

    // Writing the data into a file for later plotting
    std::string path="data.txt";
    std::ofstream out(path);

    out << planets.size() << "\n";
    auto start2 = std::chrono::high_resolution_clock::now();
    out << test.get_result() << "\n";
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2 = end2 - start2;
    std::cout << "Time to write to output-file:\t" << diff2.count() << " s" <<"\n";
    */
}
