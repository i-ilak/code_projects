#include <string>
#include <utility>
#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include "integrate.hpp"
#include "stelar_object.hpp"
#include "n_body_solver.hpp"


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
     * We simulate a part of the Solar system.
     * The initial data is given down below and corresponds to 
     * September 5th, 1994 at 00:00h.
     * The masses are taken relative to the Sun, 
     * distances are in A.U. (astronomical  units),
     * time in Earth days and 
     * the gravitational constant is G = 2.95912208286e-4.
     */

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

    mass_t mp = 7.692307692307693e-09;
    vector_t qp = vector_t (-15.5387357, -25.2225594, -3.1902382);
    vector_t vp = vector_t (0.00276725, -0.00170702, -0.00136504);

    StellarObject Sun(qSun,vSun,mSun, "Sun");
    StellarObject Jupiter(qj,vj, mj, "Jupiter");
    StellarObject Saturn(qs,vs, ms, "Saturn");
    StellarObject Uranus(qu,vu, mu, "Uranus");
    StellarObject Neptun(qn,vn, mn, "Neptun");
    StellarObject Pluto(qp,vp, mp, "Pluto");

    std::vector<StellarObject> planets = {Sun, Jupiter, Saturn, Uranus, Neptun, Pluto};
    double const G = 2.95912208286e-4;
    int const N = 40000;
    double const T=60000;

    // Running the calculation and timing it
    auto start = std::chrono::high_resolution_clock::now();
    auto res = n_body_solver(planets, T, N, G);
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

}
