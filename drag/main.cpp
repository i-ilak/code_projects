#include "../n_body/src/integrate.hpp"
#include <ostream>
#include <math.h>   //Pi


int main(){
    /* We fix the coordinate system in such a way that the point where
     * The Rock starts jumping is at (x_star,y_star). We call the point where he wants
     * to land (0,0). Note that since he's jumping down, y_star>0.
     * 
     * First define all the necessary parameters for the simulation.
     */
    double const length_ref = 1.95;                     //meters (height of The Rock)
    double const mass = 118;                            //kg     (mass of The Rock)

    // Some coefficients necessary for the physics simulation
    double const area = 1.9/2.;
    double const cw = 0.78;                              // drag coefficient
    double const rho = 1.293;                            // density of air
    double const lam = cw * area * rho /(2.0 * mass);    // Coefficient of drag force

    // Target 
    double const y_star = 4.0*length_ref;
    double const x_star = 8.0*length_ref;

    // The Rock needs to choose his velocity or his angle appropriately to hit
    // the target.
}