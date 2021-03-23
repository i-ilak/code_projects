#include "../n_body/src/integrate.hpp"
#include <ostream>
#include <math.h>   //Pi
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <tuple>

// Small helper function to get Python's np.linspace
template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

namespace constants{
    //First define all the necessary parameters for the simulation.
    double const length_ref = 1.95;                     //meters (height of The Rock)
    double const mass = 118;                            //kg     (mass of The Rock)

    // Some coefficients necessary for the physics simulation
    double const area = 1.9/2.;
    double const cw = 0.78;                             // drag coefficient
    double const rho = 1.293;                           // density of air
    double const lam = cw * area * rho /(2.0 * mass);   // Coefficient of drag force
    double const grav_const = 9.81;                     // m/s^2
} // constants


std::pair<double,double> calc_tstar_ref_and_x_ref_at_tstar_ref(double const & v0,
                                                               double const & theta, 
                                                               double const & ystar){
    double const vx = v0*std::cos(theta);
    double const vy = v0*std::sin(theta);
    double const t_star_ref = (vy + std::sqrt(vy*vy + 2.0*constants::grav_const * ystar))/constants::grav_const;
    double const x_star_ref  = vx*t_star_ref;
    return std::pair<double,double>(t_star_ref, x_star_ref);
}

phase_t rhs_with_drag(double const & t, phase_t const & z){
    phase_t rhs(4);                                 // space for the rhs
    rhs.head(2) = z.tail(2);                        // fill velocities in first 2 elements of rhs
    phase_t tmp;
    double const v_norm = std::sqrt(z(2)*z(2)+z(3)*z(3));
    tmp << -constants::lam*z(2)*v_norm, -constants::lam*z(3)*v_norm-constants::grav_const;
    rhs.tail(2) = tmp;
    return rhs;
}

std::tuple<double,double,double,int,std::pair<double,double>> try_tstar_drag(
                                           double const & v0, 
                                           double const & theta, 
                                           double const & y_star, 
                                           std::vector<double> const & tstar_search_grid){
    auto t_grid = tstar_search_grid;
    t_grid.insert(t_grid.begin(), 0.0);
    phase_t z0(4);
    z0(0)=0.0;
    z0(1)=y_star;
    z0(2)=v0*std::cos(theta);
    z0(3)=v0*std::sin(theta);
    std::cout << "got to here" <<"\n";
    auto result = explicit_midpoint(rhs_with_drag, z0, tstar_search_grid.back(), tstar_search_grid.size());
    phase_t y_drag = result.col(1);
    phase_t x_drag = result.col(0);
    int ierr=false;
    unsigned t_star_index = 0;
    std::pair<double,double> t_star_bracket;
    if(y_drag(0)<0.0){
    	ierr = -1;
    	t_star_index = 0;
    	t_star_bracket.first = 0;
    	t_star_bracket.second = 0;
    }
    if(y_drag(y_drag.size()-1)>0.0){
    	ierr = 1;
    	t_star_index = y_drag.size()-1;
    	t_star_bracket.first = 0;
    	t_star_bracket.second = 0;
    }
    else{
    	ierr = 0;
    	for(std::size_t i=0; i<y_drag.size()-1;i++){
    		if(y_drag[i+1]*y_drag[i]<=0.0){
    			t_star_bracket.first = t_grid[i];
    			t_star_bracket.second = t_grid[i+1];
    			if (std::abs(y_drag[i+1])<std::abs(y_drag[i])){
    				t_star_index = i+1;
    			}
    			else{
    				t_star_index = i;
    			}
    		}
    	}
    }
    
    return std::make_tuple(t_grid[t_star_index], x_drag[t_star_index], y_drag[t_star_index], ierr, t_star_bracket);
}

std::tuple<double,double,double,bool> calc_x_star_with_drag_at_tstar(
                                           double const & v0, 
                                           double const & theta, 
                                           double const & y_star, 
                                           double t_star_est,
                                           double const eps_y=1.0e-3,
                                           int const & number_of_t_searches=20,
                                           double const & rel_range_lower=0.8,
                                           double const & rel_range_upper=1.2,
                                           int const & num_try=5){
    bool success = false;
    double t_star_est_lower=t_star_est*rel_range_lower;
    double t_star_est_upper=t_star_est*rel_range_upper;
    // Space to save results
    int ierr=false;
    unsigned t_star_index = 0;
    std::pair<double,double> t_star_bracket;
    double y_drag;
    double x_drag;
    for(std::size_t try_=0; try_<num_try; try_++){
        auto t_star_search_grid=linspace(t_star_est_lower,t_star_est_upper,number_of_t_searches);
        auto tmp = try_tstar_drag(v0,theta,y_star,t_star_search_grid);
        t_star_est = std::get<0>(tmp);
        x_drag = std::get<1>(tmp);
        y_drag = std::get<2>(tmp);
        ierr = std::get<3>(tmp);
        t_star_bracket = std::get<4>(tmp);

        if(ierr==-1){
        	t_star_est_upper=t_star_est_lower;
        	t_star_est_lower=t_star_est_lower*rel_range_lower;
        }
        if(ierr==1){
        	t_star_est_lower=t_star_est_upper;
        	t_star_est_upper=t_star_est_upper*rel_range_upper;
        }
        else{
        	if(std::abs(y_drag)<eps_y){
        		success=true;
        		break;
        	}
        	else{
        		t_star_est_lower=t_star_bracket.first;
        		t_star_est_upper=t_star_bracket.second;
        	}
        }
    } 
    return std::make_tuple(t_star_est, x_drag, y_drag, success);
}

std::tuple<double,double,bool,int> find_v0(double const & x_star, 
                                           double v0_est, 
                                           double const & theta, 
                                           double const & y_star, 
                                           double t_star_est,
                                           double const & eps_x=1.0e-3,
                                           double const & num_try=6){
    bool success = false;
    std::vector<double> v0_history;
    std::vector<double> x_star_history;
    int number_of_tries = 0;
    double x_drag;
    double y_drag;
    for(std::size_t try_=0; try_<num_try; try_++){
        auto tmp = calc_x_star_with_drag_at_tstar(v0_est, theta, y_star, t_star_est);
        t_star_est = std::get<0>(tmp);
        x_drag = std::get<1>(tmp);
        y_drag = std::get<2>(tmp);
        success = std::get<3>(tmp);
        
        v0_history.push_back(v0_est);
        x_star_history.push_back(x_drag);
        if(success==false){                             //<-------------- probably wrong
            break;
        }
        if(std::abs(x_drag-x_star)<eps_x){
            success=true;
            number_of_tries=try_;
        }
        else{
            if(v0_history.size()<2){
                double dx = x_star - x_drag;
                double dv0 = dx/(t_star_est*std::cos(theta));
                v0_est += dv0;
            }
            else{
                int const end_v = v0_history.size();
                int const end_x = x_star_history.size();
                v0_est  = v0_history[end_v-2] 
                        + (v0_history[end_v-1]-v0_history[end_v-2])
                        *(x_star-x_star_history[end_x-2])
                        /(x_star_history[end_x-1]-x_star_history[end_x-2]);
            }
        }
    }
    return std::make_tuple(0,0,0,0);
}

int main(){
    /* We fix the coordinate system in such a way that the point where
     * The Rock starts jumping is at (x_star,y_star). We call the point where he wants
     * to land (0,0). Note that since he's jumping down, y_star>0.
     */

    // Target 
    double const y_star = 4.0*constants::length_ref;
    double const x_star = 8.0*constants::length_ref;

    /* The Rock needs to choose his velocity and his angle appropriately to hit
     * the target. We solve for the possible combinations as follows:
     *  1.) We call the time when The Rock would hit the ground t_star.
     *  2.) We calculate a reference time t_star_ref and the drop point x_star_ref
     *      neglecting air drag.
     *  3.) Using our calculated reference we make an estimate of v0 and t_star in
     *      the dragged system.
     *  4.) We compute the final v0 for each theta under air drag. 
     * 
     * Note that for the simulation of air drag we use Newtons's Drag (for more 
     * infos see Wikipedia on the topic), which should be appropriate for a human
     * flying through air. 
     */

    // We fix the number of sample points
    unsigned const theta_n = 51;
    double const theta_min = -0.1*M_PI;             // rad
    double const theta_max = 0.4*M_PI;              // rad
    unsigned const v0_n = 100;
    double const v0_min = 6.0;                      // m/s
    double const v0_max = 18.0;                     // m/s


    // Create theta and v0 grid. These will be all the possible combinations of 
    // angle and velocities that The Rock can have at the start.
    std::vector<double> v0_grid = linspace(v0_min, v0_max, v0_n);
    std::vector<double> theta_grid = linspace(theta_min, theta_max, theta_n);

    // Grid on which we evaluate calc_tstar_ref_and_x_ref_at_tstar_ref
    Eigen::MatrixXd t_star_ref(theta_n, v0_n);
    Eigen::MatrixXd x_star_ref(theta_n, v0_n);
    for(std::size_t i=0; i<theta_n; i++){
        for(std::size_t j=0; j<v0_n; j++){
            std::pair<double, double> const tmp = calc_tstar_ref_and_x_ref_at_tstar_ref(v0_grid[i], theta_grid[j], y_star);
            t_star_ref(i,j)=tmp.first;
            x_star_ref(i,j)=tmp.second;
        }
    }

    // We now make the estimate using linear interpolation
    std::vector<double> t_star_est(theta_n);
    std::vector<double> v0_est(theta_n);
    for(auto & t : t_star_est){
        t=-1.;                      // Null value
    }
    for(auto & v : v0_est){
        v=-1.;                      // Null value
    }

    for(std::size_t i=0; i<theta_n; i++){
        for(std::size_t j=0; j<v0_n-1; j++){
            if((x_star_ref(i,j+1)-x_star)*(x_star_ref(i,j)-x_star)<=0.0){
                t_star_est[i] = t_star_ref(i,j);
                v0_est[i] = v0_grid[j] 
                          + (v0_grid[j+1]-v0_grid[j])*(x_star-x_star_ref(i,j))
                          /(x_star_ref(i,j+1)-x_star_ref(i,j));
                break;
            }
        }
    }

    // We compute now v0 for each theta under air drag
    std::vector<double> theta_sols;
    std::vector<double> t_star_sols;
    std::vector<double> v0_sols;
    for(std::size_t i=0; i<theta_n; i++){
        if(v0_est[i]>0.0){
            auto sol = find_v0(x_star, v0_est[i], theta_grid[i], y_star, t_star_est[i]);
        }
    }
}