#include "ode_int_sol.hpp"
#include <H5Cpp.h>

// The following code is pretty much verbatim copied from 
// https://github.com/headmyshoulder/odeint-v2/blob/master/examples/solar_system.cpp

/* Boost libs/numeric/odeint/examples/solar_system.cpp
 Copyright 2010-2012 Karsten Ahnert
 Copyright 2011 Mario Mulansky
 Solar system example for Hamiltonian stepper
 Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */

// We will wrap the code to make it callable in a similar way as the other ode-solvers.
// We also modify the streaming_observer slightly to be able to write the data to a 
// HDF5 document in the end.

typedef point< double , 3 > point_type;
typedef boost::array< point_type , n > container_type;
typedef boost::array< double , n > mass_type;
//]

// added container_type for planets
using planet_container_type = std::vector<StellarObject>;

struct solar_system_coor
{
    const mass_type &m_masses;

    solar_system_coor( const mass_type &masses ) : m_masses( masses ) { }

    void operator()( const container_type &p , container_type &dqdt ) const
    {
        for( size_t i=0 ; i<n ; ++i )
            dqdt[i] = p[i] / m_masses[i];
    }
};
//]


//[ momentum_function
struct solar_system_momentum
{
    const mass_type &m_masses;

    solar_system_momentum( const mass_type &masses ) : m_masses( masses ) { }

    void operator()( const container_type &q , container_type &dpdt ) const
    {
        const size_t n = q.size();
        for( size_t i=0 ; i<n ; ++i )
        {
            dpdt[i] = 0.0;
            for( size_t j=0 ; j<i ; ++j )
            {
                point_type diff = q[j] - q[i];
                double d = abs( diff );
                diff *= ( gravitational_constant * m_masses[i] * m_masses[j] / d / d / d );
                dpdt[i] += diff;
                dpdt[j] -= diff;

            }
        }
    }
};
//]







//[ some_helpers
point_type center_of_mass( const container_type &x , const mass_type &m )
{
    double overall_mass = 0.0;
    point_type mean( 0.0 );
    for( size_t i=0 ; i<x.size() ; ++i )
    {
        overall_mass += m[i];
        mean += m[i] * x[i];
    }
    if( !x.empty() ) mean /= overall_mass;
    return mean;
}

//[ streaming_observer
struct streaming_observer
{
    std::vector<std::vector<double>>& res_;

    streaming_observer( std::vector<std::vector<double>>& buffer ) : res_( buffer ) { }

    template< class State >
    void operator()( const State &x , double t ) const
    {   
        container_type &q = x.first;
        std::vector<double> row;
        for( size_t i=0 ; i<q.size() ; ++i ){
            for(size_t j=0; j<3; j++){
                row.push_back(q[i][j]);
            }
        }
        res_.push_back(row);
        row.push_back(t);
    }
};


void perform_ode_int_simulation(planet_container_type const & objects,
                        double const & T, int const & N){

    using namespace std;
    using namespace boost::numeric::odeint;

    mass_type masses;
    container_type q;
    container_type p;
    for(std::size_t i=0; i<objects.size();i++){
        masses[i] = objects[i].get_mass();
        for(std::size_t j=0; j<3; j++){
            q[i][j] = objects[i].get_position()[j];
            p[i][j] = objects[i].get_velocity()[j] * masses[i];
        }
    }

    point_type qmean = center_of_mass( q , masses );
    point_type pmean = center_of_mass( p , masses );
    for( size_t i=0 ; i<n ; ++i )
    {
        q[i] -= qmean ;
        p[i] -= pmean;
    }

    //[ integration_solar_system
    typedef symplectic_rkn_sb3a_mclachlan< container_type > stepper_type;
    const double dt = T/N;

    // Create space to save data
    std::vector<std::vector<double>> results;
    
    auto start = std::chrono::high_resolution_clock::now();
    integrate_const(
            stepper_type() ,
            make_pair( solar_system_coor( masses ) , solar_system_momentum( masses ) ) ,
            make_pair( boost::ref( q ) , boost::ref( p ) ) ,
            0.0 , T, dt , streaming_observer( results ) );
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "ODEInt calculation took:\t" << diff.count() << " s" <<"\n";

    // Write stuff to HDF5 file 
    const std::string file_name = "./data/ode_int.hdf5";
    const H5std_string  FILE_NAME( file_name );
    H5::H5File file( FILE_NAME, H5F_ACC_TRUNC );
    const int   NX = results.size();                    // dataset dimensions
    const int   NY = 3;
    const int   RANK = 2;
    auto start2 = std::chrono::high_resolution_clock::now();
    for(std::size_t i=0; i<objects.size(); i++){
        const H5std_string  DATASET_NAME( objects[i].get_name() );

        int k, j;
        double data_arr[NX][NY];            // buffer for data to write
        // Note that data_arr is on the stack, i.e. we don't need to clean it up at the end!
        for (j = 0; j < NX; j++)
        {
            for (k = 0; k < NY; k++)
                data_arr[j][k] = results[j][3*i+k];
        }
        
        hsize_t dimsf[2];                   // dataset dimensions
        dimsf[0] = NX;
        dimsf[1] = NY;
        H5::DataSpace dataspace( RANK, dimsf );
        /*
        * Define datatype for the data in the file.
        * We will store little endian DOUBLE numbers.
        */
        H5::FloatType datatype( H5::PredType::NATIVE_DOUBLE );
        datatype.setOrder( H5T_ORDER_LE );

        H5::DataSet dataset = file.createDataSet( DATASET_NAME, datatype, dataspace );

        dataset.write( data_arr, H5::PredType::NATIVE_DOUBLE );
    }
    // Need to write time separately
    const H5std_string  DATASET_NAME( "Time" );
    double time_arr[NX];
    for (int j = 0; j < NX; j++){
        time_arr[j] = results[j][results[0].size()-1];
    }
    hsize_t dimsf[1];                   // dataset dimensions
    dimsf[0] = NX;
    H5::DataSpace dataspace( 1 , dimsf );
    /*
    * Define datatype for the data in the file.
    * We will store little endian DOUBLE numbers.
    */
    H5::FloatType datatype( H5::PredType::NATIVE_DOUBLE );
    datatype.setOrder( H5T_ORDER_LE );

    H5::DataSet dataset = file.createDataSet( DATASET_NAME, datatype, dataspace );

    dataset.write( time_arr, H5::PredType::NATIVE_DOUBLE );
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2 = end2 - start2;
    std::cout << "Time to write to output-file:\t" << diff2.count() << " s" <<"\n";
}