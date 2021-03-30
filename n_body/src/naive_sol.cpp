#include "naive_sol.hpp"
#include <H5Cpp.h>

void perform_simulation(planet_container_type const & objects,
                        double const & T, int const & N, double const & G,
                        int const & method){
    std::string name;
    switch (method)
    {
    case 1:
        name = "_ee";
        break;
    case 2:
        name = "_em";
        break;
    case 3:
        name = "_vv";
        break;
    }
    // Running the calculation and timing it
    auto start = std::chrono::high_resolution_clock::now();
    auto res = n_body_solver(objects, T, N, G, method);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << name.substr(1,name.size()-1) + " calculation took:\t\t" << diff.count() << " s" <<"\n";

    // Writing the data into a file for later plotting
    const std::string file_name = "./data/naive" + name + ".hdf5";
    const H5std_string  FILE_NAME( file_name );
    H5::H5File file( FILE_NAME, H5F_ACC_TRUNC );
    const int   NX = res.rows();                    // dataset dimensions
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
                data_arr[j][k] = res(j,3*i+k);
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
        time_arr[j] = res(j,res.cols()-1);
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