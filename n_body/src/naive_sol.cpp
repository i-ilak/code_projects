#include "naive_sol.hpp"

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
    std::string path= "data/data" + name + ".txt";
    std::ofstream out(path);

    out << objects.size() << "\n";
    for(auto stellar_obj : objects){
        out << stellar_obj.get_name() << ",";
    }
    out << "\n";
    auto start2 = std::chrono::high_resolution_clock::now();
    out << res << "\n";
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff2 = end2 - start2;
    std::cout << "Time to write to output-file:\t" << diff2.count() << " s" <<"\n";
}