#include <iostream>
#include <armadillo>
#include "cn.hpp"


int main(int argc, char *argv[])
{
    // System parameters and variables
    double h, dt, T, n_s;
    arma::cx_double x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0;
    std::string outfile;
    arma::cx_cube out_cube;
    arma::vec out_vec;
    arma::mat out_mat;

    // System 1
    h = 0.005; // Step size
    dt = 2.5e-5; // Time step
    T = 0.008; // Total time

    n_s = T / dt; // Number of time steps

    // Initialize system 1
    crank_nicolson system1 = crank_nicolson(1/h, dt, T);
    
    // Set the wave packet parameters and create the wave packet
    x_c = arma::cx_double{0.25, 0};
    sigma_x = arma::cx_double{0.05, 0};
    p_x = arma::cx_double{200, 0};
    y_c = arma::cx_double{0.5, 0};
    sigma_y = arma::cx_double{0.05, 0};
    p_y = arma::cx_double{0, 0};
    v_0 = arma::cx_double{0, 0};

    system1.wave_packet_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    system1.wave_packet();

    // This system has no potential
    // system1.potential(0)

    // Generate the matrices A and B
    system1.generate_A_B();

    // Evolve the system
    for (int i = 1; i < n_s; i++)
    {
        std::cout << i << std::endl;
        system1.step(i);
    }

    // Save the probability deviation data to file
    outfile = "data/system1_P_0.bin";
    out_vec = system1.probability_deviation();
    out_vec.save(outfile, arma::arma_binary);

    // Reset the system
    system1.reset_system();

    // New wave packet
    sigma_y = 0.10;
    v_0 = arma::cx_double{1e10, 0};
    system1.wave_packet_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    system1.wave_packet();

    // Evolve the system
    for (int i = 1; i < n_s; i++)
    {
        system1.step(i);
    }

    // Save the probability deviation data to file
    outfile = "data/system1_P_1.bin";
    out_vec = system1.probability_deviation();
    out_vec.save(outfile, arma::arma_binary);


    // System 2
    h = 0.005; // Step size
    dt = 2.5e-5; // Time step 
    T = 0.002; // Total time

    n_s = T / dt; // Number of time steps

    // Initialize system 2
    crank_nicolson system2 = crank_nicolson(1/h, dt, T);

    // Set the wave packet parameters and create the wave packet
    x_c = arma::cx_double{0.25, 0};
    sigma_x = arma::cx_double{0.05, 0};
    p_x = arma::cx_double{200, 0};
    y_c = arma::cx_double{0.5, 0};
    sigma_y = arma::cx_double{0.2, 0};
    p_y = arma::cx_double{0, 0};
    v_0 = arma::cx_double{1e10, 0};

    system2.wave_packet_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    system2.wave_packet();

    // Set the potential (2 slits)
    system2.potential(2);

    // Generate the matrices A and B
    system2.generate_A_B();

    // Evolve the system
    for (int i = 1; i < n_s; i++)
    {
        system2.step(i);
    }

    // Save U and V to file
    outfile = "data/system2_U.bin";
    out_cube = system2.U;
    out_cube.save(outfile, arma::arma_binary);

    outfile = "data/system2_V.bin";
    out_mat = arma::real(system2.V);
    out_mat.save(outfile, arma::arma_binary);

    system2.reset_system();


    // System 2.1 (1 slit)
    // Create the wave packet
    system2.wave_packet();

    // Set the potential (1 slit)
    system2.potential(1);

    // Generate the matrices A and B
    system2.generate_A_B();

    // Evolve the system
    for (int i = 1; i < n_s; i++)
    {
        system2.step(i);
    }

    // Save U and V to file
    outfile = "data/system2.1_U.bin";
    out_cube = system2.U;
    out_cube.save(outfile, arma::arma_binary);

    outfile = "data/system2.1_V.bin";
    out_mat = arma::real(system2.V);
    out_mat.save(outfile, arma::arma_binary);

    system2.reset_system();

    // System 2.3 (3 slits)
    // Create the wave packet
    system2.wave_packet();

    // Set the potential (3 slits)
    system2.potential(3);
    
    // Generate the matrices A and B
    system2.generate_A_B();

    // Evolve the system
    for (int i = 1; i < n_s; i++)
    {
        system2.step(i);
    }

    // Save U and V to file
    outfile = "data/system2.3_U.bin";
    out_cube = system2.U;
    out_cube.save(outfile, arma::arma_binary);

    outfile = "data/system2.3_V.bin";
    out_mat = arma::real(system2.V);
    out_mat.save(outfile, arma::arma_binary);

    return 0;
}
