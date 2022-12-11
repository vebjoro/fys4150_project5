#include <iostream>
#include <armadillo>
#include "cn.hpp"


int main(int argc, char *argv[])
{
    double h, dt, T, n_s;
    arma::cx_double x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0;
    std::string outfile;
    arma::vec out_vec;
    arma::mat out_mat;

    // System 1
    h = 0.005;
    dt = 2.5e-5;
    T = 0.008;

    n_s = T / dt;


    crank_nicolson system1 = crank_nicolson(1/h, dt, T);
    
    x_c = arma::cx_double{0.25, 0};
    sigma_x = arma::cx_double{0.05, 0};
    p_x = arma::cx_double{200, 0};
    y_c = arma::cx_double{0.5, 0};
    sigma_y = arma::cx_double{0.05, 0};
    p_y = arma::cx_double{0, 0};
    v_0 = arma::cx_double{0, 0};

    system1.init_state_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    system1.init_state();
    system1.generate_A_B();

    for (int i = 1; i < n_s; i++)
    {
        std::cout << i << std::endl;
        system1.step(i);
    }

    outfile = "data/system1_P_0.bin";
    out_vec = system1.probability_deviation();
    out_vec.save(outfile, arma::arma_binary);

    system1.reset_system();

    

    // System 1 with barrier
    sigma_y = 0.10;
    v_0 = arma::cx_double{1e10, 0};
    system1.init_state_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    system1.init_state();

    for (int i = 1; i < n_s; i++)
    {
        std::cout << i << std::endl;

        system1.step(i);
    }

    outfile = "data/system1_P_1.bin";
    out_vec = system1.probability_deviation();
    out_vec.save(outfile, arma::arma_binary);

    // // System 2
    h = 0.005;
    dt = 2.5e-5;
    T = 0.002;


    n_s = T / dt;
    crank_nicolson system2 = crank_nicolson(1/h, dt, T);

    x_c = arma::cx_double{0.25, 0};
    sigma_x = arma::cx_double{0.05, 0};
    p_x = arma::cx_double{200, 0};
    y_c = arma::cx_double{0.5, 0};
    sigma_y = arma::cx_double{0.2, 0};
    p_y = arma::cx_double{0, 0};
    v_0 = arma::cx_double{1e10, 0};

    system2.init_state_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    system2.potential(2);
    system2.init_state();


    system2.generate_A_B();

    for (int i = 1; i < n_s; i++)
    {

        system2.step(i);
    }

    outfile = "data/system2_U.bin";
    arma::cx_cube out_cube = system2.U;
    out_cube.save(outfile, arma::arma_binary);

    outfile = "data/system2_V.bin";
    out_mat = arma::real(system2.V);
    out_mat.save(outfile, arma::arma_binary);





    return 0;
}
