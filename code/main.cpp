#include <iostream>
#include <armadillo>
#include "cn.hpp"


int main(int argc, char *argv[])
{
    double h, dt, T, n_s;
    arma::cx_double x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0;
    std::string outfile;
    arma::vec out_vec;

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

    //std::cout << system1.U.slice(0) << std::endl;

    system1.generate_A_B();

    for (int i = 1; i < n_s; i++)
    {
        system1.step(i);
    }

    outfile = "data/system1_P_0.bin";
    out_vec = system1.probability_deviation();
    out_vec.save(outfile, arma::arma_binary);


    outfile = "data/system1_U.bin";
    system1.save_simulation(outfile);



    system1.reset_system();

    

    // System 1 with barrier
    sigma_y = 0.10;
    v_0 = arma::cx_double{1e10, 0};
    system1.init_state_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    system1.init_state();

    for (int i = 1; i < n_s; i++)
    {
        system1.step(i);
    }

    outfile = "data/system1_P_e10.bin";
    out_vec = system1.probability_deviation();
    out_vec.save(outfile, arma::arma_binary);

    // System 2
    h = 0.005;
    dt = 2.5e-5;
    T = 0.002;
    crank_nicolson system2 = crank_nicolson(1/h, dt, T);

    x_c = arma::cx_double{0.25, 0};
    sigma_x = arma::cx_double{0.05, 0};
    p_x = arma::cx_double{200, 0};
    y_c = arma::cx_double{0.5, 0};
    sigma_y = arma::cx_double{0.2, 0};
    p_y = arma::cx_double{0, 0};
    v_0 = arma::cx_double{1e10, 0};

    
>>>>>>> 200679572543a714645e4a85d02b4ee9c9df0475

    // arma::cx_mat out_mat = system1.U.slice(0);
    // outfile = "data/system1_U.bin";
    // out_mat.save(outfile, arma::arma_binary);

    

    // // System 2
    // h = 0.005;
    // dt = 2.5e-5;
    // T = 0.002;
    // crank_nicolson system2 = crank_nicolson(1/h, dt, T);
    
    // x_c = arma::cx_double{0.25, 0};
    // sigma_x = arma::cx_double{0.05, 0};
    // p_x = arma::cx_double{200, 0};
    // y_c = arma::cx_double{0.5, 0};
    // sigma_y = arma::cx_double{0.2, 0};
    // p_y = arma::cx_double{0, 0};
    // v_0 = arma::cx_double{1e10, 0};

    // system2.init_state_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    // system2.init_state();

    // for (int i = 0; i < n_s; i++)
    // {
    //     system2.step(i);
    // }












    // crank_nicolson system = crank_nicolson(5, 5, 1);
    // arma::cx_double x_c(0.2, 0);
    // arma::cx_double sigma_x(0.4, 0);
    // arma::cx_double p_x(1, 0);
    // arma::cx_double y_c(0.5, 0);
    // arma::cx_double sigma_y(0.4, 0);
    // arma::cx_double p_y(0, 0);
    // arma::cx_double v_0(0, 0);
    // system.init_state_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
    // arma::cx_mat U(5, 5, arma::fill::zeros);
    // arma::cx_mat U_new(5, 5, arma::fill::zeros);
    // system.init_state(U, U_new);
    // std::string out = "./data/U.bin";
    // U_new.save(out, arma::arma_binary);


    return 0;
}
