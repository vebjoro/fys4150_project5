#include <iostream>
#include <armadillo>
#include "cn.hpp"


int main(int argc, char *argv[])
{

crank_nicolson system = crank_nicolson(5, 5, 1);
arma::cx_double x_c(0.2, 0);
arma::cx_double sigma_x(0.4, 0);
arma::cx_double p_x(1, 0);
arma::cx_double y_c(0.5, 0);
arma::cx_double sigma_y(0.4, 0);
arma::cx_double p_y(0, 0);
arma::cx_double v_0(0, 0);
 system.init_state_params(x_c, sigma_x, p_x, y_c, sigma_y, p_y, v_0);
arma::cx_mat U(5, 5, arma::fill::zeros);
arma::cx_mat U_new(5, 5, arma::fill::zeros);
system.init_state(U, U_new);
std::string out = "./data/U.bin";
U_new.save(out, arma::arma_binary);


    return 0;
}
