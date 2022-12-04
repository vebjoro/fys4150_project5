#ifndef STATE_CN_
#define STATE_CN_

#include <armadillo>

int k(int i, int j, int M);
arma::cx_mat generate_matrix(int M, double h, double dt, arma::cx_mat V, char sign);

struct crank_nicolson
{
    double h;
    double M;
    double dt;
    double T;

    arma::cx_double x_c;
    arma::cx_double sigma_x;
    arma::cx_double p_x;
    arma::cx_double y_c;
    arma::cx_double sigma_y;
    arma::cx_double p_y;
    arma::cx_double v_0;

    crank_nicolson(double h, double dt, double T);

    void init_state_params(
        arma::cx_double x_c,
        arma::cx_double sigma_x,
        arma::cx_double p_x,
        arma::cx_double y_c,
        arma::cx_double sigma_y,
        arma::cx_double p_y,
        arma::cx_double v_0);

    void init_state(arma::cx_mat &V, arma::cx_mat &V_new);
};
#endif // STATE_CN_
