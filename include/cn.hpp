//testing github setup
#ifndef STATE_CN_
#define STATE_CN_

#include <armadillo>

int k(int i, int j, int L);

struct crank_nicolson
{
    double h;
    double M;
    double dt;
    double T;
    double n_time_steps;
    double N;
    arma::cx_double r;
    arma::cx_double a;
    arma::cx_double b;

    arma::cx_mat V;
    arma::cx_cube U;
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;

    arma::cx_double x_c;
    arma::cx_double sigma_x;
    arma::cx_double p_x;
    arma::cx_double y_c;
    arma::cx_double sigma_y;
    arma::cx_double p_y;
    arma::cx_double v_0;

    crank_nicolson(double h, double dt, double T);

    void wave_packet_params(
        arma::cx_double x_c,
        arma::cx_double sigma_x,
        arma::cx_double p_x,
        arma::cx_double y_c,
        arma::cx_double sigma_y,
        arma::cx_double p_y,
        arma::cx_double v_0);

    void wave_packet();

    void generate_A_B();

    void potential(int slits);

    void step(int n);

    void reset_system();

    arma::vec probability_deviation();
};
#endif // STATE_CN_
