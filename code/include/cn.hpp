#ifndef STATE_CN_
#define STATE_CN_

#include <armadillo>

int k(int i, int j, int M);
arma::cx_mat generate_matrix(int M, double h, double dt, arma::cx_mat V);
#endif // STATE_CN_
