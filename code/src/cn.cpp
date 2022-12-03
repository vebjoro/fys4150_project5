#include "cn.hpp"
#include <iostream>
#include <armadillo>

int k(int i, int j, int M)
// Translate inner matrix indicies to a single index
{
    return (i - 1) + (j - 1) * M;
};

arma::cx_mat generate_matrix(int M, double h, double dt, arma::cx_mat V)
{
    int N = (M - 2) * (M - 2);
    arma::mat real = arma::mat(N, N, arma::fill::zeros);
    arma::mat imag = arma::mat(N, N, arma::fill::zeros);

    double r_imag = dt / (2 * h * h);
    imag.diag(3) = r_imag * arma::ones(N - 3);
    imag.diag(-3) = r_imag * arma::ones(N - 3);
    imag.diag(1) = r_imag * arma::ones(N - 1);
    imag.diag(-1) = r_imag * arma::ones(N - 1);

    for (int i = (M - 2); i < N; i += M - 2)
    {
        imag(i - 1, i) = 0;
        imag(i, i - 1) = 0;
    }
    arma::cx_mat A = arma::cx_mat(real, imag);

    for (int i = 0; i < M - 2; i++)
    {
        for (int j = 0; j < M - 2; j++)
        {
            A(i, j) += V(k(i, j, M));
        }
    }
}
return A;
}
