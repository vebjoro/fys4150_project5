#include "cn.hpp"
#include <iostream>
#include <armadillo>

int k(int i, int j, int M)
// Translate matrix indicies to a single index
{
    return i + j * M;
};

arma::cx_mat generate_matrix(int M, double h, double dt, arma::cx_mat V, char sign)
{
    // Edge length
    int N = (M - 2) * (M - 2);

    // Filling the matrix with r values
    arma::cx_mat A = arma::cx_mat(N, N, arma::fill::zeros);

    double r = dt / (2 * h * h);
    // Fill vector with r values
    arma::cx_vec r_vec;
    if (sign == 'B')
    {

        r_vec = arma::cx_vec(N - 3, arma::fill::ones);
        r_vec = -r * r_vec;
        A.diag(3) = r_vec;
        A.diag(-3) = r_vec;

        r_vec = arma::cx_vec(N - 1, arma::fill::ones);
        r_vec = -r * r_vec;
        A.diag(1) = r_vec;
        A.diag(-1) = r_vec;

        arma::cx_double a = 1 + 4 * r + arma::cx_double{0, 1} * dt / arma::cx_double{2, 0};
    }
    else
    {
        r_vec = arma::cx_vec(N - 3, arma::fill::ones);
        r_vec = r * r_vec;
        A.diag(3) = r_vec;
        A.diag(-3) = r_vec;

        r_vec = arma::cx_vec(N - 1, arma::fill::ones);
        r_vec = r * r_vec;
        A.diag(1) = r_vec;
        A.diag(-1) = r_vec;
    }

    // Removing corner r values (!verify this)
    for (int i = (M - 2); i < N; i += M - 2)
    {
        A(i - 1, i) = arma::cx_double{0, 0};
        A(i, i - 1) = arma::cx_double{0, 0};
    }

    // Filling the diagonal
    arma::cx_vec diagonal = arma::cx_vec(N);
    arma::cx_double a = 1 + 4 * r + arma::cx_double{0, 1} * dt / arma::cx_double{2, 0};
    arma::cx_double b = 1 - 4 * r - arma::cx_double{0, 1} * dt / arma::cx_double{2, 0};
    if (sign == 'B')
    {
        for (int j = 0; j < M - 2; j++)
        {

            for (int i = 0; i < M - 2; i++)
            {
                diagonal(k(i, j, M - 2)) = b * V(i, j);
            }
        }
    }
    else
    {
        for (int j = 0; j < M - 2; j++)
        {

            for (int i = 0; i < M - 2; i++)
            {
                diagonal(k(i, j, M - 2)) = a * V(i, j);
            }
        }
    }
    A.diag() = diagonal;
    return A;
}
