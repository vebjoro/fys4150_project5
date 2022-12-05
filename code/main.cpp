#include <iostream>
#include <armadillo>
#include "cn.hpp"

int translate_indecies(int i, int j, arma::cx_mat &A);

arma::cx_mat construct_matrix(arma::cx_vec v, double r, int n, char m);

int main(int argc, char *argv[])
{

    //   int n = 3;
    //   double r = 0.1;
    //   arma::vec a_r = arma::ones(n*n)*333;
    //   arma::vec a_i = arma::ones(n*n)*444;
    //   arma::cx_vec a = arma::cx_vec(a_r, a_i);
    //   //std::cout <<  a << std::endl;
    //   arma::cx_mat A =  construct_matrix(a, r,  n, 'A');
    //   std::cout << A << std::endl;

    // arma::cx_mat V = arma::cx_mat(3, 3, arma::fill::randn);
    // std::cout << V << std::endl;
    // arma::cx_mat A = generate_matrix(5, 1, 1, V, 'A');
    // std::cout << A << std::endl;



    return 0;
}

// int translate_indecies(int i, int j, arma::cx_mat &A)
//   {
//   int n = A.n_cols;
//   int l = A.n_rows;
//   int k = 0;

//   k = i + l*j;
//   return k;
// }

// arma::cx_mat construct_matrix(arma::cx_vec v, double r, int n, char m)
// {
//   if (m == 'B') //tells which matrix to construct
//     {
//       r = -r;
//     }
//   int N2 = n*n;
//   int j = 0;
//   arma::cx_mat central_submatrix = arma::cx_mat(n, n);
//   arma::mat offdiag_imag = -r*arma::eye(n, n);
//   arma::mat offdiag_real = arma::mat(n, n).fill(0);
//   arma::cx_mat offdiag_submatrix = arma::cx_mat(offdiag_real, offdiag_imag);
//   arma::cx_mat big_matrix = arma::cx_mat(N2, N2);

//   central_submatrix = arma::diagmat(v(arma::span(0, n-1)));
//   arma::imag(central_submatrix.diag(-1)) = {2, 2};
//   //std::cout << typeid(arma::imag(central_submatrix.diag(-1))).name() << std::endl; //.set_imag(arma::ones(n-1));
//   //std::cout << typeid(arma::imag(central_submatrix)).name() << std::endl;;// = -r*arma::ones(n-1);
//   //central_submatrix.imag.diag(-1) = -r*arma::ones(n-1);

//   //big_matrix(arma::span(0, n-1), arma::span(0, n-1)) = central_submatrix; //sett første matrise

//   for (int i = 1; i<n; i++){
//       j = i*n;
//      //central_submatrix = arma::diagmat(v(arma::span(j, j+n-1)));
//      //central_submatrix.diag(1) = -r*arma::ones(n-1);
//      //central_submatrix.diag(-1) = -r*arma::ones(n-1);

//      //big_matrix(arma::span(j, j+n-1), arma::span(j, j+n-1)) = central_submatrix;
//      //big_matrix(arma::span(j-n, j-1), arma::span(j, j+n-1)) = offdiag_submatrix;
//      //big_matrix(arma::span(j, j+n-1), arma::span(j-n, j-1)) = offdiag_submatrix;

//   }

//   return central_submatrix; // big_matrix;

// }
