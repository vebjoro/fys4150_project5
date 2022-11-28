#include <iostream>
#include <armadillo>


int translate_indecies(int i, int j, arma::mat &A);

arma::mat construct_matrix(arma::vec a, double r, int n);

int main(int argc, char *argv[]) {


  int n = 3;
  double r = 0.1;
  arma::vec a = arma::linspace(333, 334, n*n);
  arma::mat A =  construct_matrix(a, r,  n);
  std::cout << A << std::endl;


  return 0;
}


int translate_indecies(int i, int j, arma::mat &A)
  {
  int n = A.n_cols;
  int l = A.n_rows;
  int k = 0;

  k = i + l*j;
  return k;
}

arma::mat construct_matrix(arma::vec a, double r, int n)
{
  int N2 = n*n;
  int j = 0;
  arma::mat central_submatrix = arma::mat(n, n);
  arma::mat offdiag_submatrix = -222*arma::eye(n, n);
  arma::mat big_matrix = arma::mat(N2, N2);

  //sett første matrise
  central_submatrix = arma::diagmat(a(arma::span(0, n-1)));
  central_submatrix.diag(1) = -r*arma::ones(n-1);
  central_submatrix.diag(-1) = -r*arma::ones(n-1);

  for (int i = 1; i<n; i++){
      j = i*n;
     central_submatrix = arma::diagmat(a(arma::span(j, j+n-1)));
     central_submatrix.diag(1) = -r*arma::ones(n-1);
     central_submatrix.diag(-1) = -r*arma::ones(n-1);

     big_matrix(arma::span(j, j+n-1), arma::span(j, j+n-1)) = central_submatrix;
     big_matrix(arma::span(j-n, j-1), arma::span(j-1, j+n-2)) = offdiag_submatrix;

    //std::cout << i << std::endl;
    //big_matrix(arma::span(i*n-1, i*n+1), arma::span(i*n-1, i*n+1)) = central_submatrix;
  }

  return big_matrix;

}
