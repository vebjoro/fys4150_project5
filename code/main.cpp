#include <iostream>
#include <armadillo>


int translate_indecies(int i, int j, arma::mat &A);

arma::mat construct_matrix(arma::vec a, double r, int n);

int main(int argc, char *argv[]) {


  int n = 3;
  double r = 0.1;
  arma::vec a = arma::linspace(0, 1, n);
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
  arma::mat central_submatrix = arma::diagmat(a);
  central_submatrix.diag(1) = -r*arma::ones(n-1);
  central_submatrix.diag(-1) = -r*arma::ones(n-1);

  arma::mat big_matrix = arma::mat(N2, N2);
  for (int i = 0; i<n; i++){
    std::cout << i << std::endl;
    //big_matrix(arma::span(i*n-1, i*n+1), arma::span(i*n-1, i*n+1)) = central_submatrix;
  }

  return central_submatrix;

}
