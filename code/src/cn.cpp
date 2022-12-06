#include "cn.hpp"
#include <iostream>
#include <armadillo>
#include <complex>

int k(int i, int j, int L)
// Translate matrix indicies to a single index
{
    return i + j * L;
};



// (!The above code should be implemented in below struct?)
crank_nicolson::crank_nicolson(double xy_steps, double time_step, double time)
{
    h = 1 / xy_steps;
    M = xy_steps;
    dt = time_step;
    T = time;
    N = (M - 2) * (M - 2);
    r = arma::cx_double{0,1} * dt / (2 * h * h);
    a = arma::cx_double{1, 0} + arma::cx_double{4, 0} * r; // 1 + 4r
    b = arma::cx_double{1, 0} - arma::cx_double{4, 0} * r; // 1 - 4r

    V = arma::cx_mat(M - 2, M - 2, arma::fill::zeros);
    U = arma::cx_mat(M - 2, M - 2, arma::fill::zeros);

    A = arma::cx_mat(N, N, arma::fill::zeros);
    B = arma::cx_mat(N, N, arma::fill::zeros);
    
}

void crank_nicolson::init_state_params(
    //initialize as complex numbers with only real parts
    arma::cx_double x_c_,
    arma::cx_double sigma_x_,
    arma::cx_double p_x_,
    arma::cx_double y_c_,
    arma::cx_double sigma_y_,
    arma::cx_double p_y_,
    arma::cx_double v_0_)
{
    x_c = x_c_;
    sigma_x = sigma_x_;
    p_x = p_x_;
    y_c = y_c_;
    sigma_y = sigma_y_;
    p_y = p_y_;
    v_0 = v_0_;
}

void crank_nicolson::init_state()
{
    // Initialize the state
    for (int j = 0; j < M - 2; j++)
    {
        for (int i = 0; i < M - 2; i++)
        {
            double x = i * h;
            double y = j * h;
            U(i, j) = std::exp(-std::pow(x - x_c, 2) / (arma::cx_double{2, 0} \
            * sigma_x * sigma_x) - pow(y - y_c, 2) / (arma::cx_double{2, 0} * sigma_y * sigma_y) \
            + arma::cx_double{0, 1} * p_x * (x - x_c) \
            + arma::cx_double{0, 1} * p_y * (y - y_c));
        }
    }

    // Normalize state
    arma::cx_double pp;
    double norm = 0;
    for (int j = 0; j < M - 2; j++)
    {
        for (int i = 0; i < M - 2; i++)
        {
            pp = U(i, j) * conj(U(i, j));
            norm += pp.real();
        }
    }

    for (int j = 0; j < M - 2; j++)
    {
        for (int i = 0; i < M - 2; i++)
        {
            U(i, j) = U(i, j) / sqrt(norm);
        }
    }

}

// int M, double h, double dt, arma::cx_mat V, char sign
void crank_nicolson::generate_A_B()
{
    // Fill vector with r values
    arma::cx_vec r_vec;

    r_vec = arma::cx_vec(N - 3, arma::fill::ones);
    r_vec = r * r_vec;
    A.diag(3) = r_vec;
    A.diag(-3) = r_vec;

    r_vec = arma::cx_vec(N - 1, arma::fill::ones);
    r_vec = r * r_vec;
    A.diag(1) = r_vec;
    A.diag(-1) = r_vec;

    r_vec = arma::cx_vec(N - 3, arma::fill::ones);
    r_vec = -r * r_vec;
    B.diag(3) = r_vec;
    B.diag(-3) = r_vec;

    r_vec = arma::cx_vec(N - 1, arma::fill::ones);
    r_vec = -r * r_vec;
    B.diag(1) = r_vec;
    B.diag(-1) = r_vec;


    // Removing corner r values (!verify this)
    for (int i = (M - 2); i < N; i += M - 2)
    {
        A(i - 1, i) = arma::cx_double{0, 0};
        A(i, i - 1) = arma::cx_double{0, 0};
        B(i - 1, i) = arma::cx_double{0, 0};
        B(i, i - 1) = arma::cx_double{0, 0};
    }

    update_A_B();
}

void crank_nicolson::update_A_B()
{
    arma::cx_vec diagonal_a = arma::cx_vec(N);
    arma::cx_vec diagonal_b = arma::cx_vec(N);


    for (int j = 0; j < M - 2; j++)
    {
        for (int i = 0; i < M - 2; i++)
        {
            diagonal_a(k(i, j, M - 2)) = a + (arma::cx_double{0, 1} * dt / arma::cx_double{2, 0}) * V(i, j);
            diagonal_b(k(i, j, M - 2)) = b - (arma::cx_double{0, 1} * dt / arma::cx_double{2, 0}) * V(i, j);

        }
    }

    A.diag() = diagonal_a;
    B.diag() = diagonal_b;
}

// Define potential
void crank_nicolson::potential(int slits)
{
    double x;
    double y;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            x = i * h;
            y = j * h;
            if ((y > 0.49) && (y < 0.51))
            {
                if (slits == 1)
                {
                    if ((x > 0.475) && (x < 0.525))
                    {
                        V(i, j) = arma::cx_double{0, 0};
                    }
                    else
                    {
                        V(i, j) = v_0;
                    }
                }

                else if (slits == 2)
                {
                    if (((x > 0.425) && (x < 0.475))||((x > 0.525) && (x < 0.575)))
                    {
                        V(i, j) = arma::cx_double{0, 0};
                    }
                    else
                    {
                        V(i, j) = v_0;
                    }
                }

                else if (slits == 3)
                {
                    if (((x > 3.75)&&(x<4.25))||((x > 4.75)&&(x<5.25))||((x > 5.75)&&(x<6.25)))
                    {
                        V(i, j) = arma::cx_double{0, 0};
                    }
                    else
                    {
                        V(i, j) = v_0;
                    }
                }
            }
        }
    }
}


void crank_nicolson::update_U()
{
    arma::cx_mat U_new = arma::cx_mat(M - 2, M - 2);
    arma::cx_vec u_vec = arma::cx_vec(N);
    U_new = arma::solve(A, B * U);
    U = U_new;
}