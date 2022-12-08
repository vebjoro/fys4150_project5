#include "cn.hpp"
#include <iostream>
#include <armadillo>
#include <complex>

int k(int i, int j, int L)
// Translate matrix indicies to a single index (internal points)
{
    return (i-1) + (j-1) * L;
};



// (!The above code should be implemented in below struct?)
crank_nicolson::crank_nicolson(double xy_steps, double time_step, double time)
{
    h = 1 / xy_steps;
    M = xy_steps;
    dt = time_step;
    T = time;
    n_time_steps = std::round(T / dt);
    N = (M - 2) * (M - 2);
    r = arma::cx_double{0,1} * dt / (2 * h * h);
    a = arma::cx_double{1, 0} + arma::cx_double{4, 0} * r; // 1 + 4r
    b = arma::cx_double{1, 0} - arma::cx_double{4, 0} * r; // 1 - 4r

    V = arma::cx_mat(M, M, arma::fill::zeros);
    U = arma::cx_cube(M, M, n_time_steps, arma::fill::zeros);



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
    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 1; j <= M - 2; j++)
        {

            double x = i * h;
            double y = j * h;
            U(i, j, 0) = std::exp(-std::pow(x - x_c, 2) / (arma::cx_double{2, 0} \
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
            pp = U(i, j, 0) * conj(U(i, j, 0));
            norm += pp.real();
        }
    }

    for (int j = 0; j < M - 2; j++)
    {
        for (int i = 0; i < M - 2; i++)
        {
            U(i, j, 0) = U(i, j, 0) / sqrt(norm);
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

    arma::cx_vec diagonal_a = arma::cx_vec(N);
    arma::cx_vec diagonal_b = arma::cx_vec(N);


    for (int i = 1; i <= M-2; i++)
    {
        for (int j = 1; j < M-2; j++)
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

    if (slits == 0)
    {
        V = arma::cx_mat(M, M, arma::fill::zeros);
    }

    else
    {
        V = arma::cx_mat(M, M, arma::fill::ones);
        V = V * v_0;
    }

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
                    if ((x < 0.475) && (x > 0.525))
                    {
                        V(i, j) = v_0;
                    }
                }

                else if (slits == 2)
                {
                    if (((x < 0.425) && (x > 0.475))||((x < 0.525) && (x > 0.575)))
                    {
                        V(i, j) = v_0;
                    }
                }

                else if (slits == 3)
                {
                    if (((x < 3.75)&&(x > 4.25))||((x < 4.75)&&(x > 5.25))||((x < 5.75)&&(x > 6.25)))
                    {
                        V(i, j) = v_0;
                    }
                }
            }
        }
    }
}


void crank_nicolson::step(int n)
{
    arma::cx_mat U_new = U.slice(n - 1);
    arma::cx_vec u_vec = arma::cx_vec(N);
    // Print length of u_vec


    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 1; j <= M-2; j++)
        {
            u_vec(k(i, j, M-2)) = U(i, j, n-1);

        }
    }
    

    u_vec = arma::solve(A, B * u_vec);
    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 1; j <= M - 2; j++)
        {
            U_new(i, j) = u_vec(k(i, j, M-2));
            
        }
    }
    U.slice(n) = U_new;
}

void crank_nicolson::reset_system()
{
    V = arma::cx_mat(M, M, arma::fill::zeros);
    U = arma::cx_cube(M, M, n_time_steps, arma::fill::zeros);
}

arma::vec crank_nicolson::probability()
{
    arma::vec out_vec = arma::vec(n_time_steps);
    for (int i = 0; i < n_time_steps; i++)
    {
        out_vec(i) = arma::accu(arma::abs(U.slice(i)) % arma::abs(U.slice(i)));
    }
    return out_vec;
}


