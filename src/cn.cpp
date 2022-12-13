#include "cn.hpp"
#include <iostream>
#include <armadillo>
#include <complex>
#include <cmath>


int k(int i, int j, int L)
// Translate matrix indicies to a single index (internal points)
{
    return (i-1)*L + (j-1);
};


crank_nicolson::crank_nicolson(double xy_steps, double time_step, double time)
// Initialize the system with system parameters
{

    h = 1 / xy_steps;    // Step size
    M = xy_steps;        // System size in x and y
    dt = time_step;      // Time step
    T = time;            // Total time
    n_time_steps = std::round(T / dt); // Number of time steps
    N = (M - 2) * (M - 2);             // Number of internal points
    r = arma::cx_double{0,1} * dt / (2 * h * h);
    a = arma::cx_double{1, 0} + arma::cx_double{4, 0} * r; // 1 + 4r
    b = arma::cx_double{1, 0} - arma::cx_double{4, 0} * r; // 1 - 4r

    // Initialize the matrices
    V = arma::cx_mat(M, M, arma::fill::zeros); // Potential matrix
    U = arma::cx_cube(M, M, n_time_steps, arma::fill::zeros); // State cube (x, y, n) (3D matrix)
    A = arma::sp_cx_mat(N, N); 
    B = arma::sp_cx_mat(N, N);
}


void crank_nicolson::wave_packet_params(
// Initialize the wave function parameters
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


void crank_nicolson::wave_packet()
// Initialize the initial state (U[0]), create a wave packet
{
    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 1; j <= M - 2; j++)
        {

            arma::cx_double x = arma::cx_double{i * h, 0};
            arma::cx_double y = arma::cx_double{j * h, 0};

            // Wave packet function
            U(i, j, 0) = std::exp(-(std::pow(x - x_c, 2) / (arma::cx_double{2, 0} \
            * sigma_x * sigma_x)) - (std::pow(y - y_c, 2) / (arma::cx_double{2, 0} * sigma_y * sigma_y)) \
            + arma::cx_double{0, 1} * p_x * (x - x_c) \
            + arma::cx_double{0, 1} * p_y * (y - y_c));
        }
    }
    // Finding the norm of the wave packet
    arma::cx_double pp;
    double norm = 0;
    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 0; j <= M - 2; j++)
        {
            pp = U(i, j, 0) * conj(U(i, j, 0));
            norm += pp.real();
        }
    }

    // Normalizing the wave packet
    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 1; j <= M - 2; j++)
        {
            U(i, j, 0) = U(i, j, 0) /sqrt(norm);
        }
    }


}


void crank_nicolson::generate_A_B()
// Generate the A and B matrices for the Crank-Nicolson method
{
    arma::cx_vec r_vec;

    // Diagonal elements
    r_vec = arma::cx_vec(N - (M-2), arma::fill::ones);
    r_vec = -r * r_vec;
    A.diag(M-2) = r_vec;
    A.diag(-(M-2)) = r_vec;

    r_vec = arma::cx_vec(N - 1, arma::fill::ones);
    r_vec = -r * r_vec;
    A.diag(1) = r_vec;
    A.diag(-1) = r_vec;

    r_vec = arma::cx_vec(N - (M-2), arma::fill::ones);
    r_vec = r * r_vec;
    B.diag(M-2) = r_vec;
    B.diag(-(M-2)) = r_vec;

    r_vec = arma::cx_vec(N - 1, arma::fill::ones);
    r_vec = r * r_vec;
    B.diag(1) = r_vec;
    B.diag(-1) = r_vec;

    


    // Removing r values according to the submatrix structure of A and B
    for (int i = (M - 2); i < N; i += M - 2)
    {
    
        A(i - 1, i) = arma::cx_double{0, 0};
        A(i, i - 1) = arma::cx_double{0, 0};
        B(i - 1, i) = arma::cx_double{0, 0};
        B(i, i - 1) = arma::cx_double{0, 0};
    }

    arma::cx_vec diagonal_a = arma::cx_vec(N);
    arma::cx_vec diagonal_b = arma::cx_vec(N);
    

    // Setting the diagonal elements of A and B
    for (int i = 1; i <= M-2; i++)
    {
        for (int j = 1; j <= M-2; j++)
        {

            diagonal_a(k(i, j, M - 2)) = a + (arma::cx_double{0, 1} * dt / arma::cx_double{2, 0}) * V(i, j);
            diagonal_b(k(i, j, M - 2)) = b - (arma::cx_double{0, 1} * dt / arma::cx_double{2, 0}) * V(i, j);
        }
    }
    A.diag() = diagonal_a;
    B.diag() = diagonal_b;
}


void crank_nicolson::potential(int slits)
// Create the potential, with the number of slits as an input. The potential
// is set to 0 in the region where the wave packet is allowed to propagate and
// to v_0 in the region where we have our "wall".
{
    double x;
    double y;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            x = i * h;
            y = j * h;

            if ((x > 0.49) && (x < 0.51))
            {
                if (slits == 1)
                {
                    if ((y < 0.475)||(y > 0.525))
                    {
                        V(i, j) = v_0;
                    }
                }
                else if (slits == 2)
                {
                    if ((y < 0.425) || (y > 0.475&&y < 0.525) || (y > 0.575))
                    {
                        V(i, j) = v_0;
                    }
                }
                else if (slits == 3)
                {
                    if (((y < 0.375)||(y > 0.425&&y < 0.475))||((y > 0.525&&y < 0.575)||(y > 0.625)))
                    {
                        V(i, j) = v_0;
                    }
                }
            }
        }
    }
}


void crank_nicolson::step(int n)
// Perform one step of the Crank-Nicolson method. The input is the number of
// the current step.
{
    arma::cx_mat U_new = U.slice(n - 1); // Copying the previous step
    arma::cx_vec u_vec = arma::cx_vec(N); // Vector for the solution of the linear system



    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 1; j <= M - 2; j++)
        {
            u_vec(k(i, j, M-2)) = U(i, j, n-1);

        }
    }
    
    // Solving the linear system
    u_vec = arma::spsolve(A,B * u_vec);

    // Updating the solution from the vector to the matrix
    for (int i = 1; i <= M - 2; i++)
    {
        for (int j = 1; j <= M - 2; j++)
        {
            U_new(i, j) = u_vec(k(i, j, M-2));
        }
    }

    // Add solution to the cube for storage
    U.slice(n) = U_new;
}


void crank_nicolson::reset_system()
// Reset V and U to zero matrices
{
    V = arma::cx_mat(M, M, arma::fill::zeros);
    U = arma::cx_cube(M, M, n_time_steps, arma::fill::zeros);
}


arma::vec crank_nicolson::probability_deviation()
// Calculate the deviation of the probability from 1 as a function of steps.
{
    arma::vec out_vec = arma::vec(n_time_steps);
    for (int i = 0; i < n_time_steps; i++)
    {
        out_vec(i) = arma::accu(arma::abs(U.slice(i)) % arma::abs(U.slice(i)))-1;
    }
    return out_vec;
}




