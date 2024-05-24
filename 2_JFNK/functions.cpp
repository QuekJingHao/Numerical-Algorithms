# include <iostream>
# include <string>
# include <fstream>
# include <cmath>
# include <iomanip>
# include <chrono>
# include <utility>
# include <functional>
# include <eigen3/Eigen/Dense>
# include "functions.h"

using namespace Eigen;


//  residual form of the nonlinear equations 
VectorXf F(double x, double y)
{   
    VectorXf b(2);

    b(0) = sin(x) + 3.0*cos(y) - 2.0;
    b(1) = cos(x) - sin(y) + 0.2;

    return b;
}

//*******************************************************************************************************************//

//                DO NOT CHANGE ANY OF THE FOLLOWING FUNCTIONS!!                                                     //

/********************************************************************************************************************
 * 		                     Jacobian approx
 * Approximates the Jacobian matrix using finite difference formula
*********************************************************************************************************************/
MatrixXf jacobian_approx(VectorXf X, std::function<VectorXf(double, double)>func)
{
    //  data declarations
    int n = X.rows();
    double b = 0.000001;
    MatrixXf J(n, n);
    MatrixXf e;

    J.setZero();
    e.setIdentity(n, n);  //  set an identital matrix for use later

    for (int i{0}; i<n; i++)
    {
        for (int j{0}; j<n; j++)
        {   
            double episilon_j = b * X(j) + b;
            VectorXf X_forward = X + e.col(j) * episilon_j;
            VectorXf func_forward = func(X_forward(0), X_forward(1));
            VectorXf func_backward = func(X(0), X(1));

            J(i, j) = (func_forward(i) - func_backward(i)) / episilon_j;
        }
    }

    return J;

}



/********************************************************************************************************************
 * 		                     Arnoldi iteration
 *           Performs the Arnoldi iteration on the system Ax = b
*********************************************************************************************************************/
std::pair<MatrixXf, MatrixXf> arnoldi_iteration(MatrixXf A, MatrixXf b, int n)
{
    //  declare data
    auto m = A.rows();
    MatrixXf q, h;
    VectorXf v; // you have to do declare the temp v vector

    //  set q and h as zeros
    q.setZero(m, n+1);
    h.setZero(n+1, n);

    // set the first Krylov vector
    q.col(0) = b / b.norm();

    for (int i{0}; i<n; i++)
    {   
        v = A * q.col(i);
        
        for (int j{0}; j<=i; j++)
        {   
            h(j, i) = q.col(j).dot(v);
            v = v - h(j, i) * q.col(j);
        }
            h(i+1, i) = v.norm();
            q.col(i+1) = v / h(i+1, i);
        
    }

    return std::make_pair(q, h);
}


/********************************************************************************************************************
 * 		                     Givens rotation
 *             Perform the Givens rotation on input matrix
*********************************************************************************************************************/
std::pair<MatrixXf, VectorXf> givens_rotation(MatrixXf A, VectorXf beta, int k)
{   
    /*
    @param m is the number of iterations - should be the loop counter if this function is nested inside GMRES
    note that the indexing starts from 0!
    */

    //  declare variables
    MatrixXf omega; 
    int m = A.rows();
    omega.setIdentity(k+1, k+1);

    for (int j{0}; j<k; j++)
    {
        double c = A(j, j) / sqrt(A(j, j)*A(j, j) + A(j+1, j)*A(j+1, j));
        double s = A(j+1, j) / sqrt(A(j, j)*A(j, j) + A(j+1, j)*A(j+1, j));

        omega(j, j) = c;
        omega(j, j+1) = s;
        omega(j+1, j) = -s;
        omega(j+1, j+1) = c;

        // std::cout << "The givens matrix is \n\n" << omega << "\n\n";
        //  apply the givens matrix on A and delete the very last element of the rotated matrix
        A = omega * A;
        A(j+1, j) = 0.0;

        // for the RHS beta vector
        beta = omega * beta;

        // std::cout << "Applying the Givens matrix on A, we have \n\n" << A << "\n\n";
        // std::cout << "Applying the Givens matrix on beta vector, we have \n\n" << beta << "\n\n";

        //  reset the givens matrix for the next iteration
        omega.setIdentity(k+1, k+1);
    }

    return std::make_pair(A, beta);

}



/********************************************************************************************************************
 * 		                                      GMRES
 *    Uses the Generalized Mininum Residuals (GMRES) method to solve for x in the system Ax = b.
*********************************************************************************************************************/
VectorXf gmres(MatrixXf A, MatrixXf b, VectorXf x, int max_iterations, double error)
{   
    //  data and variable declarations
    int m = A.rows();
    int n = max_iterations;

    std::pair<MatrixXf, MatrixXf> arnoldi_pair;
    MatrixXf Q;          // Q: (m) x (n+1)     Here the n refers to the steps of the GMRES iteration!!!!
    MatrixXf Q_temp;     // (m) x n
    MatrixXf H;          // H: (n+1) x n
    MatrixXf R_tilda;    // R_tilda: (n+1) x n
    MatrixXf R;          // R: n by n 
    VectorXf r;          // r: m x 1
    VectorXf beta;       // beta: (n+1) x 1 this will be the axis vector multiplied by the norm
    VectorXf g_tilda;    // g_bar: (n+1) x 1
    VectorXf g;          // g: n x 1
    VectorXf y;          // y: n x 1
    double residual;     // double


    //  print out message at the begining
    std::cout << "\t\t This program executes the Generalized Mininum Residuals (GMRES) Algorithm. \n\n";

    // start with r0, the first guessed Krylov vector
    r = b - A * x;
    double b_norm = b.norm();
    residual = r.norm() / b_norm;

    // std::cout << "The inital Krylov vector is \n\n" << r << "\n\n";

    for (int i{1}; i<=n; i++)
    {   
        //  run Arnoldi iteration
        arnoldi_pair = arnoldi_iteration(A, r/r.norm(), i);
        Q = arnoldi_pair.first;
        H = arnoldi_pair.second;


        //  Set up the RHS of the equation
        beta.setZero(i+1);
        beta(0) = r.norm();
        

        //  apply the givens rotation on the Hessenberg matrix and the RHS
        R_tilda = givens_rotation(H, beta, i).first;
        g_tilda = givens_rotation(H, beta, i).second;


        //  delete the last row of the R tilda matrix and g tilda vector
        R = R_tilda.topRows(i);
        g = g_tilda.head(i);

        //  find the residual
        residual = abs(g_tilda(i)) / b_norm;
        std::cout << "The residual is \n\n" << residual << "\n\n";

        // if the residual is smaller than the error, break out of this loop and compute the value of x
        if (residual <= error) 
        {
            break;
        }

        // this will be what is going to be multipled y and add to solution vector x
        Q_temp = Q.leftCols(i+1);
        
        /*
        std::cout << "At iteration k = \t" << i << "\t The Arnoldi matrix Q is \n\n" << Q << "\n\n";
        std::cout << "Getting rid of the last column of Q is \n\n" << Q.leftCols(i) << "\n\n";
        std::cout << "The Hessenberg matrix h is \n\n" << H << "\n\n";
        std::cout << "After applying the Givens rotation, we have \n\n" << R_tilda <<  "\n\n";
        std::cout << "Deleting the last row of R tilda, R is \n\n" << R << "\n\n";
        std::cout << "The g tilda vector is \n\n" << g_tilda << "\n\n";
        std::cout << "Deleting the last row of g tilda, g is \n\n" << g << "\n\n"; */

        std::cout << "----------------------------------------------------------- END ----------------------------------------------------------- \n\n";
        
    }

    // calculate yn
    y = R.inverse() * g;
    //std::cout << "The least squares vector y is \n\n" << y << "\n\n";
    //std::cout << "The Arnoldi matrix Q is \n\n" << Q_temp << "\n\n";

    // calculate x
    x = x + Q_temp * y;

    // std::cout << "Q temp is \n\n" << Q_temp << "\n\n";
    // std::cout << "The solution vector is \n\n" << x << "\n\n"; 

    return x;
}


/********************************************************************************************************************
 * 		                                 JFNK
 *  Performs the Jacobian-Free Newton Krylov method for solving system of nonlinear equations
*********************************************************************************************************************/
VectorXf jfnk(std::function<VectorXf(double, double)>func, VectorXf X_initial, double err, int max_iterations)
{   

    for (int i{0}; i<max_iterations; i++)
    {   
        VectorXf b_right = -func(X_initial(0), X_initial(1));

        MatrixXf J = jacobian_approx(X_initial, func);

        VectorXf v = gmres(J, b_right, X_initial, 10, err);

        X_initial = X_initial + v;

        double err = sqrt( X_initial.dot(X_initial) / X_initial.rows() );
        double err_tol = 0.000001 / v.norm();

        if (err < err_tol)
        {
            break;
        }
        
        else if (err > err_tol)
        {
            continue;
        }

        else
        {
            std::cout << "JFNK failed to converge! Check tolerance and maximum iterations. \n\n"; 
        }
    }

    return X_initial;
}
