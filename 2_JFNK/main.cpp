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

/**
 *  This program performs the Jacobian-Free Newton Krylov (JFNK) method for solving solution 
 *  vector of sparse matrices.
 * 
 *  Purpose:  Performs JFNK calculation to solve for solution vector, using Krylov subspace method.
 *            The JFNK method is widely used in nuclear reactor physics calculations, 
 *            where usually, the analytical Jacobian matrix cannot be easily obtained. 
 *            As such, one can use the finite difference approximation to get an estimated Jacobian matrix. 
 *            The classical Newton method for solving system of nonlinear equations is written as
 * 
 *            J * delta-x = -F,
 * 
 *            where J is the Jacobian matrix, delta-x is the correction vector, in which after every iteration, 
 *            one makes the assignment (x = x + delta-x), 
 *            where x is the initial guess vector of the variables in the system
 * 
 *            As mentioned above, it is difficult to get J. After approximating J, 
 *            we can exploit the Kryolov subspace method, and recognize that the system we are solving above
 *            is in the following form
 * 
 *            Ax = b.
 * 
 *            Here, 'x' is our delta x. This equation is then solved by Generalized Mininum Residuals (GMRES) method
 *            for x. Lastly, we make the correction simply by (x = x + delta x). We iterate this procedure until 
 *            convergence. 
 * 
 *  
 *   Usage:   Recommended to open this entire folder in VSCode or equivalent.
 *            The compile instructions are in the JSON file 'tasks.JSON'. 
 *            One can change the function 'VectorXf F' in the file (functions.cpp),
 *            to solve for the unknown  variables in the system of nonlinear equations.
 * 
 *            Then hit Ctrl + B to build this program.
 *            Open a bash terminal and enter .\main.
 * 
 *            NOTE: This program uses third-party packages such such as Eigen. One need to use 
 *            LINUX to install the required packages.
 * 
 *  @file     main.cpp
 *  @author   Quek Jing Hao
 *  @version  August 2022
 *  
 * 
*/

using namespace Eigen;


//  driver program
int main()
{   
    // data declarations
    int row{4};
    int cols{4};
    double err{0.00000001};

    // matrix declarations
    MatrixXf mat(row, cols), trymat, b2;
    MatrixXf b(3, 1);
    Matrix3f A;
    Matrix3f G; 
    VectorXf x(3);
    VectorXf X(2); 
    Matrix3f B;

    //  system that needs to be solved
    A << 3.0, 4.0, 5.0, 
         1.0, 1.0, 1.0,
         1.0, 3.0, 7.0;

    X << 1.0, 1.0; 

    X = jfnk(&F, X, err, 20);
    
    std::cout << "By JFNK method, \n\n" << X << "\n\n";

    //   YES! WE DID IT!

    return 0;
}
