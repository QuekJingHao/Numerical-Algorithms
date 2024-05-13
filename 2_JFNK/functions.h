#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using namespace Eigen;

/************************************************************************
 * 
 *    F is the system of nonlinear equations in residual form
 *    As F is dynamically sized, in general, you can have n variables.
 *    I.e. F can be 
 * 
 *    VectorXf F(double x, double y, double z, double u, ... double n);
 * 
 *    @param  x unknown variable to be solved for 
 *    @param  y unknown variable to be solved for 
 *    @return   vector of size number of unknown variables evaluated at 
 *            (x, y, z, ...)
 * 
 ************************************************************************/
VectorXf F(double x, double y);



/****************************************************************************
 * 
 *    Approximate the Jacobian matrix using finite difference method
 * 
 *    @param  X    input vector for which the Jacobian will be evaluated at
 *    @param  func function name of the nonlinear system arranged in 
 *                 residual form
 *    @return      a matrix of size (n by n) where n is the number of unknown 
 *                 variables
 * 
 ****************************************************************************/
MatrixXf jacobian_approx(VectorXf X, std::function<VectorXf(double, double)>func);




/****************************************************************************
 * 
 *    Performs the Arnoldi iteration on the system Ax = b
 * 
 *    @param  A   input matrix of size (m by m)
 *    @param  b   input vector of size (m by 1)
 *    @param  n   number of iteration steps you want to perform Arnoldi
 *    @return     a tuple of q and h matrices, where q of size (m by n+1)
 *                is the Arnoldi matrix and h of size (n+1 by n)
 *                is the Hessenberg matrix
 *               
 ****************************************************************************/
std::pair<MatrixXf, MatrixXf> arnoldi_iteration(MatrixXf A, MatrixXf b, int n);




/****************************************************************************
 * 
 *    Perform the Givens rotation on input matrix
 * 
 *    @param  A     input matrix of size (m by m)
 *    @param  beta  input vector of size (m by 1)
 *    @param  k     number of iteration steps you want to apply Givens
 *    @return       a tuple of (A, beta)
 *               
 ***************************************************************************/
std::pair<MatrixXf, VectorXf> givens_rotation(MatrixXf A, VectorXf beta, int k);



/****************************************************************************
 * 
 *    Uses the Generalized Mininum Residuals (GMRES) method to solve for x
 *    in the system Ax = b.
 * 
 *    @param  A              input matrix of size (m by m)
 *    @param  b              input vector of size (m by 1)
 *    @param  x              input vector of the guessed solution 
 *                           of size (m by 1)
 *    @param  max_iterations maximum number of iterations by GMRES
 *    @param  error          error tolerance
 *    @return                solution vector of size (m by 1)
 *               
 ***************************************************************************/
VectorXf gmres(MatrixXf A, MatrixXf b, VectorXf x, int max_iterations, double error);





/****************************************************************************
 * 
 *    Performs the Jacobian-Free Newton Krylov method for solving system of
 *    nonlinear equations
 * 
 *    @param  func           function name of the nonlinear system arranged in 
 *                           residual form 
 *    @param  X_initial      initial guess of solution vector of size (m by 1)
 *    @param  max_iterations maximum number of iterations by GMRES
 *    @param  error          error tolerance
 *    @return                solution vector of size (m by 1)
 *               
 ***************************************************************************/
VectorXf jfnk(std::function<VectorXf(double, double)>func, VectorXf X_initial, double err, int max_iterations);


#endif