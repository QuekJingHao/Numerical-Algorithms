#ifndef FUNCTIONS_H
#define FUNCTIONS_H

# include "mydefs.h"


/************************************************************************
 *    Computes the one-step position operator.
 *    @param  x is a single position double value
 *    @param  time_step is a single time step value i.e. dt
 *    @return A single double value of the position operator
 ************************************************************************/
std::complex<double> Position_Operator(double x, double time_step);


/************************************************************************
 *    Computes the one-step position operator - vector version
 *    @param  X is the position array of size N 
 *    @param  time_step is a single time step value i.e. dt
 *    @return A size N array of complex numbers
 ************************************************************************/
Eigen::Matrix<std::complex<double>, N, 1> Position_Operator(Eigen::Matrix<std::complex<double>, N, 1> x, double time_step);


/************************************************************************
 *    Computes the one-step momentum operator.
 *    @param  y is a single momentum double value
 *    @param  time_step is a single time step value i.e. dt
 *    @return A single double value of the momentum operator
 ************************************************************************/
std::complex<double> Momemtum_Operator(double y, double time_step);


/************************************************************************
 *    Computes the one-step momentum operator - vector version
 *    @param  y is a single momentum double value
 *    @param  time_step is a single time step value i.e. dt
 *    @return A single double value of the momentum operator
 ************************************************************************/
Eigen::Matrix<std::complex<double>, N, 1> Momemtum_Operator(Eigen::Matrix<std::complex<double>, N, 1> y, double time_step);



/************************************************************************
 *    Computes the harmonic oscillator eigenstates
 *    @param  X is a size N array of position values
 *    @param  n is the order of the Hermite polynomials
 *    @param  x0 is the value of the center of the Gaussian distribution
 *    @param  sigma is the variance (standard deviation) of the 
 *            Gaussian distribution
 *    @return An array of size N double values
 ************************************************************************/
Eigen::Matrix<double, N, 1> HO_Eigenstates(Eigen::Matrix<double, N, 1> X, int n, double x0, double sigma);



/************************************************************************
 *    What is this declaration
 *    @param  
 *    @param  
 *    @return
 ************************************************************************/
//  Perturbating potential
Eigen::Matrix<double, N, 1> V1(Eigen::Matrix<double, N, 1> X, double Amplitude, double Frequency, double t);



/************************************************************************
 *    What is this declaration
 *    @param  
 *    @param  
 *    @return
 ************************************************************************/
//  Theoretical first order transition probability
Eigen::Matrix<double, M, 1> P_FirstOrder_Theory(double w, double wnm, Eigen::Matrix<double, M, 1> Time, double Vnm);




/************************************************************************
 *    What is this declaration
 *    @param  
 *    @param  
 *    @return
 ************************************************************************/
//  State evolution
Eigen::Matrix<double, M+N, 1> State_Evolution(Eigen::Matrix<double, N, 1> X, Eigen::Matrix<double, N, 1> P, int n_i, int n_f, double Amplitude, double Frequency, double dt);



/************************************************************************
 *    What is this declaration
 *    @param  
 *    @param  
 *    @return
 ************************************************************************/
//  Export data to csv file for Python
void Export_Data_M(const std::string out_file_name, Eigen::Matrix<double, M, 1> X_Data, Eigen::Matrix<double, M, 1> Y_Data);


/************************************************************************
 *    What is this declaration
 *    @param  
 *    @param  
 *    @return
 ************************************************************************/
void Export_Data_N(const std::string out_file_name, Eigen::Matrix<double, N, 1> X_Data, Eigen::Matrix<double, N, 1> Y_Data);


/************************************************************************
 *    What is this declaration
 *    @param  
 *    @param  
 *    @return
 ************************************************************************/
Eigen::Matrix<double, M+N, 1> Concatenate(Eigen::Matrix<double, M, 1> X, Eigen::Matrix<double, N, 1> Y);


/************************************************************************
 *    What is this declaration
 *    @param  
 *    @param  
 *    @return
 ************************************************************************/
Eigen::Matrix<double, N, 1> HermitePol(Eigen::Matrix<double, N, 1> x, int order);


//  Misc functions
double add(double x, double y);
Eigen::Matrix<double, N, 1> add(Eigen::Matrix<double, N, 1> x, Eigen::Matrix<double, N, 1> y);

#endif