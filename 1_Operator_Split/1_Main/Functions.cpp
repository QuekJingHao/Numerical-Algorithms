# include <iostream>
# include <fstream>
# include <string>
# include <fftw3.h>
# include <complex>
# include <eigen3/Eigen/Dense>
# include <cmath>
# include "Headers/functions.h"
# include "Headers/mydefs.h"
# include "Headers/dft.h"


/******************************************************************************
 * 		                      Position operator
 * Returns the position operator evaulated at position x at a single time step
*******************************************************************************/
std::complex<double> Position_Operator(double x, double time_step)
{	
	std::complex<double> exponent(0, pow(x, 2.0)/2.0*time_step/2.0);
	return exp(-exponent);
}


/******************************************************************************
 * 		              Position operator - vector version
 * Returns the position operator evaulated at position x at a single time step
*******************************************************************************/
Eigen::Matrix<std::complex<double>, N, 1> Position_Operator(Eigen::Matrix<double, N, 1> x, double time_step)
{	
	Eigen::Matrix<std::complex<double>, N, 1> Ans;

	for (int i{0}; i<N; i++)
	{	
		std::complex<double> Temp (0.0, pow(x(i), 2.0)/2.0*time_step/2.0);
		Ans(i) = exp(-Temp);
	}

	return Ans;
}


/******************************************************************************
 * 		                      Momentum operator
 * Returns the momentum operator evaulated at momentum p at a single time step
*******************************************************************************/

std::complex<double> Momemtum_Operator(double y, double time_step)
{
	std::complex<double> exponent(0, pow(y, 2.0)/2.0*time_step);
	return exp(-exponent);
}


/******************************************************************************
 * 		              Momentum operator - vector version
 * Returns the momentum operator evaulated at momentum p at a single time step
*******************************************************************************/
Eigen::Matrix<std::complex<double>, N, 1> Momemtum_Operator(Eigen::Matrix<double, N, 1> y, double time_step)
{	
	Eigen::Matrix<std::complex<double>, N, 1> Ans;

	for (int i{0}; i<N; i++)
	{	
		std::complex<double> Temp (0.0, pow(y(i), 2.0)*time_step/2.0);
		Ans(i) = exp(-Temp);
	}
	return Ans;
}


								/*********************************
								 * Harmonic Oscillator Eigenstates
								**********************************/

//	Function for Harmonic Oscillator eigenstates
Eigen::Matrix<double, N, 1> HO_Eigenstates(Eigen::Matrix<double, N, 1> X, int n, double x0, double Sigma)
{	
	Eigen::Matrix<double, N, 1> return_array;

	for (int i{0}; i<N; i++)
	{
		double exponent = pow((X(i)-x0), 2.0)/(2.0*Sigma*Sigma);
		return_array(i) = std::hermite(n, X(i))*exp(-exponent);
	}

	double Sum = sqrt(return_array.dot(return_array));
	return (1.00/Sum)*return_array;

} 


							      /*********************************
							  	   * Perturbing Potential Functions
								  **********************************/

//	First perturbing potential
Eigen::Matrix<double, N, 1> V1(Eigen::Matrix<double, N, 1> X, double Amplitude, double Frequency, double t)
{
	Eigen::Matrix<double, N, 1> return_array;

	for (int i{0}; i<N; i++)
	{
		return_array(i) = Amplitude*sin(X(i))*cos(Frequency*t);
	}
	
	return return_array;
}


//	First perturbing potential - returns only 1 double 
double V1(double X, double Amplitude, double Frequency, double t)
{
	return Amplitude*sin(X)*cos(Frequency*t);
}


//	Second perturbing potential
Eigen::Matrix<double, N, 1> V2(Eigen::Matrix<double, N, 1> X, double Amplitude, double Frequency, double power, double t)
{
	Eigen::Matrix<double, N, 1> return_array;

	for (int i{0}; i<N; i++)
	{
		double term1 = pow((Amplitude/2.0), power);
		double term2 = pow(((X(i) - Amplitude)/2.0), power);

		return_array(i) = term1 - term2;
	}

	return (exp(-Frequency*t))*t*return_array;
}



// Position operator with perturbing potential function inside
Eigen::Matrix<std::complex<double>, N, 1> Position_Operator_Perturbed(Eigen::Matrix<double, N, 1> x, double Amplitude, double Frequency, double m, double time_step)
{	
	Eigen::Matrix<std::complex<double>, N, 1> Ans;

	for (int i{0}; i<N; i++)
	{	
		double VPerturb_Temp = V1(x(i), Amplitude, Frequency, m*time_step);

		std::complex<double> exponent (0.0, (pow(x(i), 2.0)/2.0 + VPerturb_Temp)*time_step/2.0);

		Ans(i) = exp(-exponent);
	}

	return Ans;
}



							/*************************************************
							 * Theoretical first order transition probability
							**************************************************/

Eigen::Matrix<double, M, 1> P_FirstOrder_Theory(double w, double wnm, Eigen::Matrix<double, M, 1> Time, double Vnm)
{
	Eigen::Matrix<double, M, 1> return_array;

	for (int i{0}; i<M; i++)
	{	
		std::complex<double> Plus_Exponent (0.0, Time(i)*(wnm + w));
		std::complex<double> Minus_Exponent (0.0, Time(i)*(wnm - w));

		std::complex<double> Plus = (exp(Plus_Exponent)-1.00)/(std::complex<double>(0.0, (wnm + w)));
		std::complex<double> Minus = (exp(Minus_Exponent)-1.00)/(std::complex<double>(0.0, (wnm - w)));

		std::complex<double> Inner_Sum = Plus + Minus;

		return_array(i) = abs(Inner_Sum*conj(Inner_Sum));

	}

	return abs(Vnm)*abs(Vnm)*return_array;
}



							/*************************************************
							 * 			  State Evolution Function
							**************************************************/

Eigen::Matrix<double, M+N, 1> State_Evolution(Eigen::Matrix<double, N, 1> X, Eigen::Matrix<double, N, 1> P, int n_i, int n_f, double Amplitude, double Frequency, double dt)
{		
	//	Declare empty array to store the transition probability
	Eigen::Matrix<double, M, 1> Transition_Probability;
	Eigen::Matrix<double, N, 1> Psi_Evolved_Probability;

	//	Declare array to return the two concatenated probabilities
	Eigen::Matrix<double, M+N, 1> Return_Array(Transition_Probability.size() + Psi_Evolved_Probability.size());

	//	Declare the UV and UT operators
	Eigen::Matrix<std::complex<double>, N, 1> UV;
	Eigen::Matrix<std::complex<double>, N, 1> UT = Momemtum_Operator(P, dt);


	//	Declare the theoretical initial and final Psi
	Eigen::Matrix<std::complex<double>, N, 1> Psi_Initial;
	Eigen::Matrix<std::complex<double>, N, 1> Psi_Final;

	//	Initialize the initial and final states with the eigenstates of the harmonic oscillator
	for (int i{0}; i<N; i++)
	{
		std::complex<double> Psi_Initial_Complex (HO_Eigenstates(X, n_i, X0, SD)(i), 0.0);
		Psi_Initial(i) = Psi_Initial_Complex;

		std::complex<double> Psi_Final_Complex (HO_Eigenstates(X, n_f, X0, SD)(i), 0.0);
		Psi_Final(i) = Psi_Final_Complex;

	}


	//	Declare the final evolved state
	Eigen::Matrix<std::complex<double>, N, 1> Psi_Evo_Final = Psi_Initial;


	//	Main part - evolve the state using fft
	for (int m{0}; m<M; m++)
	{	

		//	Declare UV matrix with the perturbing potential inside
		UV = Position_Operator_Perturbed(X, Amplitude, Frequency, static_cast<double>(m+1.00), dt);
		//	Thanks to MATLAB the bloody array starts at 1!!!!!

		Eigen::Matrix<std::complex<double>, N, 1> Psi_1 = UV.array()*Psi_Evo_Final.array();
		Eigen::Matrix<std::complex<double>, N, 1> Psi_2 = fft(Psi_1);
		Eigen::Matrix<std::complex<double>, N, 1> Psi_3 = UT.array()*Psi_2.array();
		Eigen::Matrix<std::complex<double>, N, 1> Psi_4 = ifft(Psi_3);
		Psi_Evo_Final = UV.array()*Psi_4.array();
 
		//	After every iteration, store the output transition probability as the output array
		std::complex<double> Temp_Sum = Psi_Final.dot(Psi_Evo_Final);
		Transition_Probability(m) = abs(Temp_Sum)*abs(Temp_Sum);
	}

	for (int i{0}; i<N; i++)
	{
		Psi_Evolved_Probability(i) = abs(Psi_Evo_Final(i))*abs(Psi_Evo_Final(i));
	}

	Return_Array << Transition_Probability, Psi_Evolved_Probability;

	return Return_Array;	// YES IT WORKSSSSSSSS!!!

}








							/*************************************************
							 * 		 Function to write results to file
							**************************************************/

//	Export data - each input is of size M 
void Export_Data_M(const std::string out_file_name, Eigen::Matrix<double, M, 1> X_Data, Eigen::Matrix<double, M, 1> Y_Data)
{
	//	Declare output file stream
	std::ofstream file_out{out_file_name};

	file_out << "X," << "Y" << '\n';

	for (int i{0}; i<M; i++)
	{
		file_out << X_Data(i) << ',' << Y_Data(i) << '\n';
	}

}

//	Export data - size N version 
void Export_Data_N(const std::string out_file_name, Eigen::Matrix<double, N, 1> X_Data, Eigen::Matrix<double, N, 1> Y_Data)
{
	//	Declare output file stream
	std::ofstream file_out{out_file_name};

	file_out << "X," << "Y" << '\n';

	for (int i{0}; i<N; i++)
	{
		file_out << X_Data(i) << ',' << Y_Data(i) << '\n';
	}

}








							/*************************************************
							 * Some other auxiliary functions not important
							**************************************************/

Eigen::Matrix<double, N, 1> HermitePol(Eigen::Matrix<double, N, 1> x, int order)
{	
	Eigen::Matrix<double, N, 1> Ans;

	for (int i{0}; i<N; i++)
	{
		Ans(i) = std::hermite(order, x(i));
	}

	return Ans;
}



//	Test function that returns two arrays in 1
Eigen::Matrix<double, M+N, 1> Concatenate(Eigen::Matrix<double, M, 1> X, Eigen::Matrix<double, N, 1> Y)
{
	Eigen::Matrix<double, M+N, 1> vec_joined(X.size() + Y.size());
	vec_joined << X, Y;

	return vec_joined;
}


//	Test function for addition
double add(double x, double y)
{
	return x+y;
}

//	Overloaded addition function that returns a vector
Eigen::Matrix<double, N, 1> add(Eigen::Matrix<double, N, 1> x, Eigen::Matrix<double, N, 1> y)
{
	return x + y;
} 





















