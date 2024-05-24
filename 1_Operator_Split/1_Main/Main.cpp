# include <iostream>
# include <string>
# include <fstream>
# include <cmath>
# include <iomanip>
# include <complex>
# include <chrono>
# include <eigen3/Eigen/Dense>
# include <fftw3.h>
# include "Headers/mydefs.h"
# include "Headers/functions.h"
# include "Headers/dft.h"

/*
 * 	PC4230 Quantum Mechanics III: Full C++ implementation of the QM3 project code 
 * 	to study evolution of HO eigenstate
 * 
 * 	Purpose: 	 Studies evolution of eigenstates of harmonic oscillators, under a 
 * 			 	 time-dependent perturbing potential. We will use the operator split method 
 * 			 	 to advance the evolution of the HO eigenstate at each iteration. 
 * 
 * 				 The results of the calculation is exported as a csv file to be read by a Python
 * 			     front-end script. 	
 * 	
 * @file         Main.cpp
 * @author       Quek Jing Hao
 * @version      May 2022
 * 
*/

using namespace Eigen;

int main()
{	
	//	Type alias for some complicated data types
	using complex_t = std::complex<double>;
	using Matrix_t = Matrix<double, N, 1>;

	//	Data declarations
	auto a{-20.0}; 								   // Left end point of the trap
	auto b{20.0}; 							  	   // Right end point of the trap
	auto L{b-a}; 							  	   // Width of the trap
	auto pi{3.14159265358};						   // Constant pi
	auto T{125.0*pi};							   // Time duration of the evolution     

	auto dt{T/static_cast<double>(M)};			   // Single time step   

	Matrix_t X;					 	 			   // Dimensionless spatial coordinates
	Matrix_t P;						 			   // Dimensionless momentum lattice
	Matrix<double, M, 1> Time;					   // Time vector, to plot the transition probabilities against time
	
	Matrix<double, M, 1> Transition_Probability;   // Array to store the first order transition probabilties
	Matrix<double, N, 1> Psi_Evolved_Probability;  // Array to store probability of Psi after evolution
	Matrix<double, M+N, 1> State_Evolved_Data;	   // Array to store the entire data from the State Evolution Function

	auto start = std::chrono::steady_clock::now();

	//	Initialize the X and P arrays
	for (int i{0}; i<N; i++)
	{
		X(i) = a + L*static_cast<double>(i)/static_cast<double>(N);
		P(i) = (2*pi/L)*static_cast<double>(i);

		if (i>=N/2)
		{
			P(i) = -(2*pi/L)*(static_cast<double>(N)-static_cast<double>(i));
		}

	}

	//	Initialize Time and Transition probability arrays. For the last one set everything as 0
	for (int i{0}; i<M; i++)
	{
		Time(i) = dt*(static_cast<double>(i + 1.0));
	}

	std::cout << "Executing calculation..." << '\n';

	/*******************************************************************
	*	Use the state evolution function to evolve the initial state
	********************************************************************/
	State_Evolved_Data = State_Evolution(X, P, 0, 1, A, omega, dt);

	Transition_Probability = State_Evolved_Data.head(M);
	Psi_Evolved_Probability = State_Evolved_Data.tail<N>();


	/*******************************************************************
	*	Question 1: Investigating the spatial profile of wavefunction
	********************************************************************/

	Eigen::Matrix<std::complex<double>, N, 1> Psi_Initial;
	Matrix<double, N, 1> Psi_Initial_Probability;

	//	Initialize the initial eigenstate of the harmonic oscillator 
	//	- this is ripped from the function body lol
	for (int i{0}; i<N; i++)
	{
		std::complex<double> Psi_Initial_Complex (HO_Eigenstates(X, 0, X0, SD)(i), 0.0);
		Psi_Initial(i) = Psi_Initial_Complex;
		Psi_Initial_Probability(i) = abs(Psi_Initial(i))*abs(Psi_Initial(i));

	}
	

	/************************************************************************
	*	Question 3, 4: Investigating the transition probability against time
	*************************************************************************/
	// Some constants that will enter into the theoretical expression as Vnm

	double eta = 0.690194;
	double C = sqrt(2.0/pi)*eta*A*(1.0/2.0);


	/************************************************************************
	*     Question 6: Investigating second order transition probabilities
	*************************************************************************/
	//	First compute the inner product of V20:
	Matrix_t V_Sin;
	Matrix_t V_20;
	double C_20;
	
	for (int i{0}; i<N; i++)
	{
		V_Sin(i) = A*sin(X(i));
	}

	V_20 = V_Sin.array()*HO_Eigenstates(X, 0, X0, SD).array();

	// This will enter the theoretical first order transition probability function
	C_20 = HO_Eigenstates(X, 2, X0, SD).dot(V_20)/2.00;	


	//	Change the constant C if you are looking at second order transitions
	Matrix<double, M, 1> Theoretical_Transition_Probability = P_FirstOrder_Theory(omega, 1.00, Time, C);


	/************************************************************************
	*     			  Export data as csv files for Python
	*************************************************************************/

	//	For spatial profiles
	Export_Data_N("Psi Initial.csv", X, Psi_Initial_Probability);
	Export_Data_N("Psi Evolved.csv", X, Psi_Evolved_Probability); 

	//	For transition probabilities
	Export_Data_M("Simulated Transition Probability.csv", Time, Transition_Probability);
	Export_Data_M("Theoretical Transition Probability.csv", Time, Theoretical_Transition_Probability);


	std::cout << "Calculation completed." << '\n';
	std::cout << "\n ------------------------------------------** Exploitation Successful **--------------------------------------- \n";
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Total calculation time = " << elapsed_seconds.count() << " seconds. \n\n"; // Print out calculation time

	return 0;
}
