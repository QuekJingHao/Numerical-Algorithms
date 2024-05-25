# include <iostream>
# include <string>
# include <fstream>
# include <cmath>
# include <iomanip>
# include <chrono>
# include <utility>
# include <Eigen/Dense>
# include "Constants_U235.h"


/**
 * This function solves the discrete neutron diffusion equation in a fissioning weapon.
 * 
 * 
 * Purpose: The neutron diffusion equation is solved in a spherical coordinate (i.e in 1 dimension).
 * 
 * 			Later examples, we can try to solve it in a two dimensions - N(phi, theta, t)
 * 
 * 
 * 			Okay, let us not be so ambitious. Lets just solve the equation in 1D, and plot the results.
 * 			Once we think we are good enough, then we will make it into an animation
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * **/


using namespace Eigen;

// Driver program
int main() {


	// Declarations for FTCS scheme
	double R = 1;              // Radius of the bomb core in meters
	int Nr = 100;                    // Radial grid size
	int Nt = 100;                              // Number of iterations

	double h = R / static_cast<double>(Nr);     // Radial grid size

	// The time step tau, is related to the radial size via a equation
	double r_min = h;
	double tau = (r_min / (4.0 * beta)) * sqrt((4.0 * beta) / alpha);


	double T = tau * Nt; // Total time of simulation

	std::cout << "Choice of tau : " << tau << '\n';




	// Delcarations for solution vectors and matrices
	VectorXf N         = VectorXf::Ones(Nr);           // Solution vector i.e. neutron number density at time t
	VectorXf N_initial = VectorXf::Ones(Nr);           // Initial conditions in the bomb core



	// Set initial conditions in the bomb core - use a decaying exponential funtion
	for (int i = 1; i < Nr - 1; i++) {

		// Assume that at the center of the bomb, there is a initial N0 value
		float N0 = 5e10;
		float rate = 10000.0; 

		N_initial(i) = N0 * exp(-rate * i * h);

	}


	/*-----------------------------------------------------------------------*/
	/*      FTCS Scheme                                                      */
	/*-----------------------------------------------------------------------*/

	/*
	
	
	
	// Before the first iteration, we set the solution vector with initial conditions
	N = N_initial;
	for (int j = 1; j < Nt - 1; j++) {

		// Declare another vector to store the updated values for each iteration
		auto N_update = N;

		for (int i = 1; i < Nr - 1; i++) {

			double c1 = 2.0 * tau * (beta / (i * pow(h, 2)) + beta / pow(h, 2));
			double c2 = 2.0 * tau * (alpha - ((2.0 * beta) / pow(h, 2)));
			double c3 = 2.0 * tau * (beta / pow(h, 2) - beta / (i * pow(h, 2)));

			N_update(i) = N(i) + c1 * N(i + 1) + c2 * N(i) + c3 * N(i - 1);

		}

		// After each iteration, set the solution vector with the updated values
		N = N_update;

	}	*/

	/*
	// Print out the final values of the solution vector
	for (auto element : N) {

		std::cout << element << '\n';
	}
	
	*/

	return 0;

}