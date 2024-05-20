# include <iostream>
# include <string>
# include <fstream>
# include <cmath>
# include <iomanip>
# include <chrono>
# include <utility>
# include <Eigen/Dense>


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


// Driver program
int main() {


	// Constants declarations for fundamental quantities
	double nu_neut = 0.0;      // Average neutron velocity
	double nu = 0;             // Mean number of neutrons generate per fission event


	// Derived quantities
	double lambda_f = 9.0; // Fission mean free path
	double lambda_t = 0.0;   // Transport mean free path


	// Matrix declarations



	std::cout << "Hello world!";

	return 0;

}