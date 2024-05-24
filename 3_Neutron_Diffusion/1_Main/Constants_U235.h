#ifndef CONSTANTS_U235_H
#define CONSTANTS_U235_H

/**
*	This file contains all of the constant used in the calculation
*
*	All of the values are obtained from the book
* 
*	These values are for untampered U-235 atom
* 
* 
* 
*
*
**/


// Constants declarations for fundamental quantities
double sigma_f   = 1.235e-24;                  // Fission cross-section / cm-2
double sigma_el  = 4.566e-24;                 // Elastic scattering cross-section / cm-2
double n         = 4.794e22;                 // Nuclear number density / cm-3

double nu_neut  = 2200e3;                   // Average neutron velocity (for thermal neutrons) / ms-1
double nu       = 2.637;                    // Mean number of neutrons generate per fission event



// Derived quantities
auto sigma_t   = sigma_f + sigma_el;
auto lambda_t  = 1.0 / (n * sigma_t) / 100;     // Transport mean-free path / m
auto lambda_f  = 1.0 / (n * sigma_f) / 100;     // Fission mean-free path / m


// Defined constants to simply equation
auto alpha = (nu_neut / lambda_f) * (nu - 1.0) ;  
auto beta = (lambda_t * nu_neut) / 3.0;           


#endif