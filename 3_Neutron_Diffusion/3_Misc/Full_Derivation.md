# Numerical Algorithms

## Solving the Neutron Diffusion Equation in a Fissioning Nuclear Weapon

**<u>_Objective:_</u>** In this project, we solve the neutron diffusion equation of a fissioning nuclear weapon numerically, using the Foward-Time Centered-Space (FTCS) scheme. We will solve for the neutron field in spherical coordinates.


This document will detail all of the steps used to derive the neccessary equations needed in the project.


## Introduction - the Neutron Diffusion Equation

The Neutron Diffusion equation in a fissioning nuclear weapon, as introduced in Eqution (2.18) in _The Physics of the Manhattan Project_, is given by


$$
\boxed{\frac{\partial N(\vec r, t)}{\partial t} = \frac{\nu_{neut}}{\lambda_f}N(\vec r, t) \ + \ \frac{\lambda_t \nu_{neut}}{3} \left( \nabla^2 N(\vec r, t) \right),}
$$

where $N(\vec r, t)$ represents the neutron field $N$ at a position vector $\vec r$ at time $t$. The constants that appear in the equation are defined as:


- $\nu_{neut}$: Average neutron velocity
- $\lambda_{f}$: Fission mean free path
- $\nu$: Mean number of neutrons produced per fission event
- $\lambda_{t}$: Transport mean free path


Some of these quantities are derived from:

$$
\lambda_f = \frac{1}{\sigma_f n}
$$

$$
\lambda_t = \frac{1}{\sigma_t n}
$$

Where the mean free paths are expressed in terms of their cross-sections $\sigma$. The transport cross section is in turn, the sum of the fission and elastic scattering cross sections: 

$$
\sigma_{t} = \sigma_{f} + \sigma_{el}
$$

Lastly, $n$ is a constant, the number density