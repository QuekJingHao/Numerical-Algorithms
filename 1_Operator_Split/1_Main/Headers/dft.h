#ifndef DFT_H
#define DFT_H

# include <fftw3.h>
# include <complex>
# include <eigen3/Eigen/Dense>
# include "mydefs.h"

using namespace Eigen;

Matrix<std::complex<double>, N, 1> fft(Matrix<std::complex<double>, N, 1> input);

Matrix<std::complex<double>, N, 1> ifft(Matrix<std::complex<double>, N, 1> input);

#endif