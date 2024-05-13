# include <fftw3.h>
# include <complex>
# include <eigen3/Eigen/Dense>
# include "Headers/mydefs.h"
# include "Headers/dft.h"

using namespace Eigen;

//	Write a function that executes fft - it returns a vector of complex numbers
Matrix<std::complex<double>, N, 1> fft(Matrix<std::complex<double>, N, 1> input)
{
	fftw_complex in[N], out[N];
	fftw_plan plan;

	Matrix<double, N, 1> results_real;
	Matrix<double, N, 1> results_img;
	Matrix<std::complex<double>, N, 1> fft_result;

	for (int i{0}; i<N; i++)
	{	
		in[i][0] = real(input[i]);
		in[i][1] = imag(input[i]);
	} 

	plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan); 

	for (int i{0}; i<N; i++)
	{	
		//	Make a temporary variable to store the complex numbers before putting them into the output array
		std::complex<double> mycomplex (out[i][0], out[i][1]);
		fft_result[i] = mycomplex;
	} 

	fftw_destroy_plan(plan); 
	return fft_result;
}

//	Write a function that executes ifft - the funtcion returns a vector of complex numbers
Matrix<std::complex<double>, N, 1> ifft(Matrix<std::complex<double>, N, 1> input)
{
	fftw_complex in[N], out[N];
	fftw_plan plan;

	Matrix<double, N, 1> results_real;
	Matrix<double, N, 1> results_img;
	Matrix<std::complex<double>, N, 1> fft_result;

	for (int i{0}; i<N; i++)
	{	
		in[i][0] = real(input[i]);
		in[i][1] = imag(input[i]);
	} 

	plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan); 

	for (int i{0}; i<N; i++)
	{	
		//	Remember to normalize!
		double N_d = static_cast<double>(N);
		std::complex<double> mycomplex ((1.0/N_d)*out[i][0], (1.0/N_d)*out[i][1]);
		fft_result[i] = mycomplex;
	}

	fftw_destroy_plan(plan); 
	return fft_result;

}

