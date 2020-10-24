#include "fft.h"

/**
 * This file contains the implementation of `fft_compute` and ancillary
 * functions.
 */

// Include for math functions and definition of PI
#include <math.h>
// Included to get access to `malloc` and `free`
#include <stdlib.h>

// Forward declaration/ of helper methods
void get_even(const complex* in, complex* out, const int n) {
	for(int i = 0; i < n / 2; i++) {
		// Transfer all the even indexed numbers to the output array
		out[i] = in[2 * i];
	}
}

void get_odd(const complex* in, complex* out, const int n) {
	for(int i = 0; i < n / 2; i++) {
		// Transfer all the odd indexed numbers to the output array
		out[i] = in[2 * i + 1];
	}
}

	
void fft_twiddle(const complex* in, complex* out, const int n,complex* twiddle, const int index) {
	if(n == 1) {
		out[0] = in[0];
	} else {
		const int half = n / 2;
		// First we declare and allocate arrays
		// Allocate enough room for half the input values
		complex* even = malloc(sizeof(complex) * half);
		complex* odd  = malloc(sizeof(complex) * half);
		complex* even_out = malloc(sizeof(complex) * half);
		complex* odd_out  = malloc(sizeof(complex) * half);
		// Extract even and odd indexed numbers using methods above
		get_even(in, even, n);
		get_odd(in, odd, n);
		// Recursively calculate the result for bottom and top half
		fft_twiddle(even, even_out, n / 2, twiddle,index*2);
		fft_twiddle(odd, odd_out, n / 2, twiddle, index*2);
			

		// Combine the output of the two previous recursions
		for(int i = 0; i < half; ++i) {
			const complex e = even_out[i];
			const complex o = odd_out[i];
			const complex w = twiddle[i*index];
			out[i]        = e + w * o;
			out[i + half] = e - w * o;
		}
		// Since we allocated room for variables we need to release
		// the memory!
		free(even);
		free(odd);
		free(even_out);
		free(odd_out);


	}
}


void fft_compute(const complex* in, complex* out, const int n) {
	// compute all twiddle factors
	complex* twiddle = malloc(sizeof(complex) *(n/2));
	for(int i=0; i<n/2; i++){ 
		twiddle[i] = cexp(0 - (2. * M_PI * i) / n * I);
			
	}
	const int index = 1; 
	fft_twiddle(in, out, n, twiddle, index);
	free(twiddle);
	}


