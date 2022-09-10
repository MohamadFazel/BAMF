
//#include "cuda_runtime.h"
#include "definitions.h"
//#include "kernel.h"

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

__global__ void kernel_gaussX(float * d, int Xsize, float b0, float b1, float b2, float b3, float B)
{
    
	//this kernel does gaussian filter along the X dimension.  
	int Ysize = blockDim.x;
	int idy = threadIdx.x;
	int idz = blockIdx.x;

	float w0, w1, w2, w3;
	float temp;
	int ii = 0;
	const int base = idz*Xsize*Ysize + idy*Xsize;

	//forward
	w1 = w2 = w3 = d[base];
	for (ii = 0; ii<Xsize; ii++)
	{
		w0 = d[base + ii];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}

	//backward
	w1 = w2 = w3 = d[base + Xsize - 1];
	for (ii = Xsize - 1; ii >= 0; ii--)
	{
		w0 = d[base + ii];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}
}

__global__ void kernel_gaussY(float * d, int Ysize, float b0, float b1, float b2, float b3, float B)
{
	//this kernel does gaussian filter along the Y dimension.  
	int Xsize = blockDim.x;
	int idx = threadIdx.x;
	int idz = blockIdx.x;

	float w0, w1, w2, w3;
	float temp;
	int ii = 0;
	const int base = idz*Xsize*Ysize + idx;

	//forward
	w1 = w2 = w3 = d[base];
	for (ii = 0; ii<Ysize; ii++)
	{
		w0 = d[base + ii*Xsize];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii*Xsize] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}

	//backward
	w1 = w2 = w3 = d[base + Xsize*(Ysize - 1)];
	for (ii = Ysize - 1; ii >= 0; ii--)
	{
		w0 = d[base + ii*Xsize];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii*Xsize] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}
}

__global__ void kernel_subtract(float * d_A, float * d_B)
{
	//this kernel does gaussian filter along the Y dimension.  
	int Xsize = blockDim.x;
	int Ysize = gridDim.x;
	int idx = threadIdx.x + Xsize*blockIdx.x + Xsize*Ysize*blockIdx.y;
	d_A[idx] = d_A[idx] - d_B[idx];
}

__global__ void kernel_maxX(float * d_A, float * d_B, int kernelsz, float minval)
{
	//this kernel does max finding along the X dimension.  
	int Xsize = blockDim.x;
	int Ysize = gridDim.x;
	
	//this is the pixel that we are searching around
	int idx = threadIdx.x + Xsize*blockIdx.x + Xsize*Ysize*blockIdx.y;
	int x = threadIdx.x;

	//define search only up to edges 
	int start = fmaxf(0, x - kernelsz) + Xsize*blockIdx.x + Xsize*Ysize*blockIdx.y;
	int end = fminf(Xsize - 1, x + kernelsz) + Xsize*blockIdx.x + Xsize*Ysize*blockIdx.y;

	float maxval = minval;
	float inpixel = d_A[idx];

	for (int ii = start; ii<end + 1; ii++) 
		maxval = fmaxf(maxval, d_A[ii]);

	//if any other pixel is larger set pixel idx to negative of that value, otherwise keep
	d_B[idx] = (maxval>inpixel)*-maxval + (maxval == inpixel)*maxval;

}

__global__ void kernel_maxY(float * d_A, float * d_B, int kernelsz, float minval)
{
	//this kernel does max finding in second dimension.
	int Ysize = blockDim.x;
	int Xsize = gridDim.x;

	int x = blockIdx.x;
	int y = threadIdx.x;
	int z = blockIdx.y;

	//this is the pixel that we are searching around
	int idx = x + Xsize*y + Xsize*Ysize*z;

	//define search only up to edges 
	int start = fmaxf(0, y - kernelsz);
	int end = fminf(Ysize - 1, y + kernelsz);

	float maxval = minval;
	float inpixel = d_B[idx];

	//find the maximum absolute value in the filter window
	for (int ii = start; ii<=end; ii++)
		maxval = fmaxf(maxval, fabsf(d_B[ii*Xsize + x + Xsize*Ysize*z]));

	//if our pixel under test is equal to maximum, then flag that with '1', otherwise '0'
	d_A[idx] = fabsf(maxval - inpixel)< 1e-6;
	

}

