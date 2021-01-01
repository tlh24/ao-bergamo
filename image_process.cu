// process an image using CUDA. 
#include <iostream>
#include <stdio.h>
#define MAX_THREADS 256
#define SIZE (2048*2048)

unsigned char* d_image; 
float* d_results; 
int g_imagew; 
int g_imageh; 

void cu_setup(int width, int height, int nresults){
	cudaMalloc(&d_image, width*height); 
	cudaMalloc(&d_results, nresults); 
	if(d_image == 0 || d_results == 0){
		std::cout << "could not allocate CUDA memory" << endl; 
		imagew = imageh = 0
	} else {
		g_imagew = width; 
		g_imageh = height; 
	}
}

__global__ void sumall_kernel(unsigned char *dat, float *d_result){
    __shared__ float cache[MAX_THREADS];

    int cacheIdx = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    cache[cacheIdx] = i >= SIZE ? 0. : d_vector[i];
    __syncthreads();

    if (i >= SIZE)
        return;

    int padding = blockDim.x/2;
    while (padding != 0)
    {
        if (cacheIdx < padding)
            cache[cacheIdx] += cache[cacheIdx + padding];

        __syncthreads();
        padding /= 2;
    }

    if (cacheIdx == 0)
        atomicAdd(&d_result[0], cache[0]);
}

void cu_process(unsigned char* indata){
	
}
