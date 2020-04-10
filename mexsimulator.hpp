#pragma once


#include "magnetization.hpp"

__global__ void mexsimulator(magnetization* magn, double* mxout, double* myout, 
        double* mzout, size_t* nelements,  size_t* numCols, size_t* numRows, 
        size_t* numPages,size_t* ndims, double* Gx, double* Gy, double* Gz,
        double* rfamp, double* rfphase, double* events, size_t* nEvents); 
