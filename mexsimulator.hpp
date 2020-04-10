#pragma once

#include <vector>
#include "magnetization.hpp"
#include "mex.h"

void mexsimulator(std::vector<magnetization>& magn, double* mxout, double* myout,
    double* mzout, const size_t nelements, const size_t numCols, const size_t numRows,
    const size_t numPages, const size_t ndims, const double* gradx, const double* grady, const double* gradz,
    const double* rfpulse, const double* rfphase, const double* events, const size_t nEvents);
