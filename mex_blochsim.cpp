/*==========================================================
 * mex_blochsim.cpp - C++ Bloch Simulator
 *
 *
 *
 *========================================================*/

#include <iostream>
#include "mex.h"
#include "math.h"
#include "magnetization.hpp"
#include "mexsimulator.hpp"
#include "event_manager.hpp"
#include <string>
#include <float.h>
#include <vector>

#define mymax(a,b) a > b ? a : b

 /* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {

	// compile with:
	//    mex mex_blochsim.cpp mexsimulator.cpp magnetization.cpp event_manager.cpp;
		/******************* CHECK NUMBER OF ARGUMENTS **********************/
	size_t usrDefNRHS = 1;
	size_t usrDefNLHS = 3;

	std::string RHSerrormsg; //strings to store error messages
	std::string LHSerrormsg;
	std::string sNRHS = std::to_string(usrDefNRHS); //add number of arguments to string
	std::string sNLHS = std::to_string(usrDefNLHS);
	RHSerrormsg = "MEX_BLOCHSIM requires " + sNRHS + " input arguments."; //concat. strings
	LHSerrormsg = "MEX_BLOCHSIM requires " + sNLHS + " output arguments.";
	if (nrhs != usrDefNRHS) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
			RHSerrormsg.c_str()); //.c_str returns pointer to first element.
	}
	if (nlhs != usrDefNLHS) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
			LHSerrormsg.c_str());
	}

	/************************** ALTERNATE METHOD OF READING IN ELEMENTS USING A STRUCT ************************/
	size_t nfields = mxGetNumberOfFields(prhs[0]);
	size_t ifield;
	char* fname = NULL;
	fname = new char[256]; //don't have fieldnames longer than 256 characters

	//Spatial grids
	double* xgrid;
	double* ygrid;
	double* zgrid;

	/* Acquire pointers to the input data */
	size_t numRows, numCols, numPages, nelements;

	//set some default values
	numRows = 1;
	numPages = 1;
	numCols = 1;

	//Why is it letting me do this? The const without initialization...
	const size_t* dims;
	size_t ndims;

	double* Gx;
	double* Gy;
	double* Gz;
	size_t nGx;
	size_t nGy;
	size_t nGz;

	double* rfamp;
	size_t nrfamp;
	double* rfphase;
	size_t nrfphase;

	double* events;
	size_t numelemEvents;
	size_t nEvents;

	//Get user-defined object. Need an error check on the size of it still 
	double* usrObj;

	//Get voxel widths. Either need a map, a single global value, or a triplet of values (for x, y, z widths).
	double* voxelWidthX;
	double* voxelWidthY;
	double* voxelWidthZ;

	//need to include voxel widths. Account for either a single universal value, a triplet of values
	//(uniform in space but different in each dimension), or a triplet of values for each point in space.    
	bool tripleWidth = false, globalWidth = false;
	//default voxel widths of 1 mm
	double localWidthX = .1, localWidthY = .1, localWidthZ = .1;
	size_t ndims_VoxelWidths;

	for (int ifield = 0; ifield < nfields; ifield++) {
		fname = (char*)mxGetFieldNameByNumber(prhs[0], ifield);

		if (strcmp(fname, "xgrid") == 0) {
			xgrid = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
			dims = mxGetDimensions(mxGetFieldByNumber(prhs[0], 0, ifield));
			ndims = mxGetNumberOfDimensions(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "ygrid") == 0) {
			ygrid = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "zgrid") == 0) {
			zgrid = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "Gx") == 0) {
			Gx = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
			nGx = mxGetNumberOfElements(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "Gy") == 0) {
			Gy = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
			nGy = mxGetNumberOfElements(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "Gz") == 0) {
			Gz = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
			nGz = mxGetNumberOfElements(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "rfamp") == 0) {
			rfamp = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
			nrfamp = mxGetNumberOfElements(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "rfphase") == 0) {
			rfphase = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
			nrfphase = mxGetNumberOfElements(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "events") == 0) {
			events = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
			numelemEvents = mxGetNumberOfElements(mxGetFieldByNumber(prhs[0], 0, ifield));
			nEvents = mxGetM(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "usrObj") == 0) {
			usrObj = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
		}
		else if (strcmp(fname, "VoxelWidths") == 0) {
			ndims_VoxelWidths = mxGetNumberOfDimensions(mxGetFieldByNumber(prhs[0], 0, ifield));
			const size_t* dims_VoxelWidths = mxGetDimensions(mxGetFieldByNumber(prhs[0], 0, ifield));

			switch (ndims_VoxelWidths) {
				//If only one value is specified, assign it to all dimensions.
			case 1:
				localWidthX = mxGetScalar(mxGetFieldByNumber(prhs[0], 0, ifield));
				localWidthY = localWidthX;
				localWidthZ = localWidthX;
				break;
				//if there are 2 dimensions, one of them needs to be 3
			case 2: {
				if (dims_VoxelWidths[0] != 3 || dims_VoxelWidths[1] != 3) {
					mexErrMsgIdAndTxt("MATLAB:gateway:voxelDimensions",
						"Set of voxel widths must be either scalar, 3x1 (or 1x3) vector, or same size as spatial grids");
				}
				double* allWidths = mxGetPr(mxGetFieldByNumber(prhs[0], 0, ifield));
				localWidthX = allWidths[0];
				localWidthY = allWidths[1];
				localWidthZ = allWidths[2];
				break;
			}
			//Still need to write this interface.
			case 3: {
				break;
			}
			default:
				mexErrMsgIdAndTxt("MATLAB:gateway:voxelDimensions",
					"Set of voxel widths must be either scalar, 3x1 (or 1x3) vector, or same size as spatial grids");

			}
		}
	}

	/* This is the number of magnetization vectors being simulated */
	nelements = 1;
	for (int index = 0; index < ndims; index++) {
		nelements = dims[index] * nelements;
	}

	if (nrfamp - nrfphase != 0) {
		mexErrMsgIdAndTxt("MATLAB:mex_blochsim:typeargin",
			"RF amp and phase must be same size.");
	}
	/* Get total number of points in all acquisitions of experiment so can
	 * allocate memory appropriately */

	size_t nTime = 0;

	int dummyindex;
	for (int dummy = 0; dummy < static_cast<int>(nEvents); dummy++) {
		if (events[dummy] == 4 || events[dummy] == 5) { //these are acquisition events
			dummyindex = dummy + 2 * static_cast<int>(nEvents); //accesses number of points in event
			nTime += (size_t)events[dummyindex];
		}
	}

	size_t* dimsOut = new size_t[ndims + 1];
	switch (ndims) {
	case 1:
		numCols = mymax(dims[0], 1);
		dimsOut[0] = numCols;
		dimsOut[1] = nTime;
		break;
	case 2:
		numCols = mymax(dims[0], 1);
		numRows = mymax(dims[1], 1);
		dimsOut[0] = numCols;
		dimsOut[1] = numRows;
		dimsOut[2] = nTime;

		break;
	case 3:
		numCols = mymax(dims[0], 1);
		numRows = mymax(dims[1], 1);
		numPages = mymax(dims[2], 1);
		dimsOut[0] = numCols;
		dimsOut[1] = numRows;
		dimsOut[2] = numPages;
		dimsOut[3] = nTime;

		break;
	default:
		mexErrMsgIdAndTxt("MATLAB:mex_blochsim:ndims",
			"Number of dimensions of xgrid must be 1, 2, or 3.");
	};


	/*Initialize magnetization array*/
	std::vector<magnetization> magn(nelements);

	for (int index = 0; index < nelements; index++) {
		magn[index].setBin(index, numCols, numRows, numPages);
		magn[index].setVolume(numCols, numRows, numPages);
		magn[index].setpos(index, xgrid, ygrid, zgrid);
		magn[index].setobj(usrObj[index]);
		//not ready yet.
		/*
		switch ndims_VoxelWidths{
		* case 3:
		*	magn[index].setVoxelWidth(voxelWidthsX[index],voxelWidthsY[index],voxelWidthsZ[index]);
		*	break;
		*
		* default:
		*	magn[index].setVoxelWidth(localWidthX,localWidthY,localWidthZ);
		*	break;
		*}
		*/
	}

	plhs[0] = mxCreateNumericArray(ndims + 1, dimsOut, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(ndims + 1, dimsOut, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericArray(ndims + 1, dimsOut, mxDOUBLE_CLASS, mxREAL);

	double* mxout = mxGetPr(plhs[0]);
	double* myout = mxGetPr(plhs[1]);
	double* mzout = mxGetPr(plhs[2]);

	mexsimulator(magn, mxout, myout, mzout, nelements, numCols, numRows,
		numPages, ndims, Gx, Gy, Gz, rfamp, rfphase, events, nEvents);

	delete[] dimsOut;
	delete[] fname;

	dimsOut = NULL;
	fname = NULL;

	mexPrintf("Simulation Complete! \n");
}

