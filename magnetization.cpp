#include "magnetization.hpp"
#include "math.h"
#include "mex.h"

//Default Constructor, a lone spin at the origin that is thermally relaxed.
magnetization::magnetization() : mx(0.0),my(0.0),mz(1.0),xpos(0.0),ypos(0.0),zpos(0.0),bin(0),
                    offres(0.0), volume(1){}

void magnetization::rotate(const double bx, const double by, const double bz, const double tstep) {
    double xprod[3], tempm[3];
    double dot, phi, weff;
    double M_PI = 4.0 * atan(1.0);

    weff = sqrt(bx*bx + by*by + bz*bz);
    phi = -2*M_PI*weff * tstep;
    
    if (weff != 0.0) {
        
        xprod[0] = (by*mz - bz*my)/weff;
        xprod[1] = (bz*mx - bx*mz)/weff;
        xprod[2] = (bx*my - by*mx)/weff;
        
        dot = (bx*mx + by*my + bz*mz)/weff;
        
        tempm[0] = cos(phi)*mx + sin(phi)*xprod[0] + (1-cos(phi))*dot*bx/weff;
        tempm[1] = cos(phi)*my + sin(phi)*xprod[1] + (1-cos(phi))*dot*by/weff;
        tempm[2] = cos(phi)*mz + sin(phi)*xprod[2] + (1-cos(phi))*dot*bz/weff;
        
        mx = tempm[0];
        my = tempm[1];
        mz = tempm[2];
    }
}


//Calculates the output bin based on index into nD arrary (n = 1, 2, 3)
void magnetization::setBin(const size_t index, const size_t numCols, const size_t numRows,
    const size_t numPages) {
   
    size_t z = (index + 1) % (numRows * numCols);
    size_t jj, ii, kk;
    
    if (z != 0) {
        jj = z % (numCols);
        if (jj == 0) {
            jj = numCols;
        }
        ii = ((z - jj)/(numCols)) + 1;
    }
    else{
        jj = numCols;
        ii = numRows;
    }
    ii = ii - 1;
    jj = jj - 1;
    kk = (index - numCols*(ii)-jj)/(numRows * numCols);
   
    bin = jj + (numCols)*ii + ((numRows)*(numCols))*kk;
};

void magnetization::setVolume(const size_t numRows, const size_t numCols, const size_t numPages){
    volume = numRows * numCols * numPages;
};

void magnetization::setpos(const size_t index,double* xgrid, double* ygrid, double* zgrid) {
    xpos = xgrid[index];
    ypos = ygrid[index];
    zpos = zgrid[index];
};

void magnetization::set2eq(){
    mx = 0;
    my = 0;
    mz = 1;
};

void magnetization::refocusM(){
    mx = -mx;
    mz = -mz;
};

void magnetization::setobj(const double objmz){
    mz = objmz;
};

void magnetization::setOffset(const double offset){
    offres = offset;
};  

//Helper functions
void acquire(double* mxout, double* myout, double* mzout,
    const size_t ndims, const size_t time, const double mx, const double my, const double mz,
    const size_t bin, const size_t volume) {
    size_t outputBin = bin + time * volume;
    mxout[outputBin] = mx;
    myout[outputBin] = my;
    mzout[outputBin] = mz;

};

//Check if magnetization is valid value. root sum of squares should be <= 1. By extension, sum of squares <= 1.
bool is_Valid_Magn(const magnetization* magn) {
    double total_Magn = magn->getMx() * magn->getMx();
    total_Magn += magn->getMy() * magn->getMy();
    total_Magn += magn->getMz() * magn->getMz();

    if (total_Magn > 1.0 || total_Magn < 0.0) {
        mexErrMsgIdAndTxt("MATLAB:magnetization:magnetizationValue",
            "Magnitude of magnetization must be between 0 and 1 inclusive");
        return false; //redundant
    }
    return true;
}