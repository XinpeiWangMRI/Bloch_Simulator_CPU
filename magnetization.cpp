#include "magnetization.hpp"
#include "math.h"
#include "mex.h"
#include "constants.hpp"

inline double sinc(double x) {
    if (x == 0) return 1.0;

    return sin(x) / x;
}

//Default Constructor, a lone voxel at the origin that is thermally relaxed with 1mm width in all dimensions. 
//Update so I am not using magic constants.
magnetization::magnetization() : mx(0.0),my(0.0),mz(1.0),xpos(0.0),ypos(0.0),zpos(0.0),bin(0),
                    offres(0.0), volume(1), voxelWidthX(.1), voxelWidthY(.1), voxelWidthZ(.1){}

void magnetization::rotate(const double bx, const double by, const double bz, const double tstep) 
{
    double xprod[3], tempm[3];
    double dot, phi, weff;
    

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

    //Need to update dephasing, too! 
    //local field grad is already in Hz/distance.
    localKCoord[0] += 2.0 * M_PI * tstep * localFieldGrad[0] * voxelWidthX;
    localKCoord[1] += 2.0 * M_PI * tstep * localFieldGrad[1] * voxelWidthY;
    localKCoord[2] += 2.0 * M_PI * tstep * localFieldGrad[2] * voxelWidthZ;

}

void magnetization::rotate(const double bx, const double by, const double bz, const double tstep, const bool advanceK)
{
    double xprod[3], tempm[3];
    double dot, phi, weff;


    weff = sqrt(bx * bx + by * by + bz * bz);
    phi = -2 * M_PI * weff * tstep;

    if (weff != 0.0) {

        xprod[0] = (by * mz - bz * my) / weff;
        xprod[1] = (bz * mx - bx * mz) / weff;
        xprod[2] = (bx * my - by * mx) / weff;

        dot = (bx * mx + by * my + bz * mz) / weff;

        tempm[0] = cos(phi) * mx + sin(phi) * xprod[0] + (1 - cos(phi)) * dot * bx / weff;
        tempm[1] = cos(phi) * my + sin(phi) * xprod[1] + (1 - cos(phi)) * dot * by / weff;
        tempm[2] = cos(phi) * mz + sin(phi) * xprod[2] + (1 - cos(phi)) * dot * bz / weff;

        mx = tempm[0];
        my = tempm[1];
        mz = tempm[2];
    }

    //Need to update dephasing, too! 
    //local field grad is already in Hz/distance.
    if (advanceK) {
        localKCoord[0] += 2.0 * M_PI * tstep * localFieldGrad[0] * voxelWidthX;
        localKCoord[1] += 2.0 * M_PI * tstep * localFieldGrad[1] * voxelWidthY;
        localKCoord[2] += 2.0 * M_PI * tstep * localFieldGrad[2] * voxelWidthZ;
    }

}

//If event has applied gradients.
void magnetization::rotate(const double bx, const double by, const double bz, const double tstep,
    const double appliedGrads[3])
{
    double xprod[3], tempm[3];
    double dot, phi, weff;


    weff = sqrt(bx * bx + by * by + bz * bz);
    phi = -2 * M_PI * weff * tstep;

    if (weff != 0.0) {

        xprod[0] = (by * mz - bz * my) / weff;
        xprod[1] = (bz * mx - bx * mz) / weff;
        xprod[2] = (bx * my - by * mx) / weff;

        dot = (bx * mx + by * my + bz * mz) / weff;

        tempm[0] = cos(phi) * mx + sin(phi) * xprod[0] + (1 - cos(phi)) * dot * bx / weff;
        tempm[1] = cos(phi) * my + sin(phi) * xprod[1] + (1 - cos(phi)) * dot * by / weff;
        tempm[2] = cos(phi) * mz + sin(phi) * xprod[2] + (1 - cos(phi)) * dot * bz / weff;

        mx = tempm[0];
        my = tempm[1];
        mz = tempm[2];
    }

    //Need to update dephasing, too! 
    //local field grad is already in Hz/distance.
    localKCoord[0] += 2.0 * M_PI * tstep * (gauss2Hz * appliedGrads[0] + localFieldGrad[0]) * voxelWidthX;
    localKCoord[1] += 2.0 * M_PI * tstep * (gauss2Hz * appliedGrads[1] + localFieldGrad[1]) * voxelWidthY;
    localKCoord[2] += 2.0 * M_PI * tstep * (gauss2Hz * appliedGrads[2] + localFieldGrad[2]) * voxelWidthZ;

}

//If event has an RF pulse
void magnetization::rotate(const double bx, const double by, const double bz, const double tstep,
    const double appliedGrads[3], const bool advanceK)
{
    double xprod[3], tempm[3];
    double dot, phi, weff;


    weff = sqrt(bx * bx + by * by + bz * bz);
    phi = -2 * M_PI * weff * tstep;

    if (weff != 0.0) {

        xprod[0] = (by * mz - bz * my) / weff;
        xprod[1] = (bz * mx - bx * mz) / weff;
        xprod[2] = (bx * my - by * mx) / weff;

        dot = (bx * mx + by * my + bz * mz) / weff;

        tempm[0] = cos(phi) * mx + sin(phi) * xprod[0] + (1 - cos(phi)) * dot * bx / weff;
        tempm[1] = cos(phi) * my + sin(phi) * xprod[1] + (1 - cos(phi)) * dot * by / weff;
        tempm[2] = cos(phi) * mz + sin(phi) * xprod[2] + (1 - cos(phi)) * dot * bz / weff;

        mx = tempm[0];
        my = tempm[1];
        mz = tempm[2];
    }

    //Need to update dephasing, too! 
    //local field grad is already in Hz/distance.
    if (advanceK) {
        localKCoord[0] += 2.0 * M_PI * tstep * (gauss2Hz * appliedGrads[0] + localFieldGrad[0]) * voxelWidthX;
        localKCoord[1] += 2.0 * M_PI * tstep * (gauss2Hz * appliedGrads[1] + localFieldGrad[1]) * voxelWidthY;
        localKCoord[2] += 2.0 * M_PI * tstep * (gauss2Hz * appliedGrads[2] + localFieldGrad[2]) * voxelWidthZ;
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
    //resetting to thermal equilibrium erases spin history.
    localKCoord[0] = 0;
    localKCoord[1] = 0;
    localKCoord[2] = 0;
};

void magnetization::refocusM(){
    mx = -mx;
    mz = -mz;
    reverseDephase();
};

void magnetization::setobj(const double objmz){
    mz = objmz;
};

void magnetization::setOffset(const double offset){
    offres = offset;
};  

void magnetization::setVoxelWidths(const double xWidth, const double yWidth, const double zWidth) {
    voxelWidthX = xWidth;
    voxelWidthY = yWidth;
    voxelWidthZ = zWidth;
};

void magnetization::reverseDephase(){
    localKCoord[0] = -localKCoord[0];
    localKCoord[1] = -localKCoord[1];
    localKCoord[2] = -localKCoord[2];

}
    
void magnetization::setFieldGrad(double fieldX, double fieldY, double fieldZ) {
    localFieldGrad[0] = fieldX;
    localFieldGrad[1] = fieldY;
    localFieldGrad[2] = fieldZ;
}

double magnetization::getKcoord(int index) {
    
    switch (index) {
    case 0:
        return localKCoord[0];
        break;

    case 1:
        return localKCoord[1];
        break;
    case 2:
        return localKCoord[2];
        break;
    default:
        mexErrMsgIdAndTxt("MATLAB:magnetization:getKcoord",
            "index into k-coordinate array must be 0, 1, or 2.");
        break;
    }
}

//Helper functions
void acquire(double* mxout, double* myout, double* mzout,
    const size_t ndims, const size_t time, const double mx, const double my, const double mz,
    const size_t bin, const size_t volume, const double localKcoord[3]) {
    size_t outputBin = bin + time * volume;

    mxout[outputBin] = sinc(localKcoord[0]) * sinc(localKcoord[1]) * mx;
    myout[outputBin] = sinc(localKcoord[0]) * sinc(localKcoord[1]) * my;

    //mxout[outputBin] = mx;
    //myout[outputBin] = my;
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

//Check if magnetization saw a refocusing pulse based off longitudinal magnetization
bool is_Ref_Pulse(const double newMz, const double oldMz){
    
    if (oldMz >= 0 && oldMz >= refocusFraction * newMz) { return true; }
    
    else if (oldMz < 0 && oldMz <= refocusFraction * newMz) { return true; }
    
    return false;
    
}