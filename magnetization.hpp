#include "mex.h"

#pragma once

class magnetization {
  public:
    //Default constructor
    magnetization();

    //Update mx, my, mz during events
    void rotate(const double bx, const double by, const double bz, const double tstep);

    //Set output spatial bin. Used to compute output based on time in acquire events.
    void setBin(const size_t index, const size_t numCols, const size_t numRows, const size_t numPages);
    //inline simple function which won't change when representation of magn changes.
    size_t getBin() const { return bin; };

    //Set spatial positions
    void setpos(const size_t index, double* xgrid, double* ygrid, double* zgrid);

    //Return spatial positions and spin components. Appropriate to inline since they would not likely change when representation of magnetization changes.
    double getX() const { return xpos; };
    double getY() const { return ypos; };
    double getZ() const { return zpos; };
    double getMx() const { return mx; };
    double getMy() const { return my; };
    double getMz() const { return mz; };

    //Set magnetization to equilibrium value
    void set2eq();

    //Set two magnetization components = -(themselves), this would be the case with a perfect refocusing pulse.
    void refocusM();

    //Set the equilibrium magnetization to be a different value.
    void setobj(double objmz);

    //Set/Get off-resonance frequency.
    void setOffset(const double offset);
    double getOffset() const { return offres; };

    //Set volume
    void magnetization::setVolume(const size_t numRows, const size_t numCols, const size_t numPages);
    size_t magnetization::setVolume() const { return volume; };

    //Set voxel widths
    void magnetization::setVoxelWidths(const double xWidth, const double yWidth, const double zWidth);

  private:
      //current spin components, spatial position, and offresonance frequency, respectively.
    double mx,my,mz,xpos,ypos,zpos,offres;
    double voxelWidthX, voxelWidthY, voxelWidthZ;

    //Bin is used for output
    size_t bin;
    //Volume is total number of spins, will look back to see if this is even needed now.
    size_t volume;
};

//Helper functions

//Check if magnetization is valid value. root sum of squares should be <= 1. By extension, sum of squares <= 1.
bool is_Valid_Magn(const double* magnetization);

//Update output arrays by inputting (part of) current magnetization state.
void acquire(double* mxout, double* myout, double* mzout, const size_t ndims,
    const size_t time, const double mx, const double my, const double mz,
    const size_t bin, const size_t volume);