#pragma once

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

class magnetization {
  public:
    __host__ CUDA_CALLABLE_MEMBER magnetization(double mx0 = 0, double my0 = 0, double mz0 = 1, double x = 0,
            double y=0, double z=0, int index=0, double offres = 0, int volume = 1,
            int avg = 1);
    CUDA_CALLABLE_MEMBER void rotate(double bx, double by, double bz, double tstep);
    CUDA_CALLABLE_MEMBER void display();
    CUDA_CALLABLE_MEMBER void acquire(double* mxout, double* myout, double* mzout, int ndims,
            int time);
    __host__ CUDA_CALLABLE_MEMBER void setBin(int index,int numCols=1, int numRows=1,
        int numPages=1);
    CUDA_CALLABLE_MEMBER int getBin();
    __host__ CUDA_CALLABLE_MEMBER void setVolume(int numRows, int numCols, int numPages);
    CUDA_CALLABLE_MEMBER double getX();
    CUDA_CALLABLE_MEMBER double getY();
    CUDA_CALLABLE_MEMBER double getZ();
    __host__ CUDA_CALLABLE_MEMBER void setpos(int index, double* xgrid, double* ygrid, double* zgrid);
    int avg;
    CUDA_CALLABLE_MEMBER void set2eq();
    CUDA_CALLABLE_MEMBER void refocusM();
    __host__ void setobj(double objmz);
    __host__ void setOffset(double offset);
  private:
    double mx,my,mz,xpos,ypos,zpos,offres;
    int bin;
    int volume;
};
/*
CUDA_CALLABLE_MEMBER magnetization::magnetization(double mx0, double my0, double mz0, double x,
            double y, double z, int index, double offres, int volume, int avg) 
            : mx(mx0),my(my0),mz(mz0),xpos(x),ypos(y),zpos(z),bin(index),
                    offres(offres), volume(volume), avg(avg) {}

CUDA_CALLABLE_MEMBER void magnetization::rotate(double bx,double by, double bz, double tstep);

CUDA_CALLABLE_MEMBER void magnetization::display();

CUDA_CALLABLE_MEMBER int magnetization::getBin();

CUDA_CALLABLE_MEMBER void magnetization::setBin(int index,int numCols, int numRows,
        int numPages);

CUDA_CALLABLE_MEMBER void magnetization::acquire(double* mxout, double* myout, double* mzout, 
        int ndims, int time);

CUDA_CALLABLE_MEMBER void magnetization::setVolume(int numRows, int numCols, int numPages);

CUDA_CALLABLE_MEMBER double magnetization::getX();

CUDA_CALLABLE_MEMBER double magnetization::getY();

CUDA_CALLABLE_MEMBER double magnetization::getZ();
CUDA_CALLABLE_MEMBER void magnetization::setpos(int index,double* xgrid, double* ygrid, double* zgrid);
    
*/
