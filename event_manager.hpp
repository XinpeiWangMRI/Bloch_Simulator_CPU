<<<<<<< HEAD
/*Header file for handling different events in the pulse sequence */

#pragma once

#include "magnetization.hpp"
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 
/*
 * Gradients, pulses, delays, acquisition are input as vectors of fixed 
 * length, which need not all be the same length. There is a separate event 
 * array, structured as (going across columns):
 *      event type, duration, number of steps 
 * The Bloch simulator steps through the events, and updates the 
 * magnetization vectors accordingly. The counters (ctr) below update the
 * correct indices of the input vectors for the current event.
 */

//wrap these into a structure which can be passed between functions.
static int pulsectr, gradxctr, gradyctr, gradzctr, delayctr, acquirectr;
static int rfstart, rfend, gradxstart, gradxend, gradystart, gradyend, 
        gradzstart, gradzend; 

enum typeEvent {pulse,gradient,pulseAndgrad,delay,acquire,pulseGradAcq,
                thermaleq,refocus};  

class mrEvent {
  public:
      /* no defaults, want user to supply everything intentionally */
    CUDA_CALLABLE_MEMBER mrEvent(typeEvent eventtype = pulse, double duration = 0, int nSteps = 1); 
    
    CUDA_CALLABLE_MEMBER typeEvent getEvent(); /* return type of event */
    CUDA_CALLABLE_MEMBER double gettstep(); /* return time step used for this event */
    
    /* Update the necessary array indices for the corresponding type of event */
    CUDA_CALLABLE_MEMBER void indexUpdate(int* start, int nSteps);
    CUDA_CALLABLE_MEMBER void delayEvent(magnetization* magn);
    CUDA_CALLABLE_MEMBER void pulseEvent(magnetization* magn, int rfstart, double* rfpulse,
            double* rfphase);
    CUDA_CALLABLE_MEMBER void gradEvent(magnetization* magn, int xstart, int ystart, 
            int zstart, double* gradx, double* grady, 
            double* gradz);
    CUDA_CALLABLE_MEMBER void pulsegradEvent(magnetization* magn, int xstart, int ystart, 
            int zstart, double* gradx, double* grady, 
            double* gradz, double* rfpulse, double* rfphase, int rfstart);
    CUDA_CALLABLE_MEMBER void pulseGradAcqEvent(magnetization* magn, int xstart, int ystart,
        int zstart, double* gradx, double* grady,
        double* gradz, double* mxout, double* myout, double* mzout,
        int ndims, int timestart,double* rfpulse, 
        double* rfphase,int rfstart);
    CUDA_CALLABLE_MEMBER void acquireEvent(magnetization* magn, int xstart, int ystart, 
            int zstart, double* gradx, double* grady, 
            double* gradz, double* mxout, double* myout, double* mzout,
            int ndims, int timestart);
    CUDA_CALLABLE_MEMBER int getnSteps();
    CUDA_CALLABLE_MEMBER void thermaleqEvent(magnetization* magn);
    CUDA_CALLABLE_MEMBER void refocusEvent(magnetization* magn);
  private:
    typeEvent eventtype;
    double duration; /*duration of event in seconds*/
    int nSteps;      /*number of time points used to define duration */
    double tstep = duration/nSteps;
=======
/*Header file for handling different events in the pulse sequence */

#pragma once

#include "magnetization.hpp"
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 
/*
 * Gradients, pulses, delays, acquisition are input as vectors of fixed 
 * length, which need not all be the same length. There is a separate event 
 * array, structured as (going across columns):
 *      event type, duration, number of steps 
 * The Bloch simulator steps through the events, and updates the 
 * magnetization vectors accordingly. The counters (ctr) below update the
 * correct indices of the input vectors for the current event.
 */

static int pulsectr, gradxctr, gradyctr, gradzctr, delayctr, acquirectr;
static int rfstart, rfend, gradxstart, gradxend, gradystart, gradyend, 
        gradzstart, gradzend; 

enum typeEvent {pulse,gradient,pulseAndgrad,delay,acquire,pulseGradAcq,
                thermaleq,refocus};  

class mrEvent {
  public:
      /* no defaults, want user to supply everything intentionally */
    CUDA_CALLABLE_MEMBER mrEvent(typeEvent eventtype = pulse, double duration = 0, int nSteps = 1); 
    
    CUDA_CALLABLE_MEMBER typeEvent getEvent(); /* return type of event */
    CUDA_CALLABLE_MEMBER double gettstep(); /* return time step used for this event */
    
    /* Update the necessary array indices for the corresponding type of event */
    CUDA_CALLABLE_MEMBER void indexUpdate(int* start, int nSteps);
    CUDA_CALLABLE_MEMBER void delayEvent(magnetization* magn);
    CUDA_CALLABLE_MEMBER void pulseEvent(magnetization* magn, int rfstart, double* rfpulse,
            double* rfphase);
    CUDA_CALLABLE_MEMBER void gradEvent(magnetization* magn, int xstart, int ystart, 
            int zstart, double* gradx, double* grady, 
            double* gradz);
    CUDA_CALLABLE_MEMBER void pulsegradEvent(magnetization* magn, int xstart, int ystart, 
            int zstart, double* gradx, double* grady, 
            double* gradz, double* rfpulse, double* rfphase, int rfstart);
    CUDA_CALLABLE_MEMBER void pulseGradAcqEvent(magnetization* magn, int xstart, int ystart,
        int zstart, double* gradx, double* grady,
        double* gradz, double* mxout, double* myout, double* mzout,
        int ndims, int timestart,double* rfpulse, 
        double* rfphase,int rfstart);
    CUDA_CALLABLE_MEMBER void acquireEvent(magnetization* magn, int xstart, int ystart, 
            int zstart, double* gradx, double* grady, 
            double* gradz, double* mxout, double* myout, double* mzout,
            int ndims, int timestart);
    CUDA_CALLABLE_MEMBER int getnSteps();
    CUDA_CALLABLE_MEMBER void thermaleqEvent(magnetization* magn);
    CUDA_CALLABLE_MEMBER void refocusEvent(magnetization* magn);
  private:
    typeEvent eventtype;
    double duration; /*duration of event in seconds*/
    int nSteps;      /*number of time points used to define duration */
    double tstep = duration/nSteps;
>>>>>>> fece1db5dc82c1259b70f57c771f47a4d54b8c8a
};