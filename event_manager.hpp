//Header file for handling different events in the pulse sequence 

#pragma once

#include "magnetization.hpp"
#include "constants.hpp"
#include "mex.h"

// ======================================================================
// Gradients, pulses, delays, acquisition are input as vectors of fixed 
// length, which need not all be the same length. There is a separate event 
// array, structured as (going across columns):
//      event type, duration, number of steps 
// The Bloch simulator steps through the events, and updates the 
// magnetization vectors accordingly. The counters (ctr) below update the
 // correct indices of the input vectors for the current event.
 // ======================================================================

class mrEvent {
    public:
        //Constructor, needs arguments!
        mrEvent(typeEvent eventtype = pulse, double duration = temporal_Res_Std, size_t nSteps = 1);

        // return type of event 
        typeEvent getEvent() const { return eventtype; };

        // return time step used for this event 
        double gettstep() const { return tstep; };

        //The actual types of events.
        void delayEvent(magnetization* magn);
        void pulseEvent(magnetization* magn, const size_t rfstart, const double* rfpulse,
            const double* rfphase);
        void gradEvent(magnetization* magn, const size_t xstart, const size_t ystart,
            const size_t zstart, const double* gradx, const double* grady,
            const double* gradz);
        void pulsegradEvent(magnetization* magn, const size_t xstart, const size_t ystart,
            const size_t zstart, const double* gradx, const double* grady,
            const double* gradz, const double* rfpulse, const double* rfphase, const size_t rfstart);
        void pulseGradAcqEvent(magnetization* magn, const size_t xstart, const size_t ystart,
            const size_t zstart, const double* gradx, const double* grady,
            const double* gradz, double* mxout, double* myout, double* mzout,
            const size_t ndims, const size_t timestart, const double* rfpulse,
            const double* rfphase, const size_t rfstart);
        void acquireEvent(magnetization* magn, const size_t xstart, const size_t ystart,
            const size_t zstart, const double* gradx, const double* grady,
            const double* gradz, double* mxout, double* myout, double* mzout,
            const size_t ndims, const size_t timestart);
        void thermaleqEvent(magnetization* magn);
        void refocusEvent(magnetization* magn);

        //return number of steps in event
        size_t getnSteps() const { return nSteps; }

    private:
        typeEvent eventtype;
        double duration; //duration of event in seconds
        size_t nSteps;      //number of time points used to define duration 
        double tstep = duration / nSteps;
};

//Helper Functions

// Update the necessary array indices in sharedTimeCounter for the corresponding type of event 
void indexUpdate(size_t* start, const size_t nSteps);

//Will double check for positive duration and number of steps >= 1
bool valid_Event(const mrEvent& currentEvent);