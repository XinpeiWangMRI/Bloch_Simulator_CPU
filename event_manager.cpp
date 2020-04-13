#include "event_manager.hpp"
#include "math.h"
#include "constants.hpp"

// ===========================================================================
// Gradients, pulses, delays, acquisition are input as vectors of fixed 
// length, which need not all be the same length. There is a separate event 
// array, structured as (going across columns):
//      event type, duration, number of steps 
// The Bloch simulator steps through the events, and updates the 
// magnetization vectors accordingly. The counters (ctr) below update the
// correct indices of the input vectors for the current event.
// ===========================================================================


//Constructor, needs to be fed values.
mrEvent::mrEvent(typeEvent eventtype1, double duration1, size_t nSteps1) :
    eventtype(eventtype1), duration(duration1), nSteps(nSteps1) {};
    
// Don't need to define this until off-resonance map is included 
void mrEvent::delayEvent(magnetization* magn) {
    for (size_t index = 0; index < nSteps; index++) {
        magn->rotate(0, 0, magn->getOffset(), tstep); //0 could be replaced by offres
    }
};

void mrEvent::pulseEvent(magnetization* magn, const size_t rfstart, const double* rfpulse, const double* rfphase){
    double rfx, rfy;
    double oldMz = magn->getMz();
    
    for (size_t index = rfstart; index < rfstart+nSteps; index++){
        if (rfpulse[rfstart] != rfNull){
            rfx = cos(rfphase[index+rfstart]) * rfpulse[index+rfstart];
            rfy = sin(rfphase[index+rfstart]) * rfpulse[index+rfstart];
        }
        else{
            rfx = 0;
            rfy = 0;
        }
        
        magn->rotate(rfx,rfy, magn->getOffset(),tstep); //0 could be replaced by offres
    }

    if (is_Ref_Pulse(magn->getMz(),oldMz)) magn->reverseDephase();
    
};

void mrEvent::gradEvent(magnetization* magn, const size_t xstart, const size_t ystart,
        const size_t zstart, const double* gradx, const double* grady,
        const double* gradz) {
    double bx, by, bz;
    double appliedGrad[3];

    for (int index = 0; index < nSteps; index++){
        // No pretty way to implement this. Could do a bunch of cases, but 
        // there would be 2^3 = 8 possibilities, since any number of them 
        // could be null. 
         
        if (gradx[xstart] != gradNull) {
            bx = gauss2Hz * gradx[index + xstart] * magn->getX();
            appliedGrad[0] = gradx[index + xstart];
        }
        else {
            bx = 0;
            appliedGrad[0] = 0;
        }
        if (grady[ystart] != gradNull) {
            by = gauss2Hz * grady[index + ystart] * magn->getY();
            appliedGrad[1] = grady[index + ystart];

        }
        else {
            by = 0;
            appliedGrad[1] = 0;

        }
        if (gradz[zstart] != gradNull) {
            bz = gauss2Hz * gradz[index + zstart] * magn->getZ();
            appliedGrad[2] = gradz[index + zstart];

        }
        else {
            bz = 0;
        }

        magn->rotate(0,0,bx+by+bz+ magn->getOffset(),tstep, appliedGrad);
    }
    
};

void mrEvent::pulsegradEvent(magnetization* magn, const size_t xstart, const size_t ystart,
        const size_t zstart, const double* gradx, const double* grady,
        const double* gradz, const double* rfpulse, const double* rfphase, const size_t rfstart) {

    double bx, by, bz, rfx, rfy;
    double oldMz = magn->getMz();

    double appliedGrad[3];

    for (int index = 0; index < nSteps; index++){
        // No pretty way to implement this. Could do a bunch of cases, but 
        // there would be 2^3 = 8 possibilities, since any number of them 
        // could be null. 
         
        if (rfpulse[rfstart] != rfNull){
            rfx = cos(rfphase[index+rfstart]) * rfpulse[index+rfstart];
            rfy = sin(rfphase[index+rfstart]) * rfpulse[index+rfstart];
        }
        else{
            rfx = 0;
            rfy = 0;
        }
        
        if (gradx[xstart] != gradNull) {
            bx = gauss2Hz * gradx[index + xstart] * magn->getX();
            appliedGrad[0] = gradx[index + xstart];
        }
        else {
            bx = 0;
            appliedGrad[0] = 0;
        }
        if (grady[ystart] != gradNull) {
            by = gauss2Hz * grady[index + ystart] * magn->getY();
            appliedGrad[1] = grady[index + ystart];

        }
        else {
            by = 0;
            appliedGrad[1] = 0;

        }
        if (gradz[zstart] != gradNull) {
            bz = gauss2Hz * gradz[index + zstart] * magn->getZ();
            appliedGrad[2] = gradz[index + zstart];

        }
        else {
            bz = 0;
        }

        magn->rotate(rfx,rfy,bx+by+bz+ magn->getOffset(),tstep, appliedGrad);
    };

    if (is_Ref_Pulse(magn->getMz(), oldMz)) magn->reverseDephase();

};

void mrEvent::pulseGradAcqEvent(magnetization* magn, const size_t xstart, const size_t ystart,
        const size_t zstart, const double* gradx, const double* grady,
        const double* gradz, double* mxout, double* myout, double* mzout,
        const size_t ndims, const size_t timestart, const double* rfpulse,
        const double* rfphase, const size_t rfstart) {
    
    double bx, by, bz, rfx, rfy;
    double oldMz = magn->getMz();

    //Holds applied gradients for aquisition statement.
    double appliedGrad[3];
    double kcoord[3];

    for (int index = 0; index < nSteps; index++){
        // No pretty way to implement this. Could do a bunch of cases, but 
        // there would be 2^3 = 8 possibilities, since any number of them 
        // could be null. 
         
        if (rfpulse[rfstart] != rfNull){
            rfx = cos(rfphase[index+rfstart]) * rfpulse[index+rfstart];
            rfy = sin(rfphase[index+rfstart]) * rfpulse[index+rfstart];
        }
        else{
            rfx = 0;
            rfy = 0;
        }

        if (gradx[xstart] != gradNull){
            bx = gauss2Hz * gradx[index+xstart] * magn->getX();
            appliedGrad[0] = gradx[index + xstart];
        }
        else{
            bx = 0;
            appliedGrad[0] = 0;
        }
        if (grady[ystart] != gradNull){
            by = gauss2Hz * grady[index+ystart] * magn->getY();
            appliedGrad[1] = grady[index + ystart];

        }
        else{
            by = 0;
            appliedGrad[1] = 0;

        }
        if (gradz[zstart] != gradNull){
            bz = gauss2Hz * gradz[index+zstart] * magn->getZ();
            appliedGrad[2] = gradz[index + zstart];

        }
        else{
            bz = 0;
        }

        magn->rotate(rfx,rfy,bx+by+bz + magn->getOffset(),tstep,appliedGrad);
        kcoord[0] = magn->getKcoord(0);
        kcoord[1] = magn->getKcoord(1);
        kcoord[2] = magn->getKcoord(2);

        acquire(mxout, myout, mzout, ndims, (index + timestart), magn->getMx(), magn->getMy(), magn->getMz(),
            magn->getBin(), magn->getVolume(), kcoord);

    };
    if (is_Ref_Pulse(magn->getMz(), oldMz)) magn->reverseDephase();

};

void mrEvent::acquireEvent(magnetization* magn, const size_t xstart, const size_t ystart,
        const size_t zstart, const double* gradx, const double* grady,
        const double* gradz, double* mxout, double* myout, double* mzout,
        const size_t ndims, const size_t timestart) {
    double bx, by, bz;
    double appliedGrad[3];
    double kcoord[3];

    for (int index = 0; index < nSteps; index++){
        // No pretty way to implement this. Could do a bunch of cases, but 
        // there would be 2^3 = 8 possibilities, since any number of them 
        // could be null. 
         
        if (gradx[xstart] != gradNull) {
            bx = gauss2Hz * gradx[index + xstart] * magn->getX();
            appliedGrad[0] = gradx[index + xstart];
        }
        else {
            bx = 0;
            appliedGrad[0] = 0;
        }
        if (grady[ystart] != gradNull) {
            by = gauss2Hz * grady[index + ystart] * magn->getY();
            appliedGrad[1] = grady[index + ystart];

        }
        else {
            by = 0;
            appliedGrad[1] = 0;

        }
        if (gradz[zstart] != gradNull) {
            bz = gauss2Hz * gradz[index + zstart] * magn->getZ();
            appliedGrad[2] = gradz[index + zstart];

        }
        else {
            bz = 0;
        }
        
        magn->rotate(0,0,bx+by+bz + magn->getOffset(),tstep, appliedGrad);
        kcoord[0] = magn->getKcoord(0);
        kcoord[1] = magn->getKcoord(1);
        kcoord[2] = magn->getKcoord(2);

        acquire(mxout, myout, mzout, ndims, (index + timestart), magn->getMx(), magn->getMy(), magn->getMz(),
            magn->getBin(), magn->getVolume(),kcoord);
    };
    
};

void mrEvent::thermaleqEvent(magnetization* magn){
    //beware, this also resets dephaseTime to 0!
    magn->set2eq();
};
void mrEvent::refocusEvent(magnetization* magn){
    magn->refocusM();
    //refocus event is by definition a perfect refocusing pulse, so reverse the dephase time counter.
    magn->reverseDephase();
};

//Helper Functions
void indexUpdate(size_t* start, const size_t nSteps) {
    *start = *start + nSteps;
}

//Will double check for positive duration and number of steps >= 1
bool valid_Event(const mrEvent& currentEvent) {
    if (currentEvent.getnSteps() <= 0 || currentEvent.gettstep() <= 0) {
        mexErrMsgIdAndTxt("MATLAB:event_manager:eventType",
            "Events must have positive duration and 1 or more steps");
        return false; //redundant
    }
    return true;
}