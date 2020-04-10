#pragma once

//some global constants. Maybe make a struct with these values?
const double gradNull = 10000; /* No gradient is realistically this big */
const double rfNull = -10000;  /* rf amplitude should always be positive, doing 10000 to help remember it */
const double gauss2Hz = 4258.0;
const double temporal_Res_Std = 4e-6; //temporal resolution of Varian MRI systems' gradient amplifiers (time limiting hardware)

//Not a big fan of using this as a global enumeration, thinking of other ways.
//Various types of events in MRI sequences.
const enum typeEvent {
    pulse, gradient, pulseAndgrad, delay, acquisition, pulseGradAcq,
    thermaleq, refocus
};