#include "magnetization.hpp"
#include "math.h"
#include "mexsimulator.hpp"
#include "event_manager.hpp"
#include "mex.h"
#include "constants.hpp"
#include <vector>

void mexsimulator(std::vector<magnetization>& magn, double* mxout, double* myout,
        double* mzout, const size_t nelements,  const size_t numCols, const size_t numRows,
        const size_t numPages, const size_t ndims, const double* gradx, const double* grady, const double* gradz,
        const double* rfpulse, const double* rfphase, const double* events, const size_t nEvents) {

    //Starting times in each gradient + rfwaveform
    size_t xstart = 0, ystart = 0, zstart = 0, rfstart = 0;

    //index into acquisition periods
    size_t acqstart = 0;

    // index to access values in *events 
    size_t eventdurationIndex;
    size_t eventstepsIndex;
    
    //Might need to double cast to get around inability to cast from double to enumeration.
    //The events array has some elements which are necessarily doubles and some which are integers.
    //Inputting directly as a struct may be beneficial. Then can natively store the event type as the enumeration type.
    int tempEvent = 0;

    for (int eventIndex = 0; eventIndex < nEvents; eventIndex++) {
        eventdurationIndex = eventIndex + nEvents;
        eventstepsIndex = eventIndex + 2 * nEvents;

        tempEvent = static_cast<int>(events[eventIndex]);
        mrEvent myEvent((typeEvent)tempEvent, events[eventdurationIndex],
            (size_t)events[eventstepsIndex]);

        //Only perform simulation if the event is valid, otherwise valid_Event throws an error directly to Matlab.
        if (valid_Event(myEvent)) {

            //Loop over all space.
            for (int magnIndex = 0; magnIndex < magn.size(); ++magnIndex) {
                switch (myEvent.getEvent()) {

                case pulse:
                    myEvent.pulseEvent(&magn[magnIndex], rfstart, rfpulse, rfphase);
                    break;

                case gradient:
                    myEvent.gradEvent(&magn[magnIndex], xstart, ystart, zstart, gradx,
                        grady, gradz);
                    break;

                case pulseAndgrad:
                    myEvent.pulsegradEvent(&magn[magnIndex], xstart, ystart, zstart, gradx,
                        grady, gradz, rfpulse, rfphase, rfstart);
                    break;

                case delay:
                    myEvent.delayEvent(&magn[magnIndex]);
                    break;

                case acquisition:
                    myEvent.acquireEvent(&magn[magnIndex], xstart, ystart, zstart,
                        gradx, grady, gradz, mxout, myout, mzout,
                        ndims, acqstart);
                    break;

                case pulseGradAcq:
                    myEvent.pulseGradAcqEvent(&magn[magnIndex], xstart, ystart, zstart,
                        gradx, grady, gradz, mxout, myout, mzout,
                        ndims, acqstart, rfpulse, rfphase, rfstart);
                    break;

                case thermaleq:
                    myEvent.thermaleqEvent(&magn[magnIndex]);
                    break;

                case refocus:
                    myEvent.refocusEvent(&magn[magnIndex]);
                    break;

                }
            }
            
            //Update starting point in the various waveform arrays 
            if (rfpulse[rfstart] == rfNull) {
                indexUpdate(&rfstart, 1);
            }
            else {
                indexUpdate(&rfstart, myEvent.getnSteps());
            };

            if (gradx[xstart] == gradNull) {
                indexUpdate(&xstart, 1);
            }
            else {
                indexUpdate(&xstart, myEvent.getnSteps());
            };

            if (grady[ystart] == gradNull) {
                indexUpdate(&ystart, 1);
            }
            else {
                indexUpdate(&ystart, myEvent.getnSteps());
            };

            if (gradz[zstart] == gradNull) {
                indexUpdate(&zstart, 1);
            }
            else {
                indexUpdate(&zstart, myEvent.getnSteps());
            };

            if (myEvent.getEvent() == acquisition) {
                indexUpdate(&acqstart, myEvent.getnSteps());
            };
        };
    };
};