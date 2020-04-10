%[amp,phase] = HSpulse(R,n,numPts,endFreq);
%amp is scaled to have maximum of 1, phase is in radians
%endFreq = 0 if end sweep on resonance, else sweeps through resonance

function [amp,phiRF,dphiRF] = HSpulse(R,n,numPts,endFreq)

beta = asech(.01);
tau = linspace(-1,1,numPts);

%beta = 5.29;

amp = sech(beta * tau.^n);
dphiRF = flip(cumsum(amp.^2));%(pi * R/(numPts)) * tanh(beta * tau);

if endFreq == 0
    dphiRF = dphiRF - dphiRF(end);
    dphiRF = dphiRF/max(abs(dphiRF));
    dphiRF = (R) * dphiRF;
else
    dphiRF = dphiRF/max(abs(dphiRF));
    dphiRF = dphiRF - 0.5;
    dphiRF = dphiRF/max(abs(dphiRF));
    dphiRF = .5*(R) * dphiRF;
end


phiRF = 2*pi*cumsum(dphiRF)/numPts;
%last else does nothing, so left out. Then frequency sweep ends at 0

if endFreq < 0
    phiRF = -phiRF;
end

if max(abs(diff(phiRF))) >= pi
   fprintf('WARNING: Phase step greater than pi; too few points in pulse. \n'); 
    
end
end