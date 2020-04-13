%Create unrealistically large numbers for rfNull and gradNull to indicate
%these waveforms will not be used during a specified sequence event.
classdef nullTypes < double
   enumeration
      rfNull        (-10000)
      gradNull      (10000)
   end
end