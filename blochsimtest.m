% To compile the simulator:
% mex mex_blochsim.cpp mexsimulator.cpp magnetization.cpp event_manager.cpp;
%
% Event types:
% pulse,        (RF pulse only)
% gradient,     (gradients only)
% pulseAndgrad  (RF pulse and gradients, internal check on whether
%                excitation or refocusing pulse)
% delay,        (standard sequence delay)
% acquire,      (requires gradient)
% pulseGradAcq, (simultaneous transmit/receive)
% thermaleq,    (restore thermal equilibrium and set dephasing to 0)
% refocus       (perfect refocusing pulse)

n = 128;
FOV = 19.2;
x = linspace(-FOV/2,FOV/2,n);
y = x;
z = 0;
obj = ones(size(phantom(n))); %since grid is 2D, object is too. Need 3D object for 3D simulation
%there is a package on the matlab file exchange for
%3d shepp logan phantom

[xgrid,ygrid,zgrid] = ndgrid(x,y,z);

nEvent = 2;
eventlist = zeros(nEvent,3);
%matlab won't let us assign the enumeration directly in a vector.
%it only permits single entry assignment directly to enumeration.
%e.g. eventlist(1,1) = EventTypes.pulseAndgrad is perfectly allowed.

eventlist(1,:) = [double(EventTypes.pulseAndgrad), 1000*4e-6, 1000];
eventlist(2,:) = [double(EventTypes.gradient), 500*4e-6, 500];
eventlist(3,:) = [double(EventTypes.acquisition) 1*4e-6 1];
eventlist(4,:) = [double(EventTypes.pulseAndgrad) 1000*4e-6 1000];
eventlist(5,:) = [double(EventTypes.gradient), 500*4e-6, 500];
eventlist(6,:) = [double(EventTypes.acquisition) 1*4e-6 1];

%only two pulses.
pulseTypeList = [logical(pulseTypes.excitation),logical(pulseTypes.refocus)];

Gx = [.5*ones(1000,1);
    -.5*ones(500,1);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull)];
Gy = [double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    .1*ones(1000,1); 
    0*-.1*ones(500,1);
    double(nullTypes.gradNull)];
Gz = [double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull);
    double(nullTypes.gradNull)];

[amp,phase] = HSpulse(10,1,1000,1); %r = 10, n = 1, 1000 time points, nonzero 4th entry means sweep goes through 0 offset
rfamp = [62.5*amp';
    double(nullTypes.rfNull);
    double(nullTypes.rfNull);
    30*62.5*amp';
    double(nullTypes.rfNull);
    double(nullTypes.rfNull)];
rfphase = [phase';
    double(nullTypes.rfNull);
    double(nullTypes.rfNull);
    phase';
    double(nullTypes.rfNull);
    double(nullTypes.rfNull)];

sig = 2;
B0max = 0*15000;
B0 = B0max * exp(-(xgrid.^2 + ygrid.^2 + zgrid.^2)/(2*sig^2));

fprintf('Simulating on grid distorted opposite what B0 causes. \n');
%xgrid = xgrid - B0/(4258 * .5);
%gradient only takes differences across columns, need it to be dependent on
%which grid I am looking at.
voxelX = permute(gradient(permute(xgrid,[2 1])),[2 1]);
voxelY = gradient(ygrid);
voxelZ = 0 * voxelX;

voxelWidths = cat(3,cat(3,voxelX,voxelY),voxelZ);

seqParams = struct();
seqParams.xgrid = xgrid;
seqParams.ygrid = ygrid;
seqParams.zgrid = zgrid;
seqParams.Gx = Gx;
seqParams.Gy = Gy;
seqParams.Gz = Gz;
seqParams.rfamp = rfamp;
seqParams.rfphase = rfphase;
seqParams.pulseTypeList = pulseTypeList;
seqParams.events = eventlist;
seqParams.usrObj = obj;
seqParams.B0 = B0;
[gB0y, gB0x] = gradient(B0);
gB0z = 0 * gB0x;
seqParams.B0gradients = cat(3,cat(3,gB0x,gB0y),gB0z);

%avoid having to take all the spatial gradients in C++...
seqParams.VoxelWidths = voxelWidths;

tic;
[mx,my,mz] = mex_blochsim(seqParams);
totaltime = toc;
mxy = mx + 1i*my;

figure(42);
subplot(1,2,1);
imagesc((angle(mxy(:,:,2))));
axis square;
colorbar;

subplot(1,2,2);
imagesc(abs(mxy(:,:,2)));
colorbar;
axis square;
