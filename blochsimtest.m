% mexcuda mex_blochsim.cu mexsimulator.cu magnetization.cu event_manager.cu -dynamic;
% the first run always takes about 10x longer than subsequent runs. Maybe
% in initializing the GPU?
% types of events: 0 - 7 corresponding to (in order)
% pulse,gradient,pulseAndgrad,delay,acquire,pulseGradAcq,thermaleq,refocus 
avg = 2;
n = 128;
x = linspace(-19.2/2,19.2/2,avg*n);
y = x;
z = 0;
obj = phantom(n*avg); %since grid is 2D, object is too. Need 3D object for 3D simulation
                      %there is a package on the matlab file exchange for
                      %3d shepp logan phantom

[xgrid,ygrid,zgrid] = ndgrid(x,y,z);
nEvent = 2;
eventlist = zeros(nEvent,3);
eventlist(1,:) = [2 1500*4e-6 1500];
eventlist(2,:) = [4 1*4e-6 1];
eventlist(3,:) = [2 1500*4e-6 1500];
eventlist(4,:) = [7 1*4e-6 1]; %refocusing event
eventlist(5,:) = [4 1*4e-6 1];

gradNull = 10000;
rfNull = -10000;

Gx = [.5*ones(1000,1);-.5*ones(500,1);gradNull;gradNull;gradNull];
Gy = [gradNull;gradNull;.1*ones(1000,1);-.1*ones(500,1);gradNull;gradNull];
Gz = [gradNull;gradNull;gradNull;gradNull;gradNull];
[amp,phase] = HSpulse(10,1,1000,1); %r = 10, n = 1, 1000 time points, nonzero 4th entry means sweep goes through 0 offset
rfamp = [62.5*amp';zeros(500,1);rfNull;62.5*amp';zeros(500,1);rfNull;rfNull];
rfphase = [phase';zeros(500,1);rfNull;phase';zeros(500,1);rfNull;rfNull];
cudaDevNum = 0; %should be less than 4 on gpu4gar. There is an internal check if used on other systems

seqParams = struct();
seqParams.avg = avg;
seqParams.xgrid = xgrid;
seqParams.ygrid = ygrid;
seqParams.zgrid = zgrid;
seqParams.Gx = Gx;
seqParams.Gy = Gy;
seqParams.Gz = Gz;
seqParams.rfamp = rfamp;
seqParams.rfphase = rfphase;
seqParams.events = eventlist;
seqParams.usrObj = obj;

tic;
[mx,my,mz] = mex_blochsim(seqParams,cudaDevNum);
totaltime = toc;
mxy = mx + 1i*my;

figure();
imagesc((abs(mxy(:,:,2))));

