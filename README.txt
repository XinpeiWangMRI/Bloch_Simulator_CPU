Bloch simulator as a Matlab C++ MEX plugin, runs on CUDA.

Users,

It would be nice to make this simulator as general as possible, and to have a GUI. I have not yet added B1 map functionality, B0, or relaxation, although all should be straightforward additions.

Another nice feature would be to include a switch so that it compiles to be compatible with OpenMP (or other CPU parallelizing library) when not being compiled by the nvidia compiler.

One more feature that would be nice: automatically selecting the least busy GPU card (hereafter referred to as device) to run the simulation on. For now, the user must specify which GPU device to use. The mex file checks how many devices there are and, if the user provided number is available, uses that one. Otherwise it uses device 0.

Error checks. I essentially have none. These could be incorporated to run under-the-hood in a .m script before calling the simulator once there is a GUI in place.

Help on any of the above would be greatly appreciated.

As a side note, the first run usually seems to take substantially longer than subsequent runs. Maybe in initializing the GPU? Haven't looked into it much to be honest.

See INSTALL.txt for help with setup and a simple example.
