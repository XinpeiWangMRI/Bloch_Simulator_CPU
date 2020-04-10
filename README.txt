Bloch simulator as a Matlab C++ MEX plugin, runs on CUDA.

Users,

The following features do not yet exist:
1) GUI
2) B1 map functionality
3) Relaxation
4) CPU parallelization

There is some internal error checking on input types from Matlab, vector lengths,
and to test if events are valid and magnetization vectors aren't diverging.

There is an abundant use of passing-by-const-value to ensure, even though the function
taking the parameter gets a copy, that even the copy is not modified. Because some of these values
are used in calculations of indices into arrays inside the function, accidentally modifying them
could attempt an access to undefined memory (rather, memory which may belong to a different process).

Also abundant pass-by-const-pointer. The matlab interface takes in pointers to the Matlab array arguments,
so sticking with pointers avoids casting everything in the mex gateway function. 

Help on any of the above would be greatly appreciated.

See INSTALL.txt for help with setup and a simple example.
