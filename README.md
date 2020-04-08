# PhaseUnwrapping
2D&3D phase unwrapping plugin for Matlab

To compile the 2D unwrapper only in Matlab, type
mex phaseUnwrap2D.cpp

To compile the 2D/3D unwrapper, use
mex phaseUnwrap.cpp

then for a given wrapped (double) phase image, type
unwrapped = phaseUnwrap2D(wrapped); 
OR
unwrapped = phaseUnwrap(wrapped);

to get the unwrapped image.
