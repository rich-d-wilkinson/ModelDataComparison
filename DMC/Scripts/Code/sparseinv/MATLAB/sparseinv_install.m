function sparseinv_install
%SPARSEINV_INSTALL compiles and installs the sparseinv function.
% Your current working directory must be the sparseinv directory for this
% function to work.
%
% Example:
%   sparseinv_install
%
% See also sparseinv, sparseinv_test

% Copyright 2011, Timothy A. Davis, University of Florida

help sparseinv_install
is64 = ~isempty (strfind (computer, '64')) ;
if (is64)
    mex -largeArrayDims sparseinv_mex.c sparseinv.c
else
    mex sparseinv_mex.c sparseinv.c
end
addpath (pwd)

fprintf ('Added the following directory to your path:\n') ;
disp (pwd) ;
fprintf ('Use pathtool and click "save" to save this for future sessions.\n') ;
