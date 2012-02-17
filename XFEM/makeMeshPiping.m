function makeMeshPiping()
% makeMeshWithHole()
%
% creates a rectangular mesh with a central hole
% and stores it in .mat file.
%
% The file contains the following variables
%   X   node coords
%   T   connetivity
%
%addpath( 'mesh2d' )
% Limits
xlo = -1;
ylo = -1;
xup = 1;
yup = 1;   

%
node = [ xlo  ylo
         xup  ylo
         xup  yup
         xlo  yup];
cnect = [ 1 2
          2 3
          3 4
          4 1]; 

elementSize = 0.03;
hdata.hmax = elementSize;

[X,T] = mesh2d( node, cnect, hdata );
fileName = 'meshPiping22.mat';
save( fileName, 'X', 'T' );
