function makeMeshWithHole()
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
xlo = 0;
ylo = -1.2;
xup = 1.2;
yup = 0;   
% hole 1
radius = 0.05;
x0 = 0.9; y0 = -0.7;
step = 8;
theta = 0:pi/step:(2*pi-pi/step);
x1 = x0 + cos( theta )*radius;
y1 = y0 + sin( theta )*radius;
%
node = [ x1'  y1'
        xlo  ylo
        xup  ylo
        xup  yup
        xlo  yup];
n = size( node, 1 ) - 4;
cnect = [(1:n-1)' (2:n)'
          n        1
          n+1      n+2
          n+2      n+3
          n+3      n+4
          n+4      n+1]; 
% 
elementSize = 0.02;
hdata.hmax = elementSize;
% hdata.fun = @hfun2;
[X,T] = mesh2d( node, cnect, hdata );
fileName = 'meshPaper3.mat';
save( fileName, 'X', 'T' );
%fprintf( 'nn=%g\nne=%g\n', size( X, 1 ), size( T, 1 ) )


function d = hfun1( x, y )
%
% determines the expected size of the element at the position (x,y)
%
dlo = 200;
dup = 100;
d = ones( size( y ) ) * dlo;
d(y > 400) = dup;

