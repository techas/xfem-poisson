function makeUniformMesh( fileName, limits, elementSize )
% makeUniformMesh( fileName, limits, elementSize )
%
% creates a uniform rectangular mesh and stores it in .mat file.
%
% INPUT
%   fileName
%   limits = [x_lo x_up y_lo y_up]
%   elementSize
%
% The file contains the following variables
%   X   node coords
%   T   connetivity
%
addpath( 'mesh2d' )
if nargin == 1
   xlo = 0;
   ylo = 0;
   xup = 1;
   yup = 1;   
else
   xlo = limits(1);
   ylo = limits(3);
   xup = limits(2);
   yup = limits(4);
end
if nargin < 3
   elementSize = 0.1;
end
% Imposed nodes
node = [xlo ylo;
        xup ylo;
        xup yup;
        xlo yup];
% Imposed edges
edge = [];
% Max element size
hdata.hmax = elementSize;
% Size function
% hdata.fun = @hfun1;
[X,T] = mesh2d( node, edge, hdata );
save( fileName, 'X', 'T' );
fprintf( 'nn=%g\nne=%g\n', size( X, 1 ), size( T, 1 ) )


function d = hfun1( x, y )
%
% determines the expected size of the element at the position (x,y)
%
dlo = 200;
dup = 100;
d = ones( size( y ) ) * dlo;
d(y > 400) = dup;

