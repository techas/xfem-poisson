function h = calcMinEdgeLength( X, T )
% h = calcMinEdgeLength( X, T )
%
% Calculate the minimum edge length in a mesh.
%
% INPUT
%   X    nodal coords
%   T    connectivity matrix
%
% OUTPUT
%   h    min edge length
%
dmin = distanceSquare( X, T(1,1), T(1,2) );
for I = 1:size( T, 1 )
   d1 = distanceSquare( X, T(I,1), T(I,2) );
   d2 = distanceSquare( X, T(I,2), T(I,3) );
   d3 = distanceSquare( X, T(I,3), T(I,1) );
   d = min( [d1 d2 d3] );
   if d < dmin
      dmin = d;
   end
end
h = sqrt( dmin );


function d = distanceSquare( X, i, j )
%
%   X    nodal coords
%   i,j  nodes index
%
d = (X(i,1) - X(j,1))^2 + (X(i,2) - X(j,2))^2;


