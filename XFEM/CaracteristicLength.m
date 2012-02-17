function Lcar = CaracteristicLength( X, T )
% CARACTERISTICLENGTH to compute a characteristic length for each of the
% nodes of the mesh
%
% syntax: Lcar = CaracteristicLength( X, T )
%
%  X: nodal coordinates [Nn*d matrix]
%  T: connectivity matrix [Ne*e matrix]
%
%  Lcar: characteristic lengths [Ne*1 vector]

% R. Cottereau 04/2008

% constants
[ Ne d ] = size( X );
e = size( T, 2 );

% initializations
Lcar = zeros( Ne, 1 );

% loop on the nodes of the mesh
for i1 = 1: Ne
    [ touchingElts, j1 ] = find( T == i1 );
    Nte = length( touchingElts );
    touchingElts = T( touchingElts, : )';
    x = reshape( X( touchingElts, 1), e, Nte );
    y = reshape( X( touchingElts, 2), e, Nte );
%    s = polyarea( x, y );
    % the characteristic length is the minimum square root of the area of
    % the elements in contact with the corresponding node
    Lcar( i1, 1 ) = min( sqrt( polyarea( x, y ) ) ) / d;
end
