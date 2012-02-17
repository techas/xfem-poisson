function N = interpFunctionAtX( X, T, x, y, E, opts )
% INTERPFUNCTIONATX to find the value of the interpolation functions at
% several points [x y] in elements E
%
% syntax: N = interpFunctionAtX( X, T, x, y, E, opts )
%
%  X,T: nodal coordinates and connectivity matrix of the mesh
%  x,y: points where the interpolation functions are to be evaluated [n*1
%       [n*1 vectors]
%  E:   elements in which each point is contained [n*1 vector]
%
%  N:   matrix of the values of the interpolation functions for each point
%       [n*nV matrix], with nV the number of vertices for the mesh (given
%       in opts.numberOfElementNodes)

% R. Cottereau 04/2008

% checking input
if ~( opts.numberOfElementNodes==3 & opts.isTriangle==1 )
    error('INTERPFUNCTIONATX: option not implemented yet');
end

% initialization
nT = size( x, 1 );
N = zeros( nT, opts.numberOfElementNodes );

% coordinates of the vertices
X0 = X( T(E,1), : );
X1 = X( T(E,2), : ) - X0;
X2 = X( T(E,3), : ) - X0;

% maximum size of the matrices to be built
nMEM = 40;
NbM = floor( nT ./ nMEM );
NbR = mod( nT, nMEM );
nR = nT - NbM*nMEM;

% loop on matrices of the maximum size
if NbM > 0

    n2 = 2*nMEM;
    for i1=1:NbM

        ind = (i1-1)*nMEM + [1:nMEM];

        A = sparse( n2, n2 );
        A( 1:2:end , 1:2:end ) = spdiags( X1( ind, 1 ), 0, nMEM, nMEM );
        A( 2:2:end , 1:2:end ) = spdiags( X1( ind, 2 ), 0, nMEM, nMEM );
        A( 1:2:end , 2:2:end ) = spdiags( X2( ind, 1 ), 0, nMEM, nMEM );
        A( 2:2:end , 2:2:end ) = spdiags( X2( ind, 2 ), 0, nMEM, nMEM );

        B = zeros( n2, 1 );
        B( 1:2:end ) = x( ind ) - X0( ind, 1 );
        B( 2:2:end ) = y( ind ) - X0( ind, 2 );

        xieta = reshape( A \ B, 2, nMEM )';
        N( ind, : ) = [ 1-sum(xieta,2) xieta(:,1) xieta(:,2) ];

    end
end

% invertion for the rest of the lines
if NbR > 0

    n2 = 2*nR;
    ind = NbM*nMEM + [1:nR];

    A = sparse( n2, n2 );
    A( 1:2:end , 1:2:end ) = spdiags( X1( ind, 1 ), 0, nR, nR );
    A( 2:2:end , 1:2:end ) = spdiags( X1( ind, 2 ), 0, nR, nR );
    A( 1:2:end , 2:2:end ) = spdiags( X2( ind, 1 ), 0, nR, nR );
    A( 2:2:end , 2:2:end ) = spdiags( X2( ind, 2 ), 0, nR, nR );

    B = zeros( n2, 1 );
    B( 1:2:end ) = x( ind ) - X0( ind, 1 );
    B( 2:2:end ) = y( ind ) - X0( ind, 2 );

    xieta = reshape( A \ B, 2, nR )';
    N( ind, : ) = [ 1-sum(xieta,2) xieta(:,1) xieta(:,2) ];

end

