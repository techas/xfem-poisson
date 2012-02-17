function bnd = IsOnBoundary( X, P )
% ISONBOUNDARY finds on which boundary the point P is. This functions works
% only for rectangular domains whose sides are aligned with the coordinates
% (x,y).
%
%  syntax: bnd = IsOnBoundary( X, P )
%
%  X: nodal coordinates [Nn*d matrix]
%  P: nodes to be tested [Np*d matrix]
%
%  bnd: value indicating the boundary on which the points are [Np*1 vector]
%       0 = not on boundary
%       1 = left
%       2 = up
%       3 = right
%       4 = down

% R. Cottereau 04/2008

% constants
tol = 1e-8;
nP = size( P, 1 );

% initializations
bnd = zeros( nP, 1 );

% define boundaries of the domain
x_lo = min( X(:,1) );
x_up = max( X(:,1) );
y_lo = min( X(:,2) );
y_up = max( X(:,2) );

% define boundary
bnd( find( abs( P(:,1) - x_lo ) < tol ) ) = 1;
bnd( find( abs( P(:,2) - y_up ) < tol ) ) = 2;
bnd( find( abs( P(:,1) - x_up ) < tol ) ) = 3;
bnd( find( abs( P(:,2) - y_lo ) < tol ) ) = 4;

