function [ t, x ] = InContactWithNodes( T, x0 )
% INCONTACTWITHNODES to search the elements in contact with a set of nodes
% simultaneously, and the corresponding list of nodes
%
% syntax: [ t, x ] = InContactWithNodes( T, x0 )
%
%  T:  connectivity matrix
%  x0: set of nodes to be searched for

% R. Cottereau 05/2008

% constants;
x0 = x0(:);
n0 = length(x0);
Ne = size( T, 1 );

% initialization
t = 1:Ne;

% loop on the nodes in x0
for i1 = 1:n0
    [ ind1, j1 ] = find( T( t, : ) == x0( i1 ) );
    t = t( ind1 );
end

% list of corresponding nodes
x = T( t, : );
x = unique( x(:) );