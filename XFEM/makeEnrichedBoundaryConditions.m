function [Abc,bbc] = makeEnrichedBoundaryConditions( X, T, levelSet, ...
                              ksize, enrichedNodes, caseBoundaryCondition)
%
% Build the boundary conditions for the enriched nodes.
%
% INPUT
%   X,T       the mesh
%   T         conectivity matrix
%   levelSet  nodal level set 
%   ksize     size of the K matrix
%   types     element classification vector (returned by classifyElements)
%   caseBoundaryCondition type of boundary conditions 
%                         (see MAKEBOUNDARYCONDITIONS for details)
%
% OUTPUT
%   Abc,bbc   lagrange multipliers
%
x_lo = min( X(:,1) );
x_up = max( X(:,1) );
y_lo = min( X(:,2) );
y_up = max( X(:,2) );
R2 = max( sum(X.^2,2) );
%
bottomNodes = find( X(enrichedNodes,2) == y_lo );
topNodes = find( X(enrichedNodes,2) == y_up );
leftNodes = find( X(enrichedNodes,1) == x_lo );
rigthNodes = find( X(enrichedNodes,1) == x_up );
outNodes = find( abs(sum(X.^2,2)-R2) < 1e-8 ); 

switch caseBoundaryCondition
    case 1
        C = [topNodes zeros( length( topNodes ), 1 ); ...
             leftNodes zeros( length( leftNodes ), 1 ); ...
             rigthNodes zeros( length( rigthNodes ), 1 ); ...
             bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 2
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 3
        C = [rigthNodes zeros( length( rigthNodes ), 1 ); ...
             leftNodes zeros( length( leftNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 4
        C = [topNodes zeros( length( topNodes ), 1 ); ...
             bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes )  );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 5
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 6
        C = [leftNodes zeros( length( leftNodes ), 1 ); ...
             rigthNodes zeros( length( rigthNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 7
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 8
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 9
        C = [rigthNodes zeros( length( rigthNodes ), 1 ); ...
             leftNodes zeros( length( leftNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);

    case 10
%         C = [];
%         %
%         nBC = size( C, 1 );
%         Abc = zeros( nBC, length(enrichedNodes) );
%         Abc(:,C(:,1)) = eye( nBC );
%         bbc = C(:,2);
        Abc = [];
        bbc = [];

    case 11
        C = [rigthNodes zeros( length( rigthNodes ), 1 ); ...
             leftNodes zeros( length( leftNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, length( enrichedNodes ) );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
end;



