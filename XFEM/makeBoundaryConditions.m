function [Abc,bbc] = makeBoundaryConditions( X, numberOfUnknows, ...
                                                    caseBoundaryCondition )
% [Abc,bbc] = makeBoundaryConditions( X, numberOfUnknows,
%                                                    caseBoundaryCondition)
%
% creates lagrange multipliers
%
% INPUT
%   X                nodal coords
%   numberOfUnknows
%   caseBoundaryCondition type of boundary conditions
%                         1: test cases with h=0 on all sides (to use with
%                            a rectangular mesh)
%                         2: test cases with h=0 on bottom and q=0 on
%                            other sides (to use with a rectangular mesh)
%                         3: test cases with h=0 on the left, h==5e-2 on the
%                            right and q=0 on other sides (to use with a 
%                            rectangular mesh)
%                         4: test cases with h=0 below, h==1 on top and q=0
%                            on other sides (rectangular mesh)
%
% OUTPUT
%   Abc,bbc  matrix and vector 
%
x_lo = min( X(:,1) );
x_up = max( X(:,1) );
y_lo = min( X(:,2) );
y_up = max( X(:,2) );
R2 = max( sum(X.^2,2) );
%
bottomNodes = find( X(:,2) == y_lo);
topNodes = find( X(:,2) == y_up);
leftNodes = find( X(:,1) == x_lo );
rigthNodes = find( X(:,1) == x_up );
outNodes = find( abs(sum(X.^2,2)-R2) < 1e-8 ); 

switch caseBoundaryCondition
    case 1
        C = [topNodes zeros( length( topNodes ), 1 ); ...
             leftNodes zeros( length( leftNodes ), 1 ); ...
             rigthNodes zeros( length( rigthNodes ), 1 ); ...
             bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 2
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 3
        C = [rigthNodes ones( length( rigthNodes ), 1 ); ...
             leftNodes zeros( length( leftNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 4
        C = [topNodes ones( length( topNodes ), 1 ); ...
             bottomNodes zeros( length( bottomNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 5
        SourcePoint = [.68 -.68];
        SourceStrength = 1;
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        dist1 = sqrt((X(:,1)-SourcePoint(1)).^2+(X(:,2)-SourcePoint(2)).^2);
        ind = find(dist1==min(dist1));
        C = [C; ind(1) SourceStrength];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 6
        left1 = find( (X(:,1) == x_lo) & ( 0.55 < X(:,2)) );
        left2 = find( (X(:,1) == x_lo) & (X(:,2) < 0.55) );
        C = [left1 ones( length( left1 ), 1 ); ...
             left2 2*ones( length( left2 ), 1 );
             rigthNodes zeros( length( rigthNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 7
        SourcePoint = [.6878 -.751];
        SourceStrength = 1;
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        dist1 = sqrt((X(:,1)-SourcePoint(1)).^2+(X(:,2)-SourcePoint(2)).^2);
        ind = find(dist1==min(dist1));
        C = [C; ind(1) SourceStrength];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 8
        SourcePoint = [.86 -.67];
        SourceStrength = 1;
        C = [bottomNodes zeros( length( bottomNodes ), 1 )];
        dist1 = sqrt((X(:,1)-SourcePoint(1)).^2+(X(:,2)-SourcePoint(2)).^2);
        ind = find(dist1==min(dist1));
        C = [C; ind(1) SourceStrength];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
    case 9
        C = [rigthNodes ones( length( rigthNodes ), 1 ); ...
             leftNodes zeros( length( leftNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
        
    case 10
        C = [outNodes zeros( length( outNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);

    case 11
        C = [rigthNodes ones( length( rigthNodes ), 1 ); ...
             leftNodes ones( length( leftNodes ), 1 )];
        %
        nBC = size( C, 1 );
        Abc = zeros( nBC, numberOfUnknows );
        Abc(:,C(:,1)) = eye( nBC );
        bbc = C(:,2);
end;
