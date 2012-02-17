function [Seg,SegBnd] = CrossedSegments( T, indT, levelSet, tol )
% CROSSEDSEGMENTS to get a list of segments crossed by the levelSet
%
% syntax: [Seg,SegBoundary] = CrossedSegments( T, indT, levelSet );
%
%  T:           connectivity matrix [Ne*e matrix]
%  indT:        indices of T on which the search should be performed. If
%               empty, all enriched elements are used. [Nt*1 vector]
%  levelSet:    nodal values of the level Set function [Nn*1 vector]
%
%  Seg:         list of the segments crossed by the levelSet (not ordered)
%               [Ns*2 matrix]. When the level set passes exactly above a
%               node, that node is used twice to create the segment.
%  SegBoundary: list of segments in Seg that are physically on the boundary
%               of the mesh (that are contained in only one element of T)
%               [Nb*2 matrix]

% R. Cottereau 04/2008

% when indT is empty, consider all enriched elements
if isempty(indT)
    type = classifyElements( levelSet, T, 0 );
    indT = find( type > 0 );
end

% constants
Nt = length( indT );
e = size( T, 2 );

% initializations
Seg = zeros( 2*Nt, 2 );
Te = T( indT, : );

% get the two segments that are crossed on each enriched element, and put 0
% when the element is not crossed
for i1 = 1:Nt
    
    N = Te( i1, : );
    SignLS = (levelSet(N)>tol)-(levelSet(N)<-tol);
    ProdSign = prod( SignLS );
    SumSign = sum( SignLS );
    indSeg = (i1-1)*2 + (1:2);
    
    % case when at least one node is exactly crossed by the level set
    if ProdSign == 0
        i2 = find( SignLS == 0 );
        N1 = N( 1, i2 );
        N2 = N( 1, setdiff( 1:e, i2 ) );
        % case when the level set follows one side (two nodes are crossed)
        if length( i2 ) == 2
%            Seg( indSeg, : ) = [ N1 ; N1 ]';
%            Seg( indSeg, : ) = [ N1 N1 ]';
%            Seg( indSeg, : ) = [ N1 ; N1 ];
            Seg( indSeg, : ) = [ 0 0 ; 0 0 ];
        % case when exactly one node is crossed
        elseif length( i2 ) == 1
            if SumSign == 0
                Seg( indSeg, : ) = [ N1 N1 ; N2 ];
            elseif abs(SumSign) == 2
                Seg( indSeg, : ) = [ 0 0 ; 0 0 ];
            else
                error('i do not understand that ... ')
            end
        else
            error(['the case where all nodes of an element are crossed' ...
                   ' exactly by the level set is not possible']);
        end
        
    % case when the level set does not pass on top of any node
    else
        i2 = find( SignLS == -SumSign );
        i3 = setdiff( 1:e, i2 );
        Seg1 = sort( N( [i2 i3(1)] ) );
        Seg2 = sort( N( [i2 i3(2)] ) );
        [n11 n12] = size(Seg1);
        [n21 n22] = size(Seg2);
        if n11~=1 | n12~=2 | n21~=1 | n22~=2
            N
            SignLS
            levelSet( N )
            Seg1
            Seg2
        end
        Seg( indSeg, : ) = [ Seg1 ; Seg2 ];
    end

end

% get rid of the segments that are only zeros
i1 = ( Seg(:,1) ~= 0 | Seg(:,2) ~= 0 );
Seg = Seg( i1, : );

% get rid of the repeated segments, and select those that are not repeated
% the "trick" on j1 selects the elements of j1 that are not repeated
Seg = sort( Seg, 2 );
[Seg,i1,j1] = unique(Seg,'rows');
j1 = sort(j1);
test = diff( j1 );
test = ( [ test ; 1 ] + [ 1; test] ) == 2;
j1 = j1( test );
SegBnd = Seg( j1, : );


% % alternative for the Gas Natural case only
% P1 = X( Seg(:,1), : );
% P2 = X( Seg(:,2), : );
% x_lo = min( X(:,1) );
% x_up = max( X(:,1) );
% y_lo = min( X(:,2) );
% y_up = max( X(:,2) );
% centerTube1 = [.3 -.5];
% radiusTube1 = .03;
% centerTube2 = [.69 -.69];
% radiusTube2 = .03;
% bnd1hole1 = ( abs( ( P1(1) - centerTube1(1) ).^2 + ...
%                  ( P1(2) - centerTube1(2) ).^2 - radiusTube1^2 ) < 1e-16 );
% bnd2hole1 = ( abs( ( P2(1) - centerTube1(1) ).^2 + ...
%                  ( P2(2) - centerTube1(2) ).^2 - radiusTube1^2 ) < 1e-16 );
% bnd1hole2 = ( abs( ( P1(1) - centerTube2(1) ).^2 + ...
%                  ( P1(2) - centerTube2(2) ).^2 - radiusTube2^2 ) < 1e-16 );
% bnd2hole2 = ( abs( ( P2(1) - centerTube2(1) ).^2 + ...
%                  ( P2(2) - centerTube2(2) ).^2 - radiusTube2^2 ) < 1e-16 );
% bnd1ext = ( IsOnBoundary( X, P1 ) ~= 0 );
% bnd2ext = ( IsOnBoundary( X, P2 ) ~= 0 );
% bnd = ( bnd1hole1 | bnd1hole2 | bnd1ext ) & ...
%                                        ( bnd2hole1 | bnd2hole2 | bnd2ext );
% SegBnd = Seg( find( bnd ), : )
