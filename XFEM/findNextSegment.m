function [ Seg2, ListSeg, E ] = findNextSegment( Seg1, ListSeg, T );
% FINDNEXTSEG to find, in a list of segments crossed by the level set, the 
% one following Seg1
%
% syntax: [ Seg2, ListSeg, E ] = findNextSegment( Seg1, ListSeg, T );
%
%  Seg1:     starting segment [1*2 vector]
%  ListSeg:  available list of segments, from which to choose the following
%            one (NB: Seg1 should not be contained in that list)
%            [n*2 matrix]
%  T:        connectivity matrix of the mesh [Ne*e matrix]
%
%  Seg2:     segment following Seg1 in the list ListSeg [1*2 vector]
%  ListSeg:  (in output) same as in input, but withoug Seg2
%            [(n-1)*2 matrix]
%  E:        element that contains both Seg1 and Seg2 [integer]

% R. Cottereau 04/2008

% find the elements in contact with Seg1, and their nodes
% (if everything goes well, N should contain two nodes when the segments is
% in the middle of the mesh, and one when on the boundary)
[ i1, j1 ] = find( T == Seg1(1) );
[ i2, j1 ] = find( T(i1,:) == Seg1(2) );
indE = i1( i2 );
E = T( indE, : );
N = setdiff( E(:), Seg1 );

% among the possible segments, select those that contain one of the nodes
% of Seg1
[ i1, j1 ] = find( ListSeg == Seg1(1) );
[ i2, j1 ] = find( ListSeg == Seg1(2) );
indS = union( i1, i2 );
Seg2 = ListSeg( indS, : );

% select those that contain also one of the nodes of N (if all
% goes well there should be only one)
[ i1, j1 ] = find( Seg2 == N(1) );
if isempty(i1);
    [ i1, j1 ] = find( Seg2 == N(2) );
    N = N(2);
else
    N = N(1);
end
Seg2 = Seg2( i1, : );

% in case the level set crosses exactly on one node, Seg2 is probably
% empty. Select then the segments that represent one point exactly and see
% if they are in N. If both elements in N are repeated, then the problem is
% more complex ...
% another possibility is that Seg1 be a single node ... in that case,
% search all the elements in E for a segment that is in ListSeg ... and
% hope there is only one
if isempty( Seg2 )
    N = setdiff( E(:), Seg1 );
    
    % case Seg1 = one node
    if length(N)~=2
        if Seg1(1) ~= Seg1(2)
            error('for that case Seg1 should be a single node');
        end
        N1 = Seg1(1);
        Seg2 = [];
        [ i2, j1 ] = find( E == N1 );
        for i1 = 1:size(E,1)
            Segi = E( i1, setdiff( 1:size(T,2), j1(i1) ) );
            [ i2, j2 ] = find( ListSeg == Segi(1) );
            [ i3, j2 ] = find( ListSeg(i2,:) == Segi(2) );
            Seg2 = [ Seg2 ;ListSeg( i2(i3), : ) ];
        end
       if size(Seg2,1)~=1;
            error('i do not know what to do ...');
       end
        
    % case Seg2 = one node
    else
        i1 = find( ListSeg(:,1) == ListSeg(:,2) );
        [i2,j1] = find( ListSeg(i1,:) == N(1) );
        [i3,j1] = find( ListSeg(i1,:) == N(2) );
        if isempty(i2) & isempty(i3)
            Seg2 = [];
            E = 0;
            return;
        elseif isempty(i2)
            Seg2 = ListSeg( i1(i3(1)), : );
            N = N(2);
        elseif isempty(i3)
            Seg2 = ListSeg( i1(i2(1)), : );
            N = N(1);
        else
            error('what to do ??? ...');
        end
    end
end
   
% take Seg2 out of ListSeg
[ i1, j1 ] = find( ListSeg(:,1) ~= Seg2(1) | ListSeg(:,2) ~= Seg2(2) );
ListSeg = ListSeg( i1, : );

% get the element that contains both segments Seg1 and Seg2
if Seg1(1) ~= Seg1(2) 
    [ i1, j1 ] = find( E == N );
else
    [ i1, j1 ] = find( E == Seg1(1) );
    [ i2, j2 ] = find( E == Seg2(1) );
    [ i3, j3 ] = find( E == Seg2(2) );
    i1 = intersect(intersect(i1,i2),i3);
%     [ i1, j1 ] = find( E == N(1) );
%     [ i2, j1 ] = find( E(i1,:) == N(2) );
%     i1 = i1(i2);
end
E = indE(i1);

