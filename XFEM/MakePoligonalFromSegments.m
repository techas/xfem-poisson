function [ polis, Ei, Listi ] = MakePoligonalFromSegments( X, T, ...
                                              SegsBnd, ListSeg, levelSet )
% MAKEPOLIGONALSFROMSEGMENTS to construct a set of polygonals from a list 
% of segments with a given list of boundaries
%
%  syntax: [ polis, Ei ] = MakePoligonalsFromSegments( X, T, ...
%                                             ListSeg, SegsBnd, levelSet );
%
%  X:        nodal coordinates [Nn*d matrix]
%  T:        connectivity matrix [Ne*e matrix]
%  ListSeg:  list of segments from which the poligonals can be created. Note
%            that the poligonals do not have to be closed [Nl*2 matrix]
%  SegsBnd:  list of segments from which the beginning and ending segments
%            of each polygonal must be picked [Nb*2 matrix]
%  levelSet: nodal values of the level set function [Nn*1 vector]
%
%  polis:    list of poligonals [Nb/2*1 cell array of Np*2 matrices]. Note
%            that the Np may be different for two poligonals and that they
%            must sum up to Nl+Nb
%  Ei:       list of elements crossed by the poligonals 
%            [Nb/2*1 cell array of Np*1 matrices]
%
%  To close the poligonals that are created with this function, one can use
%  CLOSEPOLIGONALS

% R. Cottereau 04/2008

% case when no boundaries are found (closed loop)
if isempty( SegsBnd )
    SegsBnd = ListSeg( [1 1], : );
    ListSeg = QuitSegFromList( SegsBnd(1,:), ListSeg );
end

% constants
Ns = size( SegsBnd, 1 ) / 2;

% initialization
polis = cell( Ns, 1 );
Ei = cell( Ns, 1 );
Listi = cell( Ns, 1 );

% loop on the available boundary elements
for i1 = 1:Ns

    % select the first segment
    Seg1 = SegsBnd( 1, : );

    % take Seg1 out of SegsBnd and ListSeg
    SegsBnd = QuitSegFromList( Seg1, SegsBnd );
    ListSeg = QuitSegFromList( Seg1, ListSeg );
    
    % case of a closed loop ... SegsBnd contains Seg1 as beginning and
    % ending segment, which both have been taken out in the previous lines
    % so that Seg1 should be added once to SegsBnd and ListSeg
    n1 = size( SegsBnd, 1 );
    if floor( n1/2 ) == n1/2;
        SegsBnd = [SegsBnd ; Seg1 ];
        ListSeg = [ListSeg ; Seg1 ];
    end
 
    % get the list of segments starting from Seg1
    [ poli1, ListSeg, SegsBnd, Ei1,tmp ] = MakePoligonalFromOneSegment( ...
                                  X, T, Seg1, SegsBnd, ListSeg, levelSet );
    polis{ i1 } = poli1;
    Ei{ i1 } = Ei1;
    Listi{ i1 } = tmp;

end

% check that all segments have been considered
if ~isempty( SegsBnd ) || ~isempty( ListSeg )
    SegsBnd
    ListSeg
    error('there are some segments left in SegsBnd or ListSeg');
end

%====================================
function List = QuitSegFromList( Seg, List )
% return the list of segments which is the same as List without the
% segments Seg
[ i1, j1 ] = find( List(:,1) ~= Seg(1) );
[ i2, j1 ] = find( List(:,2) ~= Seg(2) );
List = List( union( i1, i2 ), : );
