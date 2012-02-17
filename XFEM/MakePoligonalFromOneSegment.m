function [poli2,ListSeg,SegEnd,UsedElt,poli1]=MakePoligonalFromOneSegment( ...
                               X, T, SegBegin, SegEnd, ListSeg, levelSet )
% Starting from one segment, construct the poligonal as a sequence of
% contiguous segments
%
% syntax: [poli1,ListSeg,SegEnd]=MakePoligonalFromOneSegment( ...
%                              X, T, SegBegin, SegEnd, ListSeg, levelSet )
%
% X,T:      nodal coordinates and connectivity matrix
% SegBegin: is the first segment [1*2 vector]
% SegEnd:   is a set of possible final segments [n1*2 matrix]
% ListSeg:  is the list of segments that the level Set can cross 
%           [n2*2 vector]
% levelSet: nodal values of the level Set function
%
% poli1:   ordered list of segments composing the poligonal
% ListSeg: (in output) same as input without the segments that have been 
%          used to construct the poligonal have been taken out.
% SegEnd:  (in output) same as input but without the segment that was 
%          actually used to close the segment list
% UsedElt: elements that are crossed by the segments, in appropriate order

% R. Cottereau 04/2008

% constants
nEnd = size(SegEnd,1);
NbIterMax = 10000;

% initializations
poli1 = SegBegin;
Seg1 = SegBegin;
UsedElt = [];

% if SegBegin is present in ListSeg, take it out
[ i1, j1 ] = find( ListSeg == Seg1(1) );
[ i2, j1 ] = find( ListSeg(i1,:) == Seg1(2) );
if ~isempty(i2)
    nL = size( ListSeg, 1 );
    ListSeg = ListSeg( setdiff( 1:nL, i1(i2) ), : );
end

% loop on the available segments until final segment is reached
for i1=1:NbIterMax;
    [ Seg1, ListSeg, E ] = findNextSegment( Seg1, ListSeg, T );
    % when E=0, this means that no element was found. This is normal with
    % closed loops so that the loop must be closed by hand
    if E == 0
        Seg1 = SegBegin;
        E = InContactWithNodes( T, [ poli1( end, : ); Seg1 ] );
    end
    poli1 = [ poli1 ; Seg1 ];
    UsedElt = [ UsedElt ; E ];
    testSeg = repmat( Seg1, [nEnd 1] );
    if sum( sum( SegEnd == testSeg, 2 ) == 2 ) == 1
        break
    end
end

% test on the total number of iterations
if i1==NbIterMax;
    error('UPDATELEVELSET: maximum number of iterations reached');
end

% poli1
% % plot the segments
% figure(3);hold off;
% trimesh(T,X(:,1),X(:,2),'Color','g');hold on;
% plot( X(SegBegin,1), X(SegBegin,2), 'r-' );hold on;
% text( X(SegBegin(1),1), X(SegBegin(1),2), num2str(SegBegin(1)) );
% text( X(SegBegin(2),1), X(SegBegin(2),2), num2str(SegBegin(2)) );
% plot( X(SegEnd,1), X(SegEnd,2), 'b-' );hold on;
% text( X(SegEnd(1),1), X(SegEnd(1),2), num2str(SegEnd(1)) );
% text( X(SegEnd(2),1), X(SegEnd(2),2), num2str(SegEnd(2)) );
% for i1=1:size(poli1,1)
%     plot(X(poli1(i1,:),1),X(poli1(i1,:),2),'k--');hold on;
% text( X(poli1(i1,1),1), X(poli1(i1,1),2), num2str(poli1(i1,1)) );
% text( X(poli1(i1,2),1), X(poli1(i1,2),2), num2str(poli1(i1,2)) );
% end

% quit the last segment from the list of boundary segments
ind = find( SegEnd(:,1) == Seg1(1) & SegEnd(:,2) == Seg1(2) );
SegEnd = SegEnd( setdiff(1:nEnd,ind) , : );

% replace the ends of each segment by the position of the crossing of the
% level set
poli2 = poli1;
for i1=1:size(poli1,1)
    poli2(i1,:) = intersection( X(poli1(i1,:),:), levelSet(poli1(i1,:)) );
end
