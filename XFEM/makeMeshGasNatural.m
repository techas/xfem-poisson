function makeMeshGasNatural

% R. Cottereau 03/2008

close all;
clear all;
clc;
addpath('./mesh2d');

nodeBox = [-.2 -1.2;
           1.2 -1.2;
           1.2 0 ;
           -.2 0 ];
edgeBox = [1 2;
           2 3;
           3 4;
           4 1];
    
centerTube1 = [.25 -.5];
radiusTube1 = .025;
numberNodesTube1 = 10;
thetaTube1 = 2*pi*[numberNodesTube1:-1:1]'/numberNodesTube1;
nodeTube1 = repmat(centerTube1,numberNodesTube1,1) + ...
                            radiusTube1*[cos(thetaTube1) sin(thetaTube1)];
edgeTube1 = [[1:numberNodesTube1]' [[2:numberNodesTube1]';1]];
                        
centerTube2 = [.7 -.84];
radiusTube2 = .05;
numberNodesTube2 = 28;
thetaTube2 = 2*pi*[numberNodesTube2:-1:1]'/numberNodesTube2;
nodeTube2 = repmat(centerTube2,numberNodesTube2,1) + ...
                            radiusTube2*[cos(thetaTube2) sin(thetaTube2)];
centerTube3 = [.7 -.76];
radiusTube3 = .015;
numberNodesTube3 = 10; % has to be even
thetaTube3 = 2*pi*[numberNodesTube3:-1:1]'/numberNodesTube3;
nodeTube3 = repmat(centerTube3,numberNodesTube3,1) + ...
                            radiusTube3*[cos(thetaTube3) sin(thetaTube3)];
X31 = nodeTube3(1,1);
Y31 = nodeTube3(1,2);
RefineBox23 = .008;
box23 = [X31 Y31
         X31 Y31-(centerTube3(2)-centerTube2(2))
         X31-2*radiusTube3 Y31-(centerTube3(2)-centerTube2(2))
         X31-2*radiusTube3 Y31];
tube2Inbox23 = inpolygon( nodeTube2(:,1), nodeTube2(:,2), ...
                                                  box23(:,1), box23(:,2) );
tube2Inbox23first = find( tube2Inbox23, 1, 'first' );
tube2Inbox23last = find( tube2Inbox23, 1, 'last' );
tube3Inbox23first = numberNodesTube3/2+1;
tube3Inbox23last = numberNodesTube3;
lbox23 = nodeTube3( tube3Inbox23first, 2)-nodeTube2( tube2Inbox23first, 2);
box23base = linspace( 0, lbox23, ceil(lbox23/RefineBox23) )';
nbox23 = length(box23base)-1;
box23up = repmat(nodeTube3( tube3Inbox23first, : ), [nbox23 1] ) - ...
          [ zeros( nbox23, 1 ) flipud( box23base(2:end,1) ) ];
box23down = repmat(nodeTube3( 1, : ), [nbox23 1] ) - ...
          [ zeros( nbox23, 1 ) box23base(2:end,1) ];
indTube3 = [tube3Inbox23first:tube3Inbox23last 1];

nodeTube2 = [nodeTube2( 1:tube2Inbox23first-1, : )
             box23up
             nodeTube3( indTube3, : )
             box23down
             nodeTube2( tube2Inbox23last+1:end, : )];
nnodeTube2 = size( nodeTube2, 1 );

edgeTube2 = [[1:nnodeTube2]' [[2:nnodeTube2]';1]];

node = [nodeBox;
        nodeTube1;
        nodeTube2];

edge = [edgeBox;
        size(edgeBox,1)+edgeTube1;
        size(edgeBox,1)+size(edgeTube1,1)+edgeTube2];
        
hdata.hmax = 0.02;
[X,T] = mesh2d(node,edge,hdata);

save meshGasNaturalVeryFine.mat X T;