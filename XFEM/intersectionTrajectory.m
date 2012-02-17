% function [ Ei, nx, ny, li] = intersectionTrajectory( x0, y0, ...
function [ xi, yi, Ei, nx, ny, li] = intersectionTrajectory( x0, y0, ...
                                          vx0, vy0, X, T, levelSet, opts )
% INTERSECTIONTRAJECTORY to compute the intersection points for a series of
% particles between their trajectories and the level Set
%
%  syntax: [ x, y, Ei ] = intersectionTrajectory( x0, y0, vx0, vy0, ...
%                                                   X, T, levelSet, opts )
%
%  x0, y0:   coordinates of the particles in initial position
%            [n*1 vectors]
%  theta0:   angle of the velocity with the horizontal axis [n*1 vector]        
%  v0:       norm of the velocity of the particles (constant for all 
%            particles) [scalar]
%  X, T:     nodes and connectivity matrix of the mesh 
%            [Nnode*2 and Nelt*· matrices]
%  levelSet: nodal values of the level Set function [Nelt*1 vector]
%  opts:     options [structured array]
%
%  xi, yi:   coordinates of the intersection, for each particle, of its
%            trajectory with the level Set function [n*1 vectors]
%  Ei:       list of elements on which each impact takes place [n*1 vector]
%  nx, ny:   normal vector to the level Set at the impact points 
%            [n*1 vectors]
%  li:       length of the segment on which the impact takes place

% R. Cottereau 04/2008

% constants
n = size( x0, 1 );

% define segments cut by the level Set and corresponding control points
% where the level Set function cancels
type = classifyElements( levelSet, T, 0 );
enrichedElements = find( type > 0 );
[ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet );
[ xCut, Ei ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, levelSet );
[ xCut, Ei ] = ClosePoligonals( xCut, Ei, X, levelSet );

xCut = xCut{ 1 };
Ei = Ei{ 1 };

% trajectory of the particles
% y = y0 + tan(th0)*(x-x0) - g*(x-x0)^2/[2*(v0*sin(th0))^2]
nX = size( xCut, 1 );
yObs = zeros( n, nX );
for i1 = 1:nX
    yObs(:,i1) = trajectoryParticle( x0, y0, vx0, vy0, xCut(i1,1) ) - xCut(i1,2);
end

% get impact positions and characteristics
[ Impact, Ei, li ] = GetImpact( yObs, xCut, Ei );


% compute the intersection of the segment with the trajectory
% ax+b = y0 + tan(th0)*(x-x0) - g*(x-x0)^2/[2*(v0*sin(th0))^2]
% and the normal vector to the segment at each intersection
%keyboard 
%[ nx, ny ] = intersectionTrajectorySegment( x0, y0, vx0, vy0, ...
%                                                           xCut, Impact );
[ xi, yi, nx, ny ] = intersectionTrajectorySegment( x0, y0, vx0, vy0, ...
                                                           xCut, Impact );

% compute the normal vector

% % plot the segments
% figure(4);hold off;
% trimesh(T,X(:,1),X(:,2),'Color','k');hold on;
% scatter(xCut(:,1),xCut(:,2),50,'b','filled');hold on;
% scatter( xi, yi,50,'r');hold on;
% ind = find( Ei ~= 0 );
% e = T(Ei(ind),:);e=e(:);
% scatter(X(e,1),X(e,2),50,'k','filled');hold on;

%==========================================================================
%function [nx,ny] = intersectionTrajectorySegment( x0, y0, vx0, vy0, ...
function [x,y,nx,ny] = intersectionTrajectorySegment( x0, y0, vx0, vy0, ...
                                                            xObs, Impact )
% INTERSECTIONTRAJECTORYSEGMENT to compute the point intersection between a
% segment (given by xObs and Impact) and a trajectory (defined by the other
% inputs)
%
%  syntax:   IntSect = intersectionTrajectory(x0,y0,v0,th0,xObs,Impact); 
%
% X0:     position of the particule at t0, with X0 = [x0 y0]
% V0:     velocity of the particule at t0, with V0 = v0 [cos(th0) sin(th0)]
% xObs:   coordinates of the points of observation
% Impact: segments of impact, whose vertices are given as indices of xObs
%
% x,y:    coordinates of the impact points
% nx,ny:  normal vector to the segment at the impact points

% R. Cottereau 04/2008

% constants
g = 9.81;
n = size(x0,1);

% define the segments to be computed and for each particule, to what
% segment it corresponds (SegPart)
Segments = unique( Impact(:) );
Segments = [Segments(1:end-1) Segments(2:end)];
nSeg = size( Segments, 1 );
ind = zeros( nSeg, 1 );
for i1 = 1:nSeg
    ind( i1 ) = ~isempty( ( Impact(:,1)==Segments(i1,1) ) );
end
Segments = Segments( find(ind), :);
nSeg = size( Segments, 1 );
SegPart = zeros(n,1);
for i1 = 1:nSeg
    SegPart( ( Impact(:,1)==Segments(i1,1) ) ) = i1;
end

% compute the parameters of the segments
% y = ax + b
% d = [b;a];
n2 = 2 * nSeg;
In = speye(nSeg);
A = sparse(n2,n2);
A(1:2:n2,1:2:n2) = In;
A(2:2:n2,1:2:n2) = In;
A(1:2:n2,2:2:n2) = spdiags( xObs(Segments(:,1),1), 0, nSeg, nSeg );
A(2:2:n2,2:2:n2) = spdiags( xObs(Segments(:,2),1), 0, nSeg, nSeg );
d = reshape( A \ xObs(Segments',2) , 2, nSeg )';

% compute the normal to the segments
nx = ones( size(d,1), 1 );
ny = zeros( size(d,1), 1 );
ind = find( d(:,2)~=0 );
normN = sqrt( 1 + d(:,2).^2 );
nx( ind ) = -d( ind, 2 ) ./ normN( ind );
ny( ind ) = 1 ./ normN( ind );

% assign the values of the parameters to each particule
d = d( SegPart, : );
nx = nx( SegPart );
ny = ny( SegPart );

% resolution of second order equation to find the x-coordinate of the
% intersection
% a*(x-x0)^2 + b*(x-x0) + c = 0
a = (g/2) ./ (vx0.^2);
b = d(:,2) - vy0./vx0;
c = d(:,1) + d(:,2).*x0 - y0;
delta = b.^2 - 4.*a.*c;
x = x0 + (-b + sqrt(delta)) ./ (2*a);
x2 = x0 + (-b - sqrt(delta)) ./ (2*a);
% the following choice is only valid for the Gas Natural case
ind = find( ((x2 > x) & (x0 > x2)) | (x > x0));
x(ind) = x2(ind);

% corresponding y-coordinate
y = d(:,1) + d(:,2).*x;
if ~isreal(y)
    delta
 %   x
    error('stop')
end

% xx = [0:0.01:1];Nx = length(xx);
% yy1 = d(:,1)*ones(1,Nx) + d(:,2)*xx;
% figure(2);hold off
% plot(xx,yy1,'r')
% hold on;
% scatter(xObs(Segments(:,1),1),xObs(Segments(:,1),2),50,'r')
% scatter(xObs(Segments(end,2),1),xObs(Segments(end,2),2),50,'r')
% set(gca,'Xlim',[.5 .7],'YLim',[-.7 -.5])
% figure(2);hold on;
% scatter(x,y,50,'g')
%==========================================================================
function  [ Impact, Ei, li ] = GetImpact( yObs, xCut, Ei )
% GETIMPACT to compute the position of the impacts of the particles on the
% level set
%
%  syntax: [ Impact, Ei, li ] = GetImpact( xCut, yObs, Ei );
%
%  xCut:   positions of the control points of the level set function
%          [Nc*2 matrix]
%  yObs:   difference between the altitude of the particles at each of the 
%          points in xCut, and the altitude of the point in xCut
%          [Np*Nc matrix]
%  Ei:     element corresponding to the intervals between two points in  
%          xCut [(Nc-1)*1 vector]
%
%  Impact: position of the impact for each particle [Np*2 matrix]
%  Ei:     element on which the impact takes place for each particle
%          [Np*1 vector]
%  li:     length of the segment on which the impact takes place
%          [Np*1 vector]

% R. Cottereau 04/2008

% constants
[ n, nX ] = size( yObs );
nX = nX - 1;

% order the points in xCut such that the first node is below the position
% of the initial particle
% Note that is probably valid only for some specific cases (of interest for
% the Gas Natural case)
dir1 = sign( yObs( :, 1 ) );
dir2 = sign( yObs( :, end ) );
if any( (dir1.*dir2) ~= -1 )
    error('unknown case')
end
if any(dir1 > 0 )
    warning('maybe .. should check that')
    yObs = yObs( :, end:-1:1 );
    xCut = xCut( end:-1:1, : );
    Ei = Ei( end:-1:1 );
end

% get the number of crossings of the level set for each particle
Impact0 = diff( ( yObs > 0 ), 1, 2 );
%NbCrossings = sum( abs( Impact0 ), 2 );

% once the test is passed we suppose that the segments are ordered such 
% that the first one is below the position of the initial particle
Impact = zeros( n, 1 );
RestPart = 1:n;
for i1=nX-1:-1:1
    ind = RestPart( ( Impact0( RestPart, i1 ) == 1 ) );
    Impact( ind ) = i1;
    RestPart = setdiff( RestPart, ind );
end
Impact = [ Impact Impact+1 ];

% compute the elements in which the intersection takes place
Ei = Ei( Impact(:,1) );

% compute the length of the segment on which the impact takes place
li = xCut( Impact(:,2), : ) - xCut( Impact(:,1), : );
li = sqrt( li(:,1).^2 + li(:,2).^2 );

