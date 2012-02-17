%function [ Vi, Ei, Ni, li ] = DropSimulation( X, T, levelSet, ...
function [ Xi, Vi, Ei, Ni, li ] = DropSimulation( X, T, levelSet, ...
                                                  time, Particles, opts )
% DROPSIMULATION to simulate the projection of drops of water on a surface
% defined as the level set function
%
%  syntax: [Xi,Vi,Ei,Ni,li] = DropSimulation( X, T, levelSet, ...
%                                                         it, time, opts )
%
%  X, T:     nodal coordinates and connectivity matrix
%  levelSet: nodal values of the level set function
%  it:       time step [integer]
%  time:     time options [structured array]
%  opts:     options [structured array]
%
%  Xi:       positions of the particules at impact [n*2 matrix]
%  Vi:       velocities of the particules at impact [n*2 matrix]
%  Ei:       elements on which the impact takes place [n*1 vector]
%  Ni:       normal vector to the impact zone [n*2 matrix]
%  li:       length of the segment on which the impact takes place
%            [n*1 vector]

% R. Cottereau 04/2008

% initial conditions of the particules
[ x0, y0, vx0, vy0 ] = GenerateParticules( time, Particles );

% DrawParticuleImpacts( 1, X, T, levelSet, x0, y0, vx0, vy0);

%
% define the intersections of the trajectory of each particule with the 
% level Set
% [ Ei, nx, ny, li ] = intersectionTrajectory( x0, y0, vx0, vy0, ...
%                                                     X, T, levelSet, opts);
[ xi, yi, Ei, nx, ny, li ] = intersectionTrajectory( x0, y0, vx0, vy0, ...
                                                      X, T, levelSet, opts);

% compute the velocity of the particules at the points of intersection
[ vxi, vyi ] = velocityParticules( x0, vx0, vy0, xi );

% DrawParticuleImpacts( 3, X, T, levelSet, x0, y0, vx0, vy0, xi, yi, ...
%                                                               vxi, vyi);
% prepare output
Xi = [xi yi];
Vi = [vxi vyi];
Ni = [nx ny];

% modify normal vector when not in the same direction as the velocity
ind = find( 0 > vxi .* nx + vyi .* ny );
Ni( ind ) = -Ni( ind );
     