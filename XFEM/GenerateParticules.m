function [ x0, y0, vx0, vy0 ] = GenerateParticules( time, Particles )
% GENERATEPARTICULES to define the initial position and velocity of the
% particules
%
%  syntax: [ x0, y0, theta0, v0 ] = GenerateParticules( it, time, n0 );
%
%  it:       time step [integer]
%  time:     time options [structured array]
%  n0:       number of drops of water to consider in each particle launched
%
%  x0, y0: coordinates of the particules in initial position
%  theta0: angle of the velocity with the horizontal axis
%  v0: norm of the velocity of the particules (constant for all particles)

% R. Cottereau 04/2008

% define leak
[ n, v0, S ] = defineLeak( time.initial, time.t, Particles );

% definition of the leaking area, centered on theta = 157º
centerTube = [.9 -.7];
%centerTube = [.7 -.76];
%centerTube = [.698 -.7525];
radiusTube = 0.05;
%radiusTube = 0.015;
anglePositionLeak = 3/4*pi;
%centerTube = [.69 -.69];
%radiusTube = 0.03;
%anglePositionLeak = 153/180*pi;
openingPositionLeak = sqrt( S / pi ) / radiusTube;
angleLeakCenter = 3/4*pi;

% definition of the initial conditions
%r0 = randn(n,1);
%r0( ( abs(r0) > 1 ) ) = 0;
r0 = rand(n,1) - 0.5;
%r = randn(n,1)*1/3;
r = rand(n,1) - 0.5;
theta0 = angleLeakCenter + r0*Particles.AngleVelocity;
theta = anglePositionLeak + r*openingPositionLeak;
x0 = centerTube(1) + radiusTube * cos(theta);
y0 = centerTube(2) + radiusTube * sin(theta);
vx0 = v0 * cos(theta0);
vy0 = v0 * sin(theta0);

