function [ vx, vy ] = velocityParticules( x0, vx0, vy0, xi );
% VELOCITYPARTICULES to compute the velocity of particules at xi
%
%  syntax: [ vx, vy ] = velocityParticules( x0, y0, vx0, vy0, xi );
%
%  x0,y0:   coordinates of particules at initial position
%  vx0,vy0: velocity of particules at initial position
%  xi:      coordinate of the points where the velocity is sought
%
%  All the vectors should have the same number of lines, in which case the
%  velocity is searched in one position for each particule, or their should
%  be only one particule, the velocity of which is searched in different
%  positions

% R. Cottereau 04/2008

% constants
g = 9.81;

vx = vx0;
vy = vy0 - g * ( (xi-x0) ./ vx );
