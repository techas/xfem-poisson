function yObs = trajectoryParticle( x0, y0, vx0, vy0, xObs)
% TRAJECTORYPARTICLE to compute the value of the ordinate of the particule
% in an observation point
%
% syntax: yObs = trajectoryParticle( x0, y0, v0, th0, xObs)
% 
% X0: position of the particle at t0, with X0 = [x0 y0]
% V0: velocity of the particle at t0, with V0 = [vx0 vy0]
% xObs: x-coordinate of the points of observation
%
% Note that many particles can be considered at the same time, the number
% of lines of X0, V0 must then be the same.
%
% Another possibility is to have only one particular and several points of
% observation. X0 and V0 must then have only one line.
% 
% yObs: y-coordinate of the particle at the points of observation

% R. Cottereau 04/2008

% constants
g = 9.81;

% changing sizes when multiple observation points are asked for one
% particule only
if length(x0)==1
    I1 = ones( length(xObs), 1 );
    x0 = x0 * I1;
    y0 = y0 * I1;
    vx0 = vx0 * I1;
    vy0 = vy0 * I1;
end

x = (xObs-x0);
yObs = y0 + x.*(vy0./vx0) - (g/2)*(x./vx0).^2;

%yObs = y0 + (xObs-x0).*tan(th0) - (g/(2*v0^2))*((xObs-x0)./sin(th0)).^2;