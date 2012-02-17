function [ n, v0, S, aS, bS, Qexit ] = defineLeak( t0, t, Particles )
% DEFINELEAK to define the characteristics of the leak in the water tube
% depending on time
%
% syntax [ n, v0, S, aS, bS, Qexit ] = defineLeak( t0, t, n0 );
%
%  t0:     initial time of the experience [scalar]
%  t:      current time of the experience [scalar]
%  n0:       number of drops of water to consider in each particle launched
%
%  n:      number of particles to be generated for that time step [integer]
%  v0:     initial velocity of the particles of water [scalar]
%  S:      leaking area [scalar]
%  aS, bS: define the equation of the evolution of the diameter of the leak
%          [scalars]
%  Qexit:  outgoing flux [scalar]
%
%  v0 was given by Gas Natural

% R. Cottereau 04/2008

% test
% if t > Particles.TimeForFinalHole
%     warning('overtime?')
% end

% number of particules to be generated
n = Particles.ParticlesPerTimeStep;

% initial velocity of the particles (given by Gas Natural)
v0 = 34.64;

% leaking area
d0 = .4e-3;
df = 10e-3;
tf = Particles.TimeForFinalHole;
aS = (df-d0)/(tf-t0);
bS = - aS * t0 + d0;
d = aS * t + bS;
S = pi * d^2 / 4;

% outgoing flux
Qexit = S * v0;
