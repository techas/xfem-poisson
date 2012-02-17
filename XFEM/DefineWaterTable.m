function phreaLevel = DefineWaterTable( t0, t, waterTable, particles )
% DEFINEWATERTABLE to define at a given time the level of the water table
%
%  syntax: phreaLevel = DefineWaterTable( it, t, waterTable );
%
%  it:         time step index [integer]
%  t:          time options [structured array]
%  waterTable: water table options [structured array]
%
%  phreaLevel: water table level at time index it
%
%  the equation that it solved is
%    (Qleak - Qexit) * dt = dz
%  with Qleak = A * t^2 + B * t + C ;
%  and  Qexit = D * (z - zRef)
%
%  resulting in the differential equation
%    dz/dt + D * z = Qleak + D * zRef
%  with homogeneous solution z(t) = exp( -D * t )
%  and non-homogeneous solutions z0(t) = ( Qleak + C * zRef - A / C ) / C
%
%  the initial condition is given by z(time.initial) = waterTable.Initial
%
%  the constant of the exponential is chosen such that the critical level
%  of the water table waterTable.CriticalLevel is reached when the hole in
%  the water tube has a diameter of 10mm

% R. Cottereau 04/2008

% constants
D = waterTable.Constant;
zRef = waterTable.Reference;
h0 = waterTable.Initial;
C = waterTable.CriticalLevel;

% case of ever wet ground
if C==-Inf
    phreaLevel = Inf;

% case of ever dry ground
elseif C==Inf
    phreaLevel = -Inf;

% all other cases
else
    % particular solution
    [ n, v0, S, aS, bS ] = defineLeak( t0, t, particles );
    a = v0 * pi / 4 * [aS^2 2*aS*bS bS^2] + [0 0 D*zRef];
    a = [D 0 0;2 D 0;0 1 D] \ a';
    z0 = [t^2 t 1]*a;

    % constant for the homogeneous solution
    b = ( h0 - [t0^2 t0 1]*a ) * exp( D * t0 );

    phreaLevel = b * exp( - D * t) + z0;

    % the phreatic level is blocked on top by the pavement
    phreaLevel = min( phreaLevel, 0 );

end