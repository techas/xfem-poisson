close all
clear all
clc
% computation of the XFEM and exact solutions for a 1D Poisson problem and
% study of the error in fluxes
%
%  Description of the problem
%             -> f                   -> f 
%  u=0 ___________________|_______________________ Ku'=F
%    x=0                 x=a                    x=L
%
%  equation of the problem
%                K . u'' = f

% R. Cottereau 05/2008

% parameters of the problem
K1 = 1;               % conductivity for the 1st domain
K2 = 1;              % conductivity for the 2nd domain
g = K2 / K1;          % ratio of conductivities
L = 1;                % total length of the domain
dx = .2;              % space step for the discretization
x = [0:dx:L]';        % discretization of the space
nx = length( x );     % number of nodes
N = nx + 2 - 1;       % number of unknowns (2 XFEM unknowns - 1 imposed u )
nT = nx - 1;          % number of elements;
f0 = 0;               % volume load on the entire length of the rod
F0 = 0;               % imposed flux at the end of the rod
fa = 1;               % volume load on the second domain only
uL = 0;               % imposed displacement at the end of the rod
aloc = 0.99;           % local position in an element of the discontinuity 
a = (floor(nT/2)+aloc)*dx;% global position of the discontinuity
ar = a / L;            % global relative position of the discontinuity
E1 = [ 1 floor(nT/2) ];   % indices of the elements of the 1st domain
n1 = E1(2) - E1(1) + 1;   % number of elements in the 1st domain
Ea = ceil(nT/2);          % index of the element with discontinuity
E2 = [ floor(nT/2)+2 nT]; % indices of the elements of the 2nd domain
n2 = E2(2) - E2(1) + 1;   % number of elements in the 1st domain
K = [ K1*ones(n1,1)
      K1*aloc + K2*(1-aloc)
      K2*ones(n2,1) ];    % vector of conductivity per node
xu = [ x(1:n1+1); a; x(n1+2:end) ];  % vector for evaluation of u
nxu = length( xu );       % length of vector xu
xs = zeros( 2*nxu-2, 1 );  % vector for evaluation of fluxes
xs( 1:2:end ) = xu(1:end-1) + 1e-8;
xs( 2:2:end ) = xu(2:end) - 1e-8;
nxs = length( xs );       % length of vector xs

%==========================================================================
% EXACT SOLUTION OF THE PROBLEM
%==========================================================================
Ue = zeros( nxu, 1 );
ind1 = 1 : n1+2;
ind2 = n1+2 : nxu;
Se = zeros( nxs, 1 );
ind1s = 1 : 2*n1+2;
ind2s = 2*n1+3 : nxs;

if uL~=0
    % only for Fa == 0
    if Fa~=0;error('Fa neq 0');end;
    d = ( 1 - ar * ( 1 - g ) ) * L;
    Ue( ind1 ) = f0/(2*K1) * xu( ind1 ) .* ...
             ( xu( ind1 ) - (1-ar^2*(1-g))/d ) + uL*g/d * xu( ind1 );
    Ue( ind2 ) = f0*L/(2*K2) * ( 1- xu( ind2 )/L ) .* ...
                  ( - xu( ind2 ) - ar*(1-aloc)*(1-g)/d ) + ...
                                    uL/d * ( xu( ind2 ) - ar*(1-g)*L );
    Se( ind1s ) = f0/2 * ( 2*xs( ind1s ) - (1-ar^2*(1-g))/d ) + K2*uL/d;
    Se( ind2s ) = f0/2 * ( 2*xs( ind2s ) - (1-ar^2*(1-g))/d ) + K2*uL/d;
elseif fa~=0
    % only for K1==K2
    if K1~=K2;error('K1 neq K2!!');end
    Ue( ind1 ) = -fa/(2*K1) * (1-a) * xu( ind1 );
    Ue( ind2 ) = fa*L^2/(2*K2) * ( (1-xu( ind2 )/L).^2 -(1-a^2) );
    Se( ind1s ) = -fa*L/K1*(1-a);
    Se( ind2s ) = -fa*L/K1 * ( 1 - xs( ind2s )/L );
else
    Ue( ind1 ) = f0/(2*K1) * xu( ind1 ).^2 + (F0-f0*L)/K1 * xu( ind1 );
    Ue( ind2 ) = f0/(2*K2) * xu( ind2 ).^2 + (F0-f0*L)/K2 * xu( ind2 ) ...
                                       + a*(F0+f0*(a/2-L))*(K2-K1)/(K1*K2);
    Se( ind1s ) = f0 * xs( ind1s ) + F0-f0*L;
    Se( ind2s ) = f0 * xs( ind2s ) + F0-f0*L;
end

%==========================================================================
% APPROXIMATE FEM SOLUTION OF THE PROBLEM
%==========================================================================
% construction of the mass matrix
M = zeros( N-2 );

%---mass matrix for 1st domain
M1 = 2 * K1 * eye( n1 );
M1( n1, n1 ) = 0;
m1 = - K1 * ones( n1-1, 1 );
M1 = M1 + diag( m1, 1 ) +  diag( m1, -1 );
ind1 = 1 : n1;
M( ind1, ind1 ) = M1;

%---mass matrix for 2nd domain
M2 = 2 * K2 * eye( n2+1 );
M2( 1, 1 ) = 0;
M2( n2+1, n2+1 ) = K2;
m2 = - K2 * ones( n2, 1 );
M2 = M2 + diag( m2, 1 ) +  diag( m2, -1 );
ind2 = n1+1 : n1+n2+1;
M( ind2, ind2 ) = M2;

%---mass matrix for the enriched elements
m11 = (1+aloc) * K1 + (1-aloc) * K2;
m12 = -aloc * K1 - (1-aloc) * K2;
m22 = aloc * K1 + (2-aloc) * K2;
Me = [ m11  m12
       m12  m22];
inde = [ n1 n1+1 ];
M( inde, inde ) = Me;

% construction of the load vector
F = zeros( N-2, 1 );
F( [ 1 n1+n2+1 ] ) = - f0 * dx / 2;
F( 2 : n1+n2 ) = - f0 * dx;
F( n1+n2+1 ) = F( n1+n2+1 ) + F0;

% adding boundary conditions
M = M / dx;
if uL~=0
    BcM = [ zeros(nx-2,1) ; 1 ];
    BcF = uL;
    M = [ M BcM ; BcM' 0 ];
    F = [ F ; BcF ];
elseif fa~=0
    F( n1 ) = F( n1 ) -(1-a)^2/4 * dx * fa;
    F( n1+1 ) = F( n1+1 ) -(1-a^2)/4 * dx * fa;
    F( n1+1 : n1+n2 ) = F( n1+2 : n1+n2 ) - fa * dx;
    F( n1+n2+1 ) = F( n1+n2+1 ) - fa * dx / 2;
end

% resolution of the system
X = M \ F;

% construction of the displacement
Uh = zeros( size( xu ) );
Uh( 2:n1+1 ) = X( 1:n1 );
Uh( n1+2 ) = X( n1 ) * (1-aloc) + X( n1+1 ) * aloc;
Uh( n1+3:end ) = X( n1+1:n1+n2+1 );

% construction of the flux
Sh = zeros( size( xs ) );
Sh( 1:2:2*n1 ) = K1 * ( X( 1:n1 ) - [ 0 ; X(1:n1-1) ] );
Sh( 2:2:2*n1 ) = K1 * ( X( 1:n1 ) - [ 0 ; X(1:n1-1) ] );
Sh( 2*n1+[1:2] ) = K1 * ( X( n1+1 ) - X( n1 ) );
Sh( 2*n1+[3:4] ) = K2 * ( X( n1+1 ) - X( n1 ) );
Sh( 2*n1+5:2:end ) = K2 * ( X( n1+2:n1+n2+1 ) - X(n1+1:n1+n2) );
Sh( 2*n1+6:2:end ) = K2 * ( X( n1+2:n1+n2+1 ) - X(n1+1:n1+n2) );
Sh = Sh / dx;

%==========================================================================
% APPROXIMATE XFEM SOLUTION OF THE PROBLEM
%==========================================================================
% construction of the mass matrix
M = zeros( N );

%---mass matrix for 1st domain
M1 = 2*K1 * eye( n1 );
M1( n1, n1 ) = 0;
m1 = -K1 * ones( n1-1, 1 );
M1 = M1 + diag( m1, 1 ) +  diag( m1, -1 );
ind1 = 1 : n1;
M( ind1, ind1 ) = M1;

%---mass matrix for 2nd domain
M2 = 2*K2 * eye( n2+1 );
M2( 1, 1 ) = 0;
M2( n2+1, n2+1 ) = K2;
m2 = -K2 * ones( n2, 1 );
M2 = M2 + diag( m2, 1 ) +  diag( m2, -1 );
ind2 = n1+1 : n1+n2+1;
M( ind2, ind2 ) = M2;

%---mass matrix for the enriched elements
m11 = (1+aloc) * K1 + (1-aloc) * K2;
m12 = -aloc * K1 - (1-aloc) * K2;
m22 = aloc * K1 + (2-aloc) * K2;
m1e1 = -2 * aloc * (1-aloc) * ( K1 - K2 );
m1e2 = -2 * aloc^2 * ( K1 - K2 );
m2e1 = 2 * aloc * (1-aloc) * ( K1 - K2 );
m2e2 = 2 * aloc^2 * ( K1 - K2 );
me1e1 = 2*(1-(1-2*aloc)^3)/3 * K1 + 16*aloc^2*(1-aloc)/3 * K2;
me1e2 = 8*aloc^2*(3-4*aloc)/6 * K1 - ...
                   8*aloc^2*(1-6*aloc+9*aloc^2-4*aloc^3)/6/(1-aloc)^2 * K2;
me2e2 = 16*aloc^3/3 * K1 + 4*aloc^2*(1-2*aloc+4*aloc^2)/(3*(1-aloc)) * K2;
Me = [ m11  m12  m1e1  m1e2
       m12  m22  m2e1  m2e2
       m1e1 m2e1 me1e1 me1e2
       m1e2 m2e2 me1e2 me2e2];
inde = [ n1 n1+1 n1+n2+2 n1+n2+3];
M( inde, inde ) = Me;

% construction of the load vector
F = zeros( N, 1 );
F( [ 1 n1+n2+1 ] ) = - f0 * dx / 2;
F( 2 : n1+n2 ) = - f0 * dx;
% add the volumic load term for the enriched DOFs !!!
F( n1+n2+1 ) = F( n1+n2+1 ) + F0;

% adding boundary conditions
M = M / dx;
if uL~=0
    BcM = [ zeros(nx-2,1) ; 1; 0 ; 0 ];
    BcF = uL;
    M = [ M BcM ; BcM' 0 ];
    F = [ F ; BcF ];
elseif fa~=0
    F( n1 ) = F( n1 ) -(1-a)^2/2 * dx * fa;
    F( n1+1 ) = F( n1+1 ) -(1-a^2)/2 * dx * fa;
    F( n1+1 : n1+n2 ) = F( n1+2 : n1+n2 ) - fa * dx;
    F( n1+n2+1 ) = F( n1+n2+1 ) - fa * dx / 2;
    F( n1+n2+2 ) = -2/3 * a * (1-a) * fa;
    F( n1+n2+3 ) = -1/3 * a * (1-a) * (1+2*a) * fa;
end

% resolution of the system
X = M \ F;

% construction of the displacement
Uxh = zeros( size( xu ) );
Uxh( 2:n1+1 ) = X( 1:n1 );
Uxh( n1+2 ) = X( n1 ) * (1-aloc) + X( n1+1 ) * aloc;
Uxh( n1+3:end ) = X( n1+1:n1+n2+1 );
Uxh( n1+2 ) = Uxh( n1+2 ) + 2*aloc*( (1-aloc) * X(n1+n2+2) + ...
                                                       aloc * X(n1+n2+3) );

% construction of the flux
Sxh = zeros( size( xs ) );
Sxh( 1:2:2*n1 ) = K1 * ( X( 1:n1 ) - [ 0 ; X(1:n1-1) ] );
Sxh( 2:2:2*n1 ) = K1 * ( X( 1:n1 ) - [ 0 ; X(1:n1-1) ] );
Sxh( 2*n1+[1:2] ) = K1 * ( X( n1+1 ) - X( n1 ) );
Sxh( 2*n1+[3:4] ) = K2 * ( X( n1+1 ) - X( n1 ) );
Sxh( 2*n1+5:2:end ) = K2 * ( X( n1+2:n1+n2+1 ) - X(n1+1:n1+n2) );
Sxh( 2*n1+6:2:end ) = K2 * ( X( n1+2:n1+n2+1 ) - X(n1+1:n1+n2) );
RPhi = [2*K1 2*(1-2*aloc)*K1 -4*aloc*K2                      0
        0 4*aloc*K1   2*aloc*(1-2*aloc)/(1-aloc)*K2 -2*aloc/(1-aloc)*K2]';
Sxh( 2*n1+[1:4] ) = Sxh( 2*n1+[1:4] ) + RPhi * X( n1+n2+[2 3]);
Sxh = Sxh / dx;

%==========================================================================
% COMPUTATION OF THE ERROR BETWEEN THE EXACT AND APPROXIMATE SOLUTIONS
%==========================================================================
eSh = Sh - Se;
eSxh = Sxh - Se;

figure(1); plot( xu, [ Ue Uh Uxh ],'.-')
legend('exact','FEM','XFEM')
figure(2); plot( xs, [ Se Sh Sxh ],'.-')
legend('exact','FEM','XFEM')
figure(3); plot( xs, [ Se eSh eSxh ],'.-')
legend('exact','error FEM','error XFEM');grid on
            