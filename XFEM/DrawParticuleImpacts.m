function DrawParticuleImpacts( h, X, T, levelSet, x0, y0, vx0, vy0, ...
                                                         xi, yi, vxi, vyi);
% DRAWPARTICULESIMPACTS to draw the particules and their velocities at both
% initial and impacting positions
%
%  syntax: DrawParticuleImpacts( h, X, T, leveSet, x0, y0, vx0, vy0, ...
%                                                     [xi, yi, vxi, vyi]);
%
%  h:       handle for the figure [integer]
%  X,T:     nodal coordinates and connectivity matrix
%  levelSet: nodal values of the level set function
%  x0,y0:   initial positions of the particules [n*1 vectors]
%  vx0,vy0: velocities of the particules in their initial position 
%           [n*1 vectors]
%  xi,yi:   position of particules at impact [n*1 vectors]
%  vxi,vyi: velocities of the particules at impact [n*1 vectors]
%
%  it is also possible to draw only the initial positions and velocities by
%  taking out the last four arguments

% R. Cottereau 04/2008

% figure
figure( h );

% constants
n = 100;

% limiting the input
x0 = x0(1:n,:);
y0 = y0(1:n,:);
vx0 = vx0(1:n,:);
vy0 = vy0(1:n,:);
if nargin==12
    xi = xi(1:n,:);
    yi = yi(1:n,:);
    vxi = vxi(1:n,:);
    vyi = vyi(1:n,:);
end

% mesh
trimesh( T, X(:,1), X(:,2), 'Color', 'k' );

% particules in initial positions
hold on;
scatter( x0, y0, 50, 'b' );

% velocities in initial positions
DibujaVec( [x0 y0], X, T, [vx0 vy0], levelSet, 0, 'b-', .7 )

if nargin==12
    % particules in impacting positions
    hold on;
    scatter( xi, yi, 50, 'r' );

    % velocities in impacting positions
    DibujaVec( [xi yi], X, T, [vxi vyi], levelSet, 0, 'r-', .7 )
end
