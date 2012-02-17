
R = 1;
NR = 320;
theta = linspace(0,2*pi,NR+1);
theta = theta(1:end-1)';
nodes = R*[ cos(theta) sin(theta) ];
hdata.hmax = 2*pi*R/NR;
[ X, T ] = mesh2d( nodes, [], hdata );

% save MeshDisc5.mat X T
