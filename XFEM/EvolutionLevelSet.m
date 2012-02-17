function [ levelSet, dt ] = EvolutionLevelSet( X, T, levelSet, Vn, ...
                                                 N, Nxi, Neta, Lcar, opts )
% NORMALVELOCITYINIMPACT computes the advance of the level set front when
% impacted by particles of water
%
% syntax: levelSet = NormalVelocityInImpact( X, T, levelSet, Xi, Vi, ...
%                                         Ei, Ni, N, Nxi, Neta, dt, opts );
%
% X,T:        nodal coordinates and connectivity matrix 
% levelSet:   current nodal values of the level set function
% Xi:         position of the impacts of the particles with the level set front
% Vi:         velocities of the particles at impact
% Ei:         element impacted by the particles
% Ni:         normal vector to the impact surfaces
% li:         length of the impacted segments
% N,Nxi,Neta: interpolation function and derivatives on isoparametric
%             element
% opts:       options
%
% levelSet:   new nodal values of the level set function

% R. Cottereau 04/2008

% selection of nodes and elements to be modified
[ type, enrichedNodes ] = classifyElements( levelSet, T, 0 );
enrichedElements = find( type > 0 );
BoundaryElts = inContactWithElements( X, T, enrichedElements );
Tb = T( BoundaryElts, : );

% computation of norm of the gradient of the level set
normGradPhi = GradientLevelSetNodal( X, Tb, levelSet, N, Nxi, Neta, opts );

% choice of the time step
% ind = enrichedNodes( ( Vn( enrichedNodes ) ~= 0 ) );
% dt = 1 / 3 * min( Lcar( ind ) ./  Vn( ind ) );
% 
% dt = min( dt, 13 );
% NEW
elemSize = calcMinEdgeLength(X,T);
dt = calcDT(Vn,elemSize);
% keyboard


% Solve evolution equation
if dt ~= 0
    levelSet = levelSet - dt * ( Vn .* normGradPhi );
end