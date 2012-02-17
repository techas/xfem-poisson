function Vn = NormalVelocityInImpact( X, T, levelSet, Vi, ...
                                              Ei, Ni, li, time, Particles )
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

% constants
Nn = size(X,1);

% initializations
ind = ( Ei ~= 0 );
Ei = Ei( ind );
li = li( ind );
Vn = zeros( Nn, 2 );
nx = Ni( ind, 1 );
ny = Ni( ind, 2 );
vx = Vi( ind, 1 );
vy = Vi( ind, 2 );

% selection of nodes and elements to be modified
[ type, enrichedNodes ] = classifyElements( levelSet, T, 0 );
enrichedElements = find( type > 0 );
[ BoundaryElts, BoundaryNodes ] = inContactWithElements( X, T, ...
                                                       enrichedElements );

% normal velocity of the impacting particles (normalized by surface)
NormalVeloc = ( nx .* vx + ny .* vy ) ./ li;

% pass the normal velocity to nodes
NbEltsCont = zeros( Nn, 2 );
[ listE, ind ] = unique(Ei);
NbE = length(listE);
nE = repmat( permute([ nx(ind) ny(ind) ], [3 2 1] ), [size(T,2) 1 1] );
for i1 = 1:NbE;
    ind = ( Ei == listE(i1) );
    sVn = sum( NormalVeloc( ind ) );
    Nods = T( listE(i1), : );
    NbEltsCont( Nods, : ) = NbEltsCont( Nods, : ) + 1;
    n1 = NbEltsCont( Nods, : );
    Vn( Nods, : ) = Vn( Nods, : ) .* (n1-1)./n1 + sVn * nE(:,:,i1) ./ n1;
end

% transfer normal velocity to close-by elements and nodes
BoundaryNodes = setdiff( BoundaryNodes, enrichedNodes );
for i1 = 1:length(BoundaryNodes)
    Nod = BoundaryNodes( i1 );
    [ touchingElts, j1 ] = find( T == Nod );
    touchingElts = T( touchingElts, : );
    contactNodes = intersect( touchingElts(:), enrichedNodes );
    Vn( Nod, : ) = mean( Vn( contactNodes, : ), 1 );
end
Vn = sqrt( Vn(:,1).^2 + Vn(:,2).^2 );

% normalize (see paper for formula)
[ n, v0, S, aS, bS, Qe ] = defineLeak( time.initial, time.t, Particles );
Vn = Vn * Particles.FrontAdvance * Qe / Particles.ParticlesPerTimeStep;
