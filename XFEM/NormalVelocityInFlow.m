function Vn = NormalVelocityInFlow( X, T, pospg, N, ...
                                        Nxi, Neta, levelSet, h, hEnr, opts)
% NORMALVELOCITYINFLOW to compute the right-hand side term used to update the 
% level set function
%    d_t(Phi) = -H(q,dPhi)
%
% syntax: Vn = NormalVelocityInFlow( X, T, pospg, N, ...
%                                         Nxi, Neta, levelSet, h, hE, opts)
%
% X:            nodal coordinates
% T:            Connectivity table
% pospg, pespg: positions and weights of Gauss points
% N, Nxi, Neta: FE interpolation functions and derivatives
% levelSet:     levelSet function given as a value at each node (positive
%               or negative indicating the side)
% h, hE:        value of the unknowns at the normal and enriched DOFs
% opts:         structured array of options
%     opts.EvolutionLaw : to choose between different evolution laws
%           1: H = |q|.|grad(Phi)| 
%           2: H = (q,n).|grad(Phi)| (default)

% R. Cottereau 03/2008 (using parts by S. Zlotnik and P. Diez)

% DEFAULT VALUE OF THE OPTION
if ~isfield( opts, 'EvolutionLaw' )
    opts.EvolutionLaw = 2;
end

% number of elements, of nodes in each element, of total nodes
Nn = size(X,1);

% Get the elements in a layer around the interface
[type,enrichedNodes] = classifyElements( levelSet, T, 0 );
enrichedElements = find( type > 0 );
BoundaryLayer = inContactWithElements( X, T, enrichedElements );
Tb = T( BoundaryLayer, : );

% Compute the fluxes in the gauss point of these elements
[ Xpg, Fpg ] = Flujos( h, hEnr, X, Tb, pospg, N, Nxi, Neta, levelSet, opts);
dPhi = GradientLevelSet( levelSet, X, Tb, pospg, N, Nxi, Neta, opts );
Fg = VariableLaw( Fpg, dPhi, opts.EvolutionLaw );

% Initialization
Vn = zeros( Nn, 1 );

% Compute a cuadratic approximation over the elements touching each
% enriched node
for i1 = 1:length( enrichedNodes )
    inode = enrichedNodes(i1);
    [ EltsInContact, j1 ] = find( Tb == inode );
    Te = Tb( EltsInContact, : );
    LocXpg = InStar( Xpg, X, Te );
    LocFpg = Fg( LocXpg, : );
    LocXpg = Xpg( LocXpg, : );
    MatCoord = [ ones(size(LocXpg,1),1) LocXpg ];
    Vn( inode ) = [1 X(inode,:)] * ( (MatCoord'*MatCoord)\(MatCoord'*LocFpg) );
end
Vn = max( Vn, 0 );
                       
%==========================================================================
function Fg = VariableLaw( Fluxg, dPhi, Law )
% VARIABLELAW to choose the 
% Law:    1  for d_t(Phi) = - |q|.|d_x(Phi)| = -H
%         2  for d_t(Phi) = - (q,n).|d_x(Phi)| = -H
switch Law
    
    % H = |q|.|d_x(Phi)|
    case 1
        Fg = sqrt( sum( Fluxg.*Fluxg, 2 ) );
        
    % H = (q,n).|d_x(Phi)|
    case 2
        Fg = abs( sum( Fluxg .* dPhi, 2 ) );
                
    otherwise
        error('unknown evolution law');
end
