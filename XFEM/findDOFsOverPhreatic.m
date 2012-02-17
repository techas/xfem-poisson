function [ ind, indE ] = findDOFsOverPhreatic( X, T, N, levelSet, ...
                                                     phreaticLevel, opts );
% FINDDOFSOVERPHREATIC to find the indices of the DOFs that only touch
% elements whos Gauss point are over the water table
%
% syntax: [ ind, indE ] = findDOFsOverPhreatic( X, T, N, ph );
%
%  X:    nodal coordinates [Nn*d matrix]
%  T:    connectivity matrix [Ne*e matrix]
%  N:    FE interpolation function values at the Gauss points [Ng*e matrix]
%  ph:   phreatic level [scalar]
%
%  ind:  indices of the DOFs that touch elements that have elements below
%        the water table and hence should be considered for the
%        computations
%  indE: same as ind, but for enriched DOFs. the numbering of indE is
%        global in the sense that it starts at the end of the numbering of
%        ind.

% R. Cottereau 05/2008

% constants
global useEnrichment
numberOfNodes = size( X, 1 );

% get enriched nodes
[ type, enrichedNodes ] = classifyElements( levelSet, T, opts.tolerance );

% logical indicating nodes that are below the phreatic level
below = ( X(:,2) < phreaticLevel );

% find elements that are completely or partially below the phreatic level
ElementsBelow = below( T(:,1) ) & below( T(:,2) ) & below( T(:,3) );
ElementsPartiallyBelow = ~ElementsBelow & ...
                   ( below( T(:,1) ) | below( T(:,2) ) | below( T(:,3) ) );
               
% transform to indices
ElementsBelow = find( ElementsBelow );
ElementsPartiallyBelow = find( ElementsPartiallyBelow );

% add elements that have at least one Gauss point below the water table
for i1 = 1:length( ElementsPartiallyBelow )
    Te = T( ElementsPartiallyBelow( i1 ), : );
    Ye = X( Te, 2 );
    Yg = Isopar( Ye, N );
    if any( Yg < phreaticLevel )
        ElementsBelow = [ ElementsBelow ; ElementsPartiallyBelow( i1 ) ];
    end
end

% find nodes that should be considered in the FE computation
ind = T( ElementsBelow, : );
ind = unique( ind(:) );

% find enriched nodes that should be considered in the FE computation
if useEnrichment
    indE = intersect( ind, enrichedNodes );
    indE = vectorFind( enrichedNodes, indE ) + numberOfNodes;
    indE = indE(:);
else
    indE = [];
end
