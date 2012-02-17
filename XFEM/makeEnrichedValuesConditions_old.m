function [Abc,bbc] = makeEnrichedValuesConditions( X, T, ...
                                                 Nxi, Neta, levelSet, opts)
% only valid for triangular elements !!!
if ~opts.isTriangle
    error( ['this function (makeEnrichedValuesConditions) was ' ...
           'implemented only for triangles !!' ])
end

% constants
global cond
Nn = size( X, 1 );
B0 = [ Nxi(1,:); Neta(1,:) ];
lambda = [0 1];
Nl = length(lambda);

% enriched nodes
[ type, enrichedNodes ] = classifyElements( levelSet, T, opts.tolerance );
enrichedElements = find( type > 0 );
Ne = length( enrichedElements );
Nne = length( enrichedNodes );

% initializations
Ru = zeros( Ne, Nn );
Ra = zeros( Ne*Nl, Nne );

% Loop on enriched elements
for i1 = 1:Ne
    
    % enriched element
    Te = T( enrichedElements(i1), : );
    LSe = levelSet( Te );

    % intersection points and re-ordering
    ind0 = find( [prod(LSe([2 3])) prod(LSe([1 3])) prod(LSe([1 2]))] > 0 );
    ind = [ ind0 setdiff( 1:3, ind0 ) ];
    Te = Te( ind );
    Xe = X( Te, : );
    LSe = LSe( ind );
    aLSe = abs(LSe);
    NP1 = [ LSe(2); -LSe(1); 0 ] / (LSe(2)-LSe(1));
    NP2 = [ LSe(3); 0; -LSe(1) ] / (LSe(3)-LSe(1));

    % computation of r_u, independent of the choice of P1 and P2
    Be = ( B0 * Xe ) \ B0;
    r_u  = Be' * Be * LSe;
        
    % computation of r_a, for the enriched nodes
    phi = aLSe - (cond(1)+cond(2))/(cond(1)-cond(2))*LSe;
    r_a = (NP1'*aLSe*r_u+(phi'*r_u)*NP1)*(1-lambda) ...
        + (NP2'*aLSe*r_u+(phi'*r_u)*NP2)*lambda;
        
    % introduce the values in the matrix
    Ru( i1, Te ) = r_u';
    Ra( i1:Ne:end, vectorFind( enrichedNodes, Te ) ) = r_a';
    
end

% global construction
Abc = [ repmat(Ru,[Nl 1]) Ra ];
bbc = zeros( Nl*Ne, 1);

