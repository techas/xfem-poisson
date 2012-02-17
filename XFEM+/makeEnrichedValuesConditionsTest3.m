function R = makeEnrichedValuesConditionsTest3( X, T, ...
   Nxi, Neta, levelSet, opts)
% R = makeEnrichedValuesConditions( X, T, ...
%                                                 Nxi, Neta, levelSet, opts)
% only valid for triangular elements !!!
%
global cond useEnrichment debugData

if ~opts.elementType==2
    error( ['this function (makeEnrichedValuesConditions) was ' ...
           'implemented only for triangles !!' ])
end

theNormals = [];
Nn = size( X, 1 );
B0 = [ Nxi(1,:); Neta(1,:) ];

% look for enriched elements and nodes
[ type, enrichedNodes ] = classifyElements( levelSet, T, opts.tolerance );
enrichedElements = find( type > 0 );

if 0
   %% Test: 
   % boundaary elements with a neumann-type BC are NOT enforced to 
   % have zero flux jump across the interface
   
   % Brutal code next (only works for square domains and triangular elements)
   xmin = min(X(:,1));
   xmax = max(X(:,1));
   sideNodes = find(X(:,1) == xmin | X(:,1) == xmax);
   [~,it1,~] = intersect(T(:,1), sideNodes);
   [~,it2,~] = intersect(T(:,2), sideNodes);
   [~,it3,~] = intersect(T(:,3), sideNodes);
   sideElements = unique([it1; it2; it3]);
   n = length(sideElements);
   enrichedElements = setdiff(enrichedElements, sideElements);
   disp( n-length(enrichedElements) );
end



if 1
   [ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet, opts.tolerance );
   [ polis, Ei ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, levelSet );
end

% initializations
Nne = length( enrichedNodes );
if useEnrichment
   ksize = Nn + Nne;
else
   ksize = Nn;
end
Nl = 2*length(enrichedElements);
R = zeros( Nl, ksize );

% Loop on enriched elements
for i1 = 1:length(enrichedElements)
   
   % enriched element
   Te = T( enrichedElements(i1), : );
   LSe = levelSet( Te );
   
   % intersection points
   ind0 = find( [prod(LSe([2 3])) prod(LSe([1 3])) prod(LSe([1 2]))] > 0 );
   
   if isempty( ind0 )
      %% CASE LEVEL SET CROSSES ONE NODE
      % re-ordering
      ind0 = find( LSe==0 );
      ind = circshift( [1 2 3]', 1-ind0 );
      Te = Te( ind );
      % dofs location
      lin = (i1-1)*2 + (1:2);
      enrTe = Nn + vectorFind( enrichedNodes, Te );
      % ordered element data
      Xe = X( Te, : );
      LSe = levelSet( Te );
      aLSe = abs(LSe);
      % p1 , p2
      P1 = Xe(1,:);
      P2 = intersection( Xe([2 3],:), LSe([2 3]) );
      LInt = norm( P2 - P1 );
      % shape functions @ p1, p2
%       NP1 = [ 1 0 0 ]'
%       NP2 = [ 0; LSe(3); -LSe(2) ] / (LSe(3)-LSe(2))
   else
      %% CASE LEVEL SET DOES NOT CROSS ANY NODE
      % re-ordering
      ind = circshift( [1 2 3]', 1-ind0 );
      Te = Te( ind );
      % dofs location
      lin = (i1-1)*2 + (1:2);
      enrTe = Nn + vectorFind( enrichedNodes, Te );
      % ordered element data
      Xe = X( Te, : );
      LSe = levelSet( Te );
      aLSe = abs(LSe);
      % p1, p2
      P1 = intersection( Xe([1 2],:), LSe([1 2]) );
      P2 = intersection( Xe([1 3],:), LSe([1 3]) );
      LInt = norm( P2 - P1 );
%       NP1 = [ LSe(2); -LSe(1); 0 ] / (LSe(2)-LSe(1))
%       NP2 = [ LSe(3); 0; -LSe(1) ] / (LSe(3)-LSe(1))      
   end
   
   % shape functions @ p1, p2
   [xi,eta] = invMap( Xe, [P1(1) P2(1)], [P1(2) P2(2)] );
   triangle = 2;
   nNodes = 3;
   outN = shapeFunction(triangle,nNodes,[xi',eta']);
   NP1 = outN(1,:)';
   NP2 = outN(2,:)';

   % grad x-y
   B1 = B0(:,ind);
   Be = ( B1 * Xe ) \ B1;
   
   % analytic integration of r_u (standard dofs)
   r_u = Be' * Be * LSe;
   r_u = (cond(2)-cond(1)) * r_u * LInt / sqrt(LSe'*r_u);
   
   R(lin,Te) = R(lin,Te) + ( [1/2 1/2]' * r_u' );
   
   % only for debug
   theNormals(i1,:) = Be*LSe/norm(Be*LSe);
   
   % computation of R, enriched part
   if useEnrichment
      phi = aLSe - ((cond(1)+cond(2))/(cond(1)-cond(2))*LSe);
      r_a1 = ( (NP1' * aLSe) * r_u ) + ( (phi' * r_u) * NP1 );
      r_a2 = ( (NP2' * aLSe) * r_u ) + ( (phi' * r_u) * NP2 );
      R(lin,enrTe) = R(lin,enrTe) + ...
         (([1/3 1/6]' * r_a1') + ([1/6 1/3]' * r_a2'));
   end
end



%% debug and plot stuff
if 1
   debugData.enrichedElements = enrichedElements;
   debugData.enrichedNodes = enrichedNodes;
   debugData.polis = polis;
   debugData.Ei = Ei;
   debugData.theNormals = theNormals;
end
