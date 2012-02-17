function R = makeEnrichedValuesConditionsLeastSquares( X, T, N, ...
                                                 Nxi, Neta, levelSet, opts)
%function [Abc,bbc] = makeEnrichedValuesConditionsLeastSquares( X, T, N, ...
%                                                 Nxi, Neta, levelSet, opts)
% only valid for triangular elements !!!
if ~opts.elementType==2
    error( ['this function (makeEnrichedValuesConditions) was ' ...
           'implemented only for triangles !!' ])
end

% constants
global cond useEnrichment debugData

theNorms = {};
Nn = size( X, 1 );
B0 = [ Nxi(1,:); Neta(1,:) ];

%


% connectivity between nodes of the interface and global numbering of DOFs
[ type, enrichedNodes ] = classifyElements( levelSet, T, opts.tolerance );
enrichedElements = find( type > 0 );
Nne = length( enrichedNodes );
if useEnrichment
   ksize = Nn + Nne;
else
   ksize = Nn;
end
[ Seg, SegsBnd ] = CrossedSegments( T, enrichedElements, levelSet, opts.tolerance );
[ polis, Ei, Segi ] = MakePoligonalFromSegments( X, T, SegsBnd, Seg, levelSet );

Nc = length(polis);
Ni = zeros(Nc,1);
DDLgammai = cell(Nc,1);
for i1=1:Nc
    ni = size( Segi{i1}, 1 );
    DDLgammai{i1} = sum(Ni)+[1:(ni-1);2:ni]';
    if all( Segi{i1}(1,:)==Segi{i1}(end,:) )
        Ni(i1) = ni-1;
        DDLgammai{i1}( end, 2 ) = DDLgammai{i1}( 1, 1 );
    else
        Ni(i1) = ni;
    end
end
Nl = sum(Ni);

% initializations
R = zeros( Nl, ksize );
% conductivity jump
deltaNu = cond(2)-cond(1);



% Loop on connectivity of matrix
for i0 = 1:Nc
   Ne = length(Ei{i0});
   E = Ei{i0};
   Seg = Segi{i0};
   DDLgamma = DDLgammai{i0};
   
   % Loop on enriched elements
   for i1 = 1:Ne
      
      % enriched element
      Te = T( E(i1), : );
      LSe = levelSet( Te );
      
      % intersection points
      ind0 = find( [prod(LSe([2 3])) prod(LSe([1 3])) prod(LSe([1 2]))] > 0 );
      
      % CASE LEVEL SET CROSSES ONE NODE
      if isempty( ind0 )
         ind0 = find( LSe==0 );
         ind = circshift( [1 2 3]', 1-ind0 );
         Te = Te( ind );
         if Seg(i1,1)==Te(1)
            lin = DDLgamma( i1, : );
         else
            lin = DDLgamma( i1, [2 1] );
         end
         enrTe = Nn + vectorFind( enrichedNodes, Te );
         Xe = X( Te, : );
         LSe = levelSet( Te );
         aLSe = abs(LSe);
         NP1 = [ 1 0 0 ]';
         NP2 = [ 0; LSe(3); -LSe(2) ] / (LSe(3)-LSe(2));
         
         % computation of R, constant part
         B1 = B0(:,ind);
         Be = ( B1 * Xe ) \ B1;
         r_u  = Be' * Be * LSe;         
         
         theNorms{i0}(i1,:) = Be*LSe/norm(Be*LSe);  
            
         P1 = Xe(1,:);
         P2 = intersection( Xe([2 3],:), LSe([2 3]) );
         LInt = norm( P2 - P1 );
         r_u = (cond(2)-cond(1)) * r_u * LInt / sqrt(LSe'*r_u);
         
            
         R( lin, Te ) = R( lin, Te ) + ( [1/2 1/2]' * r_u' );
         
         % CASE LEVEL SET DOES NOT CROSS ANY NODE
      else
         % re-ordering
         ind = circshift( [1 2 3]', 1-ind0 );
         Te = Te( ind );
         [i2,j2] = find( Seg(i1+(0:1),:)==Te(2) );
         [i3,j3] = find( Seg(i1+(0:1),:)==Te(3) );
         if i2<i3
            lin = DDLgamma( i1, : );
         else
            lin = DDLgamma( i1, [2 1] );
         end
         % dof numbering
         enrTe = Nn + vectorFind( enrichedNodes, Te );
         % Element data
         Xe = X( Te, : );
         LSe = levelSet( Te );
         aLSe = abs(LSe);
         P1 = intersection( Xe([1 2],:), LSe([1 2]) );
         P2 = intersection( Xe([1 3],:), LSe([1 3]) );
         % needed???
         LInt = norm( P2 - P1 );
         % normalized level set
         phi = aLSe - ((cond(1)+cond(2))/(cond(1)-cond(2))*LSe);
         % mid point ???
         Pm = (P1+P2)/2;
         % shape functions @ P1, P2 
         NP1 = [ LSe(2); -LSe(1); 0 ] / (LSe(2)-LSe(1));
         NP2 = [ LSe(3); 0; -LSe(1) ] / (LSe(3)-LSe(1));
         % shape functions @ integration points
         NPM = NLi*[NP1 NP2]';
         % ridge @ P1, P2
         RP1 = NP1(1)*aLSe(1)  +  NP1(2)*aLSe(2);
         RP2 = NP2(1)*aLSe(1)  +  NP2(3)*aLSe(3);
         % ridge @ integration points
         RPM = NLi*[RP1;RP2];         
         % grad N in global coords
         B1 = B0(:,ind);
         Be = ( B1 * Xe ) \ B1;
         
         % the normal
         normal = Be*LSe/norm(Be*LSe);
         
         % [grad_xy N1 * n, grad_xy N2 * n, grad_xy N3 * n]
         gradNn = normal' * Be;

         % test function v (matrix)
         %   each row corresponds to a int. point
         %   each col to a function N_i
         vu = deltaNu * repmat(gradNn',1,2);
         va = deltaNu * repmat(gradNn',1,2)  +  deltaNu * (gradNn' * (RPM + NPM*phi)' );
         
         keyboard
         
         R( lin, Te ) = R( lin, Te ) + ( [1/2 1/2]' * r_u' );
         
         % debug data
         theNorms{i0}(i1,:) = normal;
      end
      
      % computation of R, enriched part
      if useEnrichment
         phi = aLSe - ((cond(1)+cond(2))/(cond(1)-cond(2))*LSe);
         r_a1 = ( (NP1' * aLSe) * r_u ) + ( (phi' * r_u) * NP1 );
         r_a2 = ( (NP2' * aLSe) * r_u ) + ( (phi' * r_u) * NP2 );
         R( lin, enrTe ) = R( lin, enrTe ) + ...
            (([1/3 1/6]' * r_a1') + ([1/6 1/3]' * r_a2'));
      end
   end
end



debugData.enrichedElements = enrichedElements;
debugData.enrichedNodes = enrichedNodes;
debugData.polis = polis;
debugData.Ei = Ei;
debugData.theNorms = theNorms;
