function [K,f,enrichedNodes] = createMatrix( ...
    X, T, levelSet, pospg, ...
                  pespg, N, Nxi, Neta, caseLoad, tolerance ) 
% [K,f] = createMatrix( X, T, levelSet, pospg, pespg, N, Nxi, Neta ) 
%
% creates the matriz and the source term f 
% 
% INPUT
%   X,T           the mesh
%   levelSet
%   pospg,pespg   position and weigth of gauss points 
%   N,Nxi,Neta    shape functions and its derivetives
%   tolerance     Level set tolerance to detect if the node is on the 
%                 interface
%
% OUTPUT
%   K
%   f
%
global useEnrichment debugData

[nelem,nnode] = size( T ); 
numberOfNodes = size( X, 1 );
%
% Look for the enriched elements and nodes
[type,enrichedNodes] = classifyElements( levelSet, T, tolerance );
standardElements = find( type == 0 );
enrichedElements = find( type > 0 );
if useEnrichment
   ksize = numberOfNodes + length( enrichedNodes );
else
   ksize = numberOfNodes;
end

K = zeros( ksize, ksize );
f = zeros( ksize, 1 );

% fill the matrix for the standard elements
for I = 1:length( standardElements )
   Te = T(standardElements(I),:); 
   Xe = X(Te,:);
   LSe = levelSet(Te);
   material = mean(materialpg( LSe((abs(LSe)==max(abs(LSe)))), 0 ));
%   material = materialpg( LSe(1), tolerance );
   [Ke,fe] = elementMatrixStandard( Xe, nnode, pospg, pespg, N, Nxi, ...
                                 Neta, material, caseLoad ); 
   K(Te,Te) = K(Te,Te) + Ke; 
   f(Te) = f(Te) + fe;
end

% fill the matrix for the enriched elements
for I = 1:length( enrichedElements )
   Te = T(enrichedElements(I),:);
   Xe = X(Te,:); 
   LSe = levelSet(Te);
   % Find indices of enriched dofs
   if useEnrichment
      Re = vectorFind( enrichedNodes, Te );
      Te = [Te, Re+numberOfNodes];
   end
   % 
%    fprintf('\n%i',enrichedElements(I))
   [Ke,fe] = elementMatrixEnriched( Xe, LSe, nnode, pospg, pespg, ...
                        N, Nxi, Neta, caseLoad, tolerance );
   K(Te,Te) = K(Te,Te) + Ke; 
   f(Te) = f(Te) + fe; 
end 


debugData.enrichedElements = enrichedElements;
debugData.enrichedNodes = enrichedNodes;
