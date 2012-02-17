function [typesOfElements,enrichedNodes] = classifyElements( levelSet, T, tol )
% [typesOfElements,enrichedNodes] = classifyElements( levelSet, T, tol )
%
% returns the type of ech element with the following criterion
%    type = 0    standard element
%    type = +n   enriched element crossed by the n-th level set
%    type = -1   enriched element crossed by more than one level set
%
% INPUT
%   levelSet   nodal level set fields. Each column is a different level set
%   T          conectivity matrix
%   tol        tolerance 
%
% OUTPUT
%   typesOfElements 
%   enrichedNodes
%
%ngeom = 3;
%typesOfElements = zeros( size( T, 1 ), 1 );
%enrichedNodes = [];
typesOfElements = crossedByLevelSet( levelSet, T, tol );

%for ielem = 1:size( T, 1 ) 
%   [bool,nls,idls] = crossedByLevelSet( levelSet(T(ielem,1:ngeom),:), tol );
%   if bool
%      if nls == 1
%         typesOfElements(ielem) = idls;
%      else
%         typesOfElements(ielem) = -1;
%      end
%      enrichedNodes = [enrichedNodes; T(ielem,1:ngeom)'];
%   end
%end
%enrichedNodes = unique( enrichedNodes );
%enrichedNodes = sort( enrichedNodes );
enrichedNodes = T(typesOfElements,:);
%enrichedNodes = T(find(typesOfElements),:);
enrichedNodes = unique( enrichedNodes(:) );