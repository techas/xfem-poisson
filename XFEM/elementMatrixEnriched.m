function [Ke,fe] = elementMatrixEnriched( Xe, levelSet, numberOfNodes, ...
   pospg, pespg, N, Nxi, Neta, caseLoad, tol )
% [Ke,fe] = elementMatrixEnriched( Xe, numberOfNodes, pospg, pespg, ...
%   N, Nxi, Neta )
%
% builds the K and f matrices for a triangular element crossed by the level
% set. Integrates by splitting the element in single material parts.
%
% INPUT
%   Xe             nodal coords
%   levelSet       nodal level set 
%   numberOfNodes  number of element nodes
%   pospg,pespg    position and weigth of gauss points 
%   N,Nxi,Neta     shape functions and its derivetives
%
% OUTPUT
%   K
%   f
%
global useEnrichment debug
tol = 0;

% Look how the element is crossed by the level set
numberOfNegativeNodes = length( find( levelSet < -tol ) );
numberOfZeroNodes = length( find( abs( levelSet ) < tol ) );
code = 10*numberOfZeroNodes + numberOfNegativeNodes;

switch code
   case 11
      % The element is crossed by the level set, 
      % and split the element in 2 triangles
      %
      %             level set 
      %                 .  
      %                 .
      %                 o nodosolo
      %                /.\
      %               / . \
      %              /  .  \
      %             /ma2.ma1\
      %  nodopar(2)o----+----o nodopar(1)
      %                 . newNode
      %                 .

      % Node shared by the 2 triangles
      nodosolo = find( abs( levelSet ) < tol );
      nodopar = find( abs( levelSet ) > tol );

      % Compute the intersection point between the element edge and the level set
      newNode = intersection( Xe(nodopar,:), levelSet(nodopar) );
      Xp = [Xe; newNode];

      % Local numbering of each triangle
      switch(nodosolo)
         case 1
            index = [1 2 4; 1 4 3];
         case 2
            index = [1 2 4; 4 2 3];
         case 3
            index = [1 4 3; 4 2 3];
      end
      % Material 
      if levelSet(nodopar(1)) > 0
         material = [1 2];
      else
         material = [2 1];
      end

   case {1,2}
      % The element is split in 3 triangles
      %
      %                  nodosolo   
      %                 o 
      %                / \
      %     newNode1  /   \  newNode2               mat1
      % .............x.....x.......... level set ...........
      %             /       \                      mat2 mat3
      %            o---------o
      %   nodopar(1)          nodopar(2)

      % Find nodes with negative level set
      negNodes = find( levelSet < -tol );
      switch length( negNodes )
         case 1
            nodosolo = negNodes;
            nodopar = find( levelSet >= 0 );
            material = [1 1 2];
         case 2
            nodosolo = find( levelSet >= 0 );
            nodopar = negNodes;
            material = [2 2 1];
         otherwise
            error( 'elementMatrixEnriched' )
      end   
      % Compute the intersection nodes
      X1 = [Xe(nodosolo,:); Xe(nodopar(1),:)];
      X2 = [Xe(nodosolo,:); Xe(nodopar(2),:)];
      LS1 = [levelSet(nodosolo) levelSet(nodopar(1))];
      LS2 = [levelSet(nodosolo) levelSet(nodopar(2))];
      newNode1 = intersection( X1, LS1 );
      newNode2 = intersection( X2, LS2 );
      Xp = [Xe; newNode1; newNode2];
      switch( nodosolo )
         case 1
            index = [5 4 3; 4 2 3; 1 4 5 ];
         case 2
            index = [1 4 3; 4 5 3; 4 2 5 ];
         case 3
            index = [1 2 4; 4 2 5; 4 5 3 ];
      end
      cortes = [ nodopar(1) nodosolo ; 
                 nodopar(2) nodosolo ];
   
   case {0, 10, 20}
      % The element is NOT crossed and the level set is positive
      levelSet
      Xe
      error( 'Error: this element does not need enrichment' );
      
   case {3, 12, 21}
      % The element is NOT crossed and the level set is negative
      error( 'Error: this element does not need enrichment' );
      
   case 30
      disp( Xe )
      disp( levelSet )
      error( 'Error: one element has 3 nodes with zero level set' );
     
   otherwise
      disp( Xe )
      disp( levelSet )
      error( 'Error: unknow code %i', code );
end

% Integrate each sub element
if useEnrichment
   ksize = 2*numberOfNodes;
else
   ksize = numberOfNodes;
end
Ke = zeros( ksize, ksize ); 
fe = zeros( ksize, 1 ); 

area = 0;
for I = 1:size( index, 1 )
   [Kep,fep] = elementMatrixEnrichedPart( Xe, Xp(index(I,:),:), levelSet, ...
        numberOfNodes, pospg, pespg, N, Nxi, Neta, material(I), caseLoad );
   Ke = Ke + Kep;
   fe = fe + fep;
   
   %%
   if material(I) == 2
      area = area + miArea(Xp(index(I,:),:));
   end
   %%
end   
debug = [debug; area/(miArea(Xe)-area)];

%
% vaules of code (*)
%
% numZero   numNeg    code    comment
% -------   -------    ---    --------------------------------------------
%    0         0        0     Element in the positive side of the level set
%                             Use properties 1
%    0         3        3     Element in the negative side of the level set
%                             Ese properties 2
%    0         1        1     Element crossed by the level set. 
%                             The levelSet does not touch any node.
%                             Split the element in 3 triangles
%                             The triangular side is negative
%    0         2        2     Element crossed by the level set. 
%                             The levelSet does not touch any node.
%                             Split the element in 3 triangles
%                             The triangular side is positive
%    1         0       10     Element in the positive side of the level set
%                             The level set touchs one node
%                             Use properties 1
%    1         2       12     Element in the negative side of the level set
%                             The level set touchs one node
%                             Use properties 2
%    1         1       11     Element crossed by the level set
%                             Split the element in 2 triangles 
%                             (the level set touchs a node and crosses 
%                             one edge)
%    2         1       21     The level set goes over one edge
%                             The element is in the negative side
%                             Use properties 2
%    2         0       20     The level set goes over one edge
%                             The element is in the positive side
%                             Use properties 1
%    3         0       --     Error: the level set is zero in al the nodes
%                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
