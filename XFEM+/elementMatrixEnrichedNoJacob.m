function Ke = elementMatrixEnrichedNoJacob( Xe, levelSet, pespg, Nx, Ny )


% Look how the element is crossed by the level set
numberOfNegativeNodes = length( find( levelSet < 0 ) );
numberOfZeroNodes = length( find( abs( levelSet ) < 0 ) );
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
      nodosolo = find( abs( levelSet ) < 0 );
      nodopar = find( abs( levelSet ) > 0 );

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
      negNodes = find( levelSet < 0 );
      switch length( negNodes )
         case 1
            nodosolo = negNodes;
            nodopar = find( levelSet >= 0 );
          case 2
            nodosolo = find( levelSet >= 0 );
            nodopar = negNodes;
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
   
   case {0, 10, 20}
      % The element is NOT crossed and the level set is positive
      error( 'Error: this element does not needs enrichment' );
      
   case {3, 12, 21}
      % The element is NOT crossed and the level set is negative
      error( 'Error: this element does not needs enrichment' );
      
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
Ke = 0; 
for I = 1:size( index, 1 )
   Kep = elementMatrixStandardNoJacob( Xp(index(I,:),:), pespg, Nx, Ny );
   Ke = Ke + Kep;
end   

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
