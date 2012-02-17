function [type,Xe,Te,LSe,lin] = reOrderLocalNodes(Xe,Te,LSe,Seg,i1,DDLgamma)
% [type,Xe,Te,LSe,lin] = reOrderLocalNodes(Xe,Te,LSe,Seg,i1,DDLgamma)
%
%
% type = 1 -> level set crosses one node
% type = 2 -> level set does not cross any node
%
ind0 = find( [prod(LSe([2 3])) prod(LSe([1 3])) prod(LSe([1 2]))] > 0 );
if isempty( ind0 )
   %% Case Ls crosses one node
   type = 1;
   % re-ordering to have crossed edge 2-3
   ind0 = find( LSe==0 );
   ind = circshift( [1 2 3]', 1-ind0 );
   Xe = Xe(ind,:);
   LSe = LSe(ind,:);   
   % 
   if Seg(i1,1)==Te(1)
      lin = DDLgamma( i1, : );
   else
      lin = DDLgamma( i1, [2 1] );
   end
   
else
   %% Case Ls does not cross cross any node
   type = 2;
   % re-ordering to have crossed edges 1-2 and 1-3
   ind = circshift( [1 2 3]', 1-ind0 );
   Xe = Xe(ind,:);
   LSe = LSe(ind,:);
   %
   [i2,j2] = find( Seg(i1+(0:1),:)==Te(2) );
   [i3,j3] = find( Seg(i1+(0:1),:)==Te(3) );
   if i2<i3
      lin = DDLgamma( i1, : );
   else
      lin = DDLgamma( i1, [2 1] );
   end
   
end
