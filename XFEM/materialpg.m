function mat = materialpg( lspg, tol )
% mat = materialpg( lspg, tol )
%
% Calculate the material phase in each gauss point
%
% INPUT
%   lspg  value of all level sets in the gaus points. each column has a
%         different level set
%   tol   tolerance
%
% OUTPUT
%   mat   material phase of each gauss point
%
if nargin == 1
   tol = 0;
end
[ngaus,nls] = size( lspg );
mat = zeros( ngaus, 1 );

for igaus = 1:ngaus
   for ils = 1:nls
      if lspg(igaus,ils) >= - tol
         mat(igaus) = ils;
         break
      end
   end
end
% if no material has been assigned, the point is in the negative part of
% the last level set and corresponds to the last material
mat(mat==0) = nls + 1;
