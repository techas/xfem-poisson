function [R,Rxi,Reta] = buildCompositeRidge( LSe, N, Nxi, Neta )
% [R,Rxi,Reta] = buildCompositeRidge( LSe, N, Nxi, Neta )
%
% Build the ridge function inside an element crossed by one or more level
% sets. 
%
% INPUT
%    LSe         nodal level set values. matrix. 
%                each column is a different level set.
%                
%    N,Nxi,Neta  shape functions and its derivatives
%
% OUTPUT
%   R            ridge function in gauss points
%   Rxi,Reta     derivatives of the rigde function
%
[R,Rxi,Reta] = buildSimpleRidge( LSe(:,1), [], [], [], [], N, Nxi, Neta );
parentR = R;
parentRxi = Rxi;
parentReta = Reta;
for I = 2:size( LSe, 2 )
   [pR,pRxi,pReta] = buildSimpleRidge( LSe(:,I), LSe(:,I-1), ...
      parentR, parentRxi, parentReta, N, Nxi, Neta );
   R = R + pR;
   Rxi = Rxi + pRxi;
   Reta = Reta + pReta;
   parentR = pR;
   parentRxi = pRxi;
   parentReta = pReta;
end