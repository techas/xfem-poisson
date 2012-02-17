function [R,Rxi,Reta] = buildSimpleRidge( levelSet, parentLevelSet, ...
   parentR, parentRxi, parentReta, N, Nxi, Neta )
% [R,Rxi,Reta] = ridge( levelSet, parentLevelSet, parentR, ...
%   parentRxi, parentReta, N, Nxi, Neta )
%
% Build the ridge function and its derivatives for one level set when 
% multi level set are used. For linear triangular elements. 
%
% The ridge used is defined in
% G.Legrain, N. Moes, A. Huerta. Stability of incompressible 
%    formulation enriched with X-FEM (preprint Elsevier Science June 2005)
%
% INPUT
%    levelSet        nodal level set values
%    parentLevelSet  nodal level set values of the level set
%    parentR         ridge function of the parent level set
%    parentRxi       and derivatives
%    parentReta
%    N,Nxi,Neta      shape functions and its derivatives
%
% OUTPUT
%   R           ridge function in gauss points
%   Rxi         derivatives of the rigde function
%   Reta
%
s = double( N*levelSet > 0 );
s(s==0) = -1;
R = N*abs( levelSet ) - abs( N*levelSet );
Rxi = Nxi*abs( levelSet ) - s.*(Nxi*levelSet) ;
Reta = Neta*abs( levelSet ) - s.*(Neta*levelSet) ;
if ~isempty( parentLevelSet )
   A = N*parentLevelSet;
   B = N*abs( parentLevelSet );
   %
   C = parentR ./ B;
   C(A <= 0) = 1;
   Cxi = parentRxi.*B - parentR.*(Nxi*parentLevelSet);
   Cxi(A <= 0) = 0;
   Ceta = parentReta.*B - parentR.*(Neta*parentLevelSet);
   Ceta(A <= 0) = 0;
   %
   R = R.*C;
   Rxi = Rxi.*C + R.*Cxi;
   Reta = Reta.*C + R.*Ceta;
end