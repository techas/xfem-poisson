function [R,Rxi,Reta] = buildRidge( levelSet, N, Nxi, Neta, force )
% [R,Rxi,Reta] = buildRidge( levelSet, N, Nxi, Neta )
%
% Build the ridge function and its derivatives. 
% For linear triangular elements. 
%
% The ridge used is defined in
% G.Legrain, N. Moes, A. Huerta. Stability of incompressible 
%    formulation enriched with X-FEM (preprint Elsevier Science June 2005)
%
% INPUT
%    levelSet        nodal level set values
%    N,Nxi,Neta      shape functions and its derivatives
%
% OUTPUT
%   R           ridge function in gauss points
%   Rxi         derivatives of the rigde function
%   Reta
%
%s = double( N*levelSet > 0 );
%s(s==0) = -1;

% modification to force the computation to use one side to compute fluxes
if ( nargin>4 && ~isnan(force))
    s = force;
else
    s = sign( N*levelSet );
end
alevelSet = abs( levelSet );
%NalevelSet = N * abs( levelSet );
%aNlevelSet = abs( N * levelSet );

% R = ( NalevelSet - aNlevelSet ) ./ NalevelSet;
% Rxi = ( aNlevelSet .* (Nxi * alevelSet) - ...
%              NalevelSet .* s .* (Nxi*levelSet) ) ./ NalevelSet.^2;
% Reta = ( aNlevelSet .* (Neta * alevelSet) - ...
%              NalevelSet .* s .* (Neta*levelSet) ) ./ NalevelSet.^2;
%
R = N*alevelSet - abs( N*levelSet );
Rxi = Nxi*alevelSet - s.*(Nxi*levelSet) ;
Reta = Neta*alevelSet - s.*(Neta*levelSet) ;
