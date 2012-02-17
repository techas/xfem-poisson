function [Ke,fe] = elementMatrixStandard( Xe, numberOfNodes, pospg, ...
                  pespg, N, Nxi, Neta, material, caseLoad ) 
% [Ke,fe] = elementMatrixStandard( Xe, numberOfNodes, pospg, pespg, ...
%   N, Nxi, Neta ) 
%
% creates an elemental matrix
%
% INPUT
%   Xe             nodal coords
%   numberOfNodes  number of element nodes
%   pospg,pespg   position and weigth of gauss points 
%   N,Nxi,Neta    shape functions and its derivetives
%
% OUTPUT
%   K
%   f
%
global cond

numberOfGaussPoints = size( pospg, 1 ); 

Ke = zeros( numberOfNodes, numberOfNodes ); 
fe = zeros( numberOfNodes, 1 ); 

Xpg = Isopar( Xe, N );

for igaus = 1:numberOfGaussPoints
    Xpg1 = Xpg( igaus, : );
    dN = [ Nxi(igaus,:) ; Neta(igaus,:) ];
    jacob = dN*Xe ;
    dvolu = pespg(igaus) * det( jacob );
    res = jacob \ dN;
    Nx = res(1,:);
    Ny = res(2,:);
    Ke = Ke + (Nx'*Nx + Ny'*Ny)*dvolu*cond(material);
    if caseLoad~=2
        fe = fe + N(igaus,:)' * (sourceTerm( Xpg1, caseLoad )*dvolu );
    end
end

