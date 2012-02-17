function [Ke,fe] = elementMatrixEnrichedPart( Xe, Xp, levelSet, ...
   numberOfNodes, pospg, pespg, N, Nxi, Neta, material, caseLoad ) 
% [Ke,fe] = elementMatrixEnrichedPart( Xe, Xp, numberOfNodes, pospg, pespg, ...
%    N, Nxi, Neta, material )  
%
% integrates a single material part of an element
%
% INPUT
%   Xe             nodal coords of the element
%   Xp             nodal coords of the part to integrate
%   numberOfNodes  number of element nodes
%   pospg,pespg    position and weigth of gauss points 
%   N,Nxi,Neta     shape functions and its derivetives
%
% OUTPUT
%   K
%   f
%
global cond useEnrichment 

numberOfGaussPoints = size( pospg, 1 ); 
if useEnrichment
   ksize = 2*numberOfNodes;
else
   ksize = numberOfNodes;
end

Ke = zeros( ksize, ksize ); 
fe = zeros( ksize, 1 ); 

% Integration points and shape function for the small element
theta = N*Xp;
[xi,eta] = invMap( Xe, theta(:,1), theta(:,2) );
pgs = [xi, eta];
triangle = 2;
nPoints = 3;
[Ns,Ns_xi,Ns_eta] = shapeFunction( triangle, nPoints, pgs );

%plotDebug( Xe, Xp, Ns );

% Ridge function
if useEnrichment
   [Rs,Rsxi,Rseta] = buildRidge( levelSet, Ns, Ns_xi, Ns_eta );
   % AAAAAAAHHHHHHHHHHHH!!!!
%    M = N .* repmat(Rs, 1, size( N, 2 ) );
   Ms = Ns .* repmat(Rs, 1, size( N, 2 ) );
end

for igaus = 1:numberOfGaussPoints
    dN = [ Nxi(igaus,:) ; Neta(igaus,:)];
    dNs = [ Ns_xi(igaus,:) ; Ns_eta(igaus,:)];
    %
    jacob = dNs * Xe;
    jacobS = dN * Xp;
    dvolu = pespg(igaus) * det( jacobS );
    res = jacob \ dNs;
    Nx = res(1,:);
    Ny = res(2,:);
    if useEnrichment
        resR = jacob\[Rsxi(igaus); Rseta(igaus)];
        Mx = Nx * Rs(igaus) + Ns(igaus,:) * resR(1);
        My = Ny * Rs(igaus) + Ns(igaus,:) * resR(2);
        Tx = [Nx Mx];
        Ty = [Ny My];
        if caseLoad~=2
            tmpNXp = N(igaus,:)*Xp;
            fe = fe + [Ns(igaus,:) Ms(igaus,:)]' * ...
                (sourceTerm( tmpNXp , caseLoad)*dvolu);
        end
    else
        Tx = Nx;
        Ty = Ny;
        if caseLoad~=2
            tmpNsXe = Ns(igaus,:)*Xe;
            fe = fe + Ns(igaus,:)' * ...
                (sourceTerm( tmpNsXe, caseLoad)*dvolu);
        end
    end
    Ke = Ke + (Tx'*Tx + Ty'*Ty)*dvolu*cond(material);
end

