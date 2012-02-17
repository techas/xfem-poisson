function [ eL2, eH1, eH121,ecomp1,eH122,ecomp2,exL2,exH1 ] = ComputeEnergy( X, T, h, hE, levelSet, r0, f, opts )
% computes the norm of the error using different norms

% global
global useEnrichment
global cond
nu1 = cond(2); nu2 = cond(1);
condsav = cond;

% constants
numberOfElements = size(T,1);
numberOfGauss = 16;
[pospg,pespg] = quadrature( opts.elementType, numberOfGauss );
[N,Nxi,Neta] = shapeFunction( opts.elementType, ...
                                       opts.numberOfElementNodes, pospg );
[pospg1D,pespg1D] = quadrature1D;
Nxig = Nxi(1:5,:);
Netag = Neta(1:5,:);

% find enriched elements
[type,enrichedNodes] = classifyElements( levelSet, T, opts.tolerance );
standardElements = find( type==0 );
enrichedElements = find( type>0 );

% value of the approximate solution at gauss points
[ Xpg, qH ] = Flujos( h, hE, X, T, pospg, N, Nxi, Neta, levelSet, opts );
rpg = sqrt( Xpg(:,1).^2+Xpg(:,2).^2 );
hH = zeros( size(rpg) );
for i1 = 1:size(T,1)
    Te = T(i1,:);
    Xe = X( Te, : );
    LSe = levelSet( Te );
    ind = inpolygon( Xpg(:,1), Xpg(:,2), Xe(:,1), Xe(:,2) );
    hH( ind ) = N * h( Te );
    if (useEnrichment && type(i1)>0)
        R = buildRidge( LSe, N, Nxi, Neta );
        Re = vectorFind( enrichedNodes, Te );
        hH( ind ) = hH( ind )+R.*(N*hE( Re ));
    end
    indT = max(find(ind,1,'first'))-1;
    ind = indT+find( (N*LSe)>0 );
    qH( ind, : ) = qH( ind, : )/nu2;
    ind = indT+find( (N*LSe)<0 );
    qH( ind, : ) = qH( ind, : )/nu1;
end

% value of the approximate flux along the interface
Neg = length(enrichedElements);
qHg1 = zeros( 5, 2, Neg );
qHg2 = zeros( 5, 2, Neg );
Xg = zeros( 5, 2, Neg );
Lg = zeros( 1, Neg );
for i1 = 1:Neg
    Te = T(enrichedElements(i1),:);
    Xe = X( Te, : );
    LSe = levelSet( Te );
    [Xp1,Xp2,Rp1,Rp2,Ng1,Ng2] = FindIntersectionNodes( Xe, LSe );
    Lg(i1) = norm(Xp2-Xp1);
    Ng = (pospg1D*(Ng2-Ng1)+ones(size(pespg1D))*(Ng2+Ng1))/2;
    cond = [1 1];
    [Xg(:,:,i1),qHg1(:,:,i1)]=...
          Flujos( h, hE, X, Te, pospg1D, Ng, Nxig, Netag, levelSet, ...
                                              opts, -1*ones(size(pospg1D)) );
    [Xg(:,:,i1),qHg2(:,:,i1)]=...
          Flujos( h, hE, X, Te, pospg1D, Ng, Nxig, Netag, levelSet, ...
                                               opts, ones(size(pospg1D)) );
end
cond = condsav;

% value of the exact solution at gauss points
% eex = pi*f^2/8*(r0^4/cond(2)+(1-r0^4)/cond(1));
epsnu = 1-nu2/nu1;
u0 = -f*r0^2/4/nu2*(1/r0^2-epsnu);
hex = zeros( size(rpg) );
qex = zeros( size(Xpg) );
ind = rpg <= r0; 
hex(ind) = u0 + rpg(ind).^2*f/4/nu1;
qex(ind,:) = -Xpg(ind,:)*f/2/nu1;
ind = rpg > r0; 
hex(ind) = u0 + r0^2*f/4/nu2*(rpg(ind).^2/r0^2-epsnu);
qex(ind,:) = -Xpg(ind,:)*f/2/nu2;

% value of exact fluxes along the interface
qexg1 = -Xg*f/2/nu1;
qexg2 = -Xg*f/2/nu2;

% errors at gauss points
dq = qH-qex;
dh = hH-hex;
dqg1 = qHg1-qexg1;
dqg1n = squeeze( (dqg1(:,1,:).*Xg(:,1,:)+dqg1(:,2,:).*Xg(:,2,:)) ./ ...
        sqrt(Xg(:,1,:).^2+Xg(:,2,:).^2) ).^2;
dqg1t = squeeze( dqg1(:,1,:).^2 + dqg1(:,2,:).^2 );
dqg2 = qHg2-qexg2;
dqg2n = squeeze( (dqg2(:,1,:).*Xg(:,1,:)+dqg2(:,2,:).*Xg(:,2,:)) ./ ...
        sqrt(Xg(:,1,:).^2+Xg(:,2,:).^2) ).^2;
dqg2t = squeeze( dqg2(:,1,:).^2 + dqg2(:,2,:).^2 );

% integrating over the standard elements
eH1 = zeros( numberOfElements, 1 );
exH1 = zeros( numberOfElements, 1 );
eL2 = zeros( numberOfElements, 1 );
exL2 = zeros( numberOfElements, 1 );
for I = 1:length(standardElements)
   Te = T(standardElements(I),:); 
   Xe = X(Te,:);
   ind = inpolygon( Xpg(:,1), Xpg(:,2), Xe(:,1), Xe(:,2) );
   dqe = dq(ind,:);
   dhe = dh(ind);
   qexe = qex(ind,:);
   hexe = hex(ind);
   eH1(standardElements(I)) = ...
             elementMatrixStandardNoJacob( Xe, pespg, dqe(:,1), dqe(:,2) );
   exH1(standardElements(I)) = ...
           elementMatrixStandardNoJacob( Xe, pespg, qexe(:,1), qexe(:,2) );
   eL2(standardElements(I)) = ...
          elementMatrixStandardNoJacob( Xe, pespg, dhe, zeros(size(dhe)) );
   exL2(standardElements(I)) = ...
         elementMatrixStandardNoJacob( Xe, pespg, hexe, zeros(size(dhe)) );
end 

% fill the matrix for the enriched elements
for I = 1:length( enrichedElements )
   Te = T(enrichedElements(I),:);
   Xe = X(Te,:); 
   LSe = levelSet( Te );
   ind = inpolygon( Xpg(:,1), Xpg(:,2), Xe(:,1), Xe(:,2) );
   dqe = dq(ind,:);
   dhe = dh(ind);
   eH1(enrichedElements(I)) = ...
        elementMatrixEnrichedNoJacob( Xe, LSe, pespg, dqe(:,1), dqe(:,2) );
   exH1(enrichedElements(I)) = ...
      elementMatrixEnrichedNoJacob( Xe, LSe, pespg, qexe(:,1), qexe(:,2) );
   eL2(enrichedElements(I)) = ...
     elementMatrixEnrichedNoJacob( Xe, LSe, pespg, dhe, zeros(size(dhe)) );
   exL2(enrichedElements(I)) = ...
    elementMatrixEnrichedNoJacob( Xe, LSe, pespg, hexe, zeros(size(dhe)) );
end

% integration for elements along the interface
eH121 =  (dqg1n'*pespg1D).*Lg';
ecomp1 =  (dqg1t'*pespg1D).*Lg';
eH122 =  (dqg2n'*pespg1D).*Lg';
ecomp2 =  (dqg2t'*pespg1D).*Lg';

% function quadrature1D
function [pospg,pespg] = quadrature1D
pospg = [-sqrt(5+2*sqrt(10/7))/3
         -sqrt(5-2*sqrt(10/7))/3
          0
          sqrt(5-2*sqrt(10/7))/3
          sqrt(5+2*sqrt(10/7))/3];
pespg = [ (322-13*sqrt(70))/900
          (322+13*sqrt(70))/900
          128/225
          (322+13*sqrt(70))/900
          (322-13*sqrt(70))/900];
% function  findintersectionnodes
function [Xp1,Xp2,Rp1,Rp2,Ng1,Ng2] = FindIntersectionNodes( Xe, LSe )
negNodes = find( LSe < 0 );
switch length( negNodes )
    case 1
        nodosolo = negNodes;
        nodopar = find( LSe >= 0 );
    case 2
        nodosolo = find( LSe >= 0 );
        nodopar = negNodes;
    otherwise
        error( 'elementMatrixEnriched' )
end
X1 = [Xe(nodosolo,:); Xe(nodopar(1),:)];
X2 = [Xe(nodosolo,:); Xe(nodopar(2),:)];
LS1 = [LSe(nodosolo) LSe(nodopar(1))]; aLS1 = abs(LS1);
LS2 = [LSe(nodosolo) LSe(nodopar(2))]; aLS2 = abs(LS2);
Xp1 = intersection( X1, LS1 );
Xp2 = intersection( X2, LS1 );
ind1 = ( aLS1==max(aLS1) );
alpha1 = norm(Xp1-X1(ind1,:))./norm(X1(~ind1,:)-X1(ind1,:));
Rp1 = aLS1(ind1)-alpha1*(aLS1(ind1)-aLS1(~ind1));
ind2 = ( aLS2==max(aLS2) );
alpha2 = norm(Xp2-X2(ind2,:))./norm(X2(~ind2,:)-X2(ind2,:));
Rp2 = aLS2(ind2)-alpha2*(aLS2(ind2)-aLS2(~ind2));
A = [ Xe ones(3,1) ]\[1 0;0 1;0 0];
xp1 = [Xp1 1] * A;
xp2 = [Xp2 1] * A;
Ng1 = [xp1(1) xp1(2) 1-xp1(1)-xp1(2)];
Ng2 = [xp2(1) xp2(2) 1-xp2(1)-xp2(2)];

        