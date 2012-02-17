function [F] = VectorF(X,T,pospg,pespg,N,Nxi,Neta,Fpg,levelSet,opts)

%numero de elementos y numero de nodos en cada elemento
[nelem,nnode] = size(T);
%numero de nodos total
npoin = size(X,1);
ngaus = size(pospg,1);

%Dimensionamiento
F = sparse(npoin,1);

% Look for the enriched elements and nodes
[type,enrichedNodes] = classifyElements( levelSet, T, opts.tolerance );
standardElements = find( type == 0 );
enrichedElements = find( type > 0 );
nstand = length(standardElements);

% Loop on standard elements
for ielem = 1:length( standardElements )
    Te = T(standardElements(ielem),:);
    Xe = X(Te,:);
    Fpge = Fpg((ielem-1)*ngaus+[1:ngaus],:);
    [Fe] = VFel(Xe,nnode,pospg,pespg,N,Nxi,Neta,Fpge);
    F(Te)=F(Te)+Fe;
end
% Loop on enriched elements
for ielem = 1:length( enrichedElements )
    Te = T(enrichedElements(ielem),:);
    Xe = X(Te,:);
    Fpge = Fpg(nstand+(ielem-1)*ngaus+[1:ngaus],:);
    [Fe] = VFel(Xe,nnode,pospg,pespg,N,Nxi,Neta,Fpge);
    F(Te)=F(Te)+Fe;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fe] = VFel(Xe,nnode,pospg,pespg,N,Nxi,Neta,Fpge)
%Inicializacion
Fe = zeros(nnode,1);
%numero de puntos de Gauss
ngaus = size(pospg,1);

% BUCLE EN PUNTOS DE GAUSS
for igaus = 1:ngaus
    %Funciones de forma y derivadas en el punto de Gauss (de los nnode nodos)
    N_igaus = N(igaus,:);
    Nxi_igaus = Nxi(igaus,:);
    Neta_igaus = Neta(igaus,:);
    %Jacobiano
    Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))
             Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
    %Peso gauss correspondente
    dvolu=pespg(igaus)*det(Jacob);

    %Contribucion del punto de Gauss al vector
    Fe  = Fe  - (N_igaus'*Fpge(igaus))*dvolu;
end
