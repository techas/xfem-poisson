function gradLevelSet = GradientLevelSet(levelSet,X,T,pospg,...
                                                           N,Nxi,Neta,opts)
% GRADIENTLEVELSET to evaluate the gradient of the level set function at
% the Gauss points
%
%

nelem  = size(T,1);
ngaus = size(pospg,1);
npg = nelem*ngaus;

gradLevelSet = zeros(npg,2);

%numero de puntos de Gauss
ngaus = size(pospg,1);
% Look for the enriched elements and nodes
[type,enrichedNodes] = classifyElements( levelSet, T, opts.tolerance );
standardElements = find( type == 0 );
enrichedElements = find( type > 0 );

%contador para el numero global de punto de Gauss
ipg = 1;

for ielem = 1:length( standardElements )
    Te = T(standardElements(ielem),:);
    Xe = X(Te,:);
    %bucle en puntos de Gauss
    for igaus = 1:ngaus
%        N_igaus = N(igaus,:);
        Nxi_igaus = Nxi(igaus,:);
        Neta_igaus = Neta(igaus,:);
        Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))
            Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
        T_xi  = Nxi_igaus*levelSet(Te);
        T_eta = Neta_igaus*levelSet(Te);
        gradLevelSet(ipg,:) = Jacob\[T_xi;T_eta];
        ipg = ipg + 1;
    end
end
for ielem = 1:length( enrichedElements )
    Te = T(enrichedElements(ielem),:);
    Xe = X(Te,:);
    %bucle en puntos de Gauss
    for igaus = 1:ngaus
        %
%        N_igaus = N(igaus,:);
        %
        Nxi_igaus = Nxi(igaus,:);
        Neta_igaus = Neta(igaus,:);
        Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))
            Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
        T_xi  = Nxi_igaus*levelSet(Te);
        T_eta = Neta_igaus*levelSet(Te);
        %
        gradLevelSet(ipg,:) = Jacob\[T_xi;T_eta];
        ipg = ipg + 1;
    end
end

