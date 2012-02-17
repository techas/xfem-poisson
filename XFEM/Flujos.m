function [Xpg,Fpg] = Flujos(Temp,TempE,X,T,pospg,...
    N,Nxi,Neta,levelSet,opts,force)
% [Xpg,Fpg] = Flujos(Temp,X,T,pospg,pespg,N,Nxi,Neta)
% Calcula el flujo asociado a un campo de temperaturas Temp en
% los puntos de Gauss
% Xpg:  coordenadas de los puntos de Gauss (puntos base del campo resultado)
% Fpg:  vector de flujo en cada punto de Gauss
% Temp: vector con las temperaturas en cada nodo
% X:            coordenadas nodales
% T:            conectividades (elementos)
% pospg: posicion los puntos de Gauss en el
%        elemento de referencia
% N,Nxi,Neta:   funciones de forma y derivadas con coordenadas locales
%               en los puntos de Gauss del elemento de referencia
%
% La instruccion DibujaVec(xpg,Fpg) dibuja el flujo
%
global cond useEnrichment

nelem  = size(T,1);
ngaus = size(pospg,1);
npg = nelem*ngaus;

Xpg = zeros(npg,2);
Fpg = zeros(npg,2);

if nargin<11
    force = NaN(ngaus,1);
end

%numero de puntos de Gauss
%ngaus = size(pospg,1);
% Look for the enriched elements and nodes
[type,enrichedNodes] = classifyElements( levelSet, T, opts.tolerance );
standardElements = find( type == 0 );
enrichedElements = find( type > 0 );

%contador para el numero global de punto de Gauss
ipg = 1;

for ielem = 1:length( standardElements )
    Te = T(standardElements(ielem),:);
    Xe = X(Te,:);
    LSe = levelSet(Te);
    material = mean(materialpg( LSe((abs(LSe)==max(abs(LSe)))), 0 ));
%    material = materialpg( LSe(1), 0 );
    %bucle en puntos de Gauss
    for igaus = 1:ngaus
        N_igaus = N(igaus,:);
        Nxi_igaus = Nxi(igaus,:);
        Neta_igaus = Neta(igaus,:);
        Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))
            Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
        T_xi  = Nxi_igaus*Temp(Te);
        T_eta = Neta_igaus*Temp(Te);
        res = Jacob\[T_xi;T_eta];
        Temp_x = res(1);
        Temp_y = res(2);
        Fpg(ipg,:) = -cond(material)*[Temp_x,Temp_y];
        Xpg(ipg,:) = Isopar(Xe,N_igaus);
        ipg = ipg + 1;
    end
end

for ielem = 1:length( enrichedElements )
    Te = T(enrichedElements(ielem),:);
    Xe = X(Te,:);
    LSe = levelSet(Te);
    if useEnrichment
        Re = vectorFind( enrichedNodes, Te );
        %      Te = [Te, Re+numberOfNodes];
    end

    %bucle en puntos de Gauss
    for igaus = 1:ngaus
        %
        N_igaus = N(igaus,:);
        LSepg = N_igaus *  LSe ;
        material = materialpg( LSepg, 0 );
        %
        Nxi_igaus = Nxi(igaus,:);
        Neta_igaus = Neta(igaus,:);
        Jacob = [Nxi_igaus*(Xe(:,1))	Nxi_igaus*(Xe(:,2))
            Neta_igaus*(Xe(:,1))	Neta_igaus*(Xe(:,2))];
        T_xi  = Nxi_igaus*Temp(Te);
        T_eta = Neta_igaus*Temp(Te);
        f_gaus = force(igaus);

        if useEnrichment
            [R,Rxi,Reta] = buildRidge( LSe, N_igaus, Nxi_igaus, Neta_igaus, f_gaus );
%            M = N_igaus.* R; %repmat(R, 1, size( N_igaus, 2 ) );
            Mxi  = Nxi_igaus  * R + N_igaus * Rxi ;
            Meta = Neta_igaus * R + N_igaus * Reta;
            T_xi   = T_xi  + Mxi  * TempE(Re);
            T_eta  = T_eta + Meta * TempE(Re);
        end
        %
        res = Jacob\[T_xi;T_eta];
        Temp_x = res(1);
        Temp_y = res(2);
        Fpg(ipg,:) = -cond(material)*[Temp_x,Temp_y];
        Xpg(ipg,:) = Isopar(Xe,N_igaus);
        ipg = ipg + 1;
    end
end

