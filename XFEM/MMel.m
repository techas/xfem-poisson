function [Me] = MMel(Xe,nnode,pospg,pespg,N,Nxi,Neta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matriz de Masa elemental
%
% Xe:           coordenadas de los nodos del elemento
% nnode:        numero de nodos del elemento
% pospg, pespg: posicion y pesos de los puntos de Gauss en el
%               elemento de referencia
% N,Nxi,Neta:   funciones de forma y derivadas (con coordenadas locales)
%               en los puntos de Gauss del elemento de referencia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inicializacion
Me = zeros(nnode,nnode);

%numero de puntos de Gauss
ngaus = size(pospg,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUCLE EN PUNTOS DE GAUSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %Contribucion del punto de Gauss a la matriz
    Me  = Me  + (N_igaus'*N_igaus)*dvolu;
end




