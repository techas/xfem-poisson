function [M] = MatrixM(X,T,pospg,pespg,N,Nxi,Neta) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [K,f] = Cremat(X,T,pospg,pespg,N,Nxi,Neta) 
% Matriz de rigidez K y término independiente f (fuente) 
% 
% X:            coordenadas nodales 
% T:            conectividades (elementos) 
% pospg, pespg: posición y pesos de los puntos de Gauss en el 
%               elemento de referencia 
% N,Nxi,Neta:   funciones de forma y derivadas con coordenadas locales 
%               en los puntos de Gauss del elemento de referencia 
% w:            frecuencia
% g:            gravidad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de elementos y numero de nodos en cada elemento 
[nelem,nnode] = size(T); 
%numero de nodos total  
npoin = size(X,1); 
 
% Dimensionamiento 
M = sparse(npoin,npoin); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUCLE EN ELEMENTOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ielem = 1:nelem 
 % Xe: coordenadas de los nodos del elemento 
 Te = T(ielem,:); 
 Xe = X(Te,:); 

 % Matriz y vector elementales
 [Me] = MMel(Xe,nnode,pospg,pespg,N,Nxi,Neta); 

 M(Te,Te)=M(Te,Te)+Me; 
 clear Me;
end 
 
 
 
