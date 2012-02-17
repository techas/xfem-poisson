function [F] = VectorF(X,T,pospg,pespg,N,Nxi,vn)

%numero de elementos y numero de nodos en cada elemento
[nelem,nnode] = size(T);
%numero de nodos total 
npoin = size(X,1);

%Dimensionamiento
F = sparse(npoin,1); 

% Loop in elements
for ielem = 1:nelem
 % Xe: coordinates of nodes of element
 Te = T(ielem,:);
 Xe = X(Te,:);
 [Fe] = VFel(Xe,nnode,pospg,pespg,N,Nxi,vn);
 F(Te)=F(Te)+Fe; 
 clear  Fe;
end

