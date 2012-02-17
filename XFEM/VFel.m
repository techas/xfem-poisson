function [Fe] = VFel(Xe,nnode,pospg,pespg,N,Nxi,vn)

Fe = zeros(nnode,1); 

% Numero de puntos de Gauss
ngaus = size(pospg,1);

x1 = Xe(1,:); x2= Xe(nnode,:);
Ax = norm(x2-x1);

%normal
n1 =[-(x2(2)-x1(2))   ;  (x2(1)-x1(1)) ]/Ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUCLE EN PUNTOS DE GAUSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for igaus = 1:ngaus 
  
  dvolu=pespg(igaus)*Ax/2;
  N_igaus = N(igaus,:);
  
  %Contribucion del punto de Gauss al vector
  Fe = Fe + (N_igaus'*vn)*dvolu;

end
