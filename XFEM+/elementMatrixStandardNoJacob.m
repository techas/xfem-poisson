function Ke = elementMatrixStandardNoJacob( Xe, pespg, Nx, Ny ) 

J = [1 0 -1;0 1 -1]*Xe;
dvolu = pespg*abs(det(J));
Ke = (Nx'.*dvolu)*Nx + (Ny'.*dvolu)*Ny;
