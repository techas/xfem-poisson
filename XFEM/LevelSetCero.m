%  
% [xi,eta] = LevelSetCero([1 -1 -1])
%
function ...
    [xi,eta,Arista] = LevelSetCero(levelset)
    f1 = levelset(3); f2 = levelset(1); f3 = levelset(2);
    
    caso1 = (f2*f3 >= 0 && f1~=0);
    caso2 = (f1*f3 >= 0);
    caso3 = (f1*f2 >= 0);
    
    if caso1 == 1
         xi = [f1/(f1-f2); 0]; Arista1 = [1 3];
        eta = [0; f1/(f1-f3)]; Arista2 = [2 3];
    elseif caso2 == 1
         xi = [f1/(f1-f2); f3/(f3-f2)]; Arista1 = [1 3];
        eta = [         0; f2/(f2-f3)]; Arista2 = [1 2];      
    elseif caso3 ==1
         xi = [         0; f3/(f3-f2)]; Arista1 = [2 3];
        eta = [f1/(f1-f3); f2/(f2-f3)]; Arista2 = [1 2];      
    else
        error('Caso no identificado')
    end
    Arista = [Arista1; Arista2];