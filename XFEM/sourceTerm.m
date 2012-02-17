function f = sourceTerm( x , caseLoad) 
%  f = sourceTerm( x, caseLoad) 
%
%  x: Gauss point coordinates
%  caseLoad: type of loading case
%            1: uniform source term f=1
%            2: uniform source term f=0
%            3: source term f=1 within a circle at (.5 .5) with a radius
%               .1, and f=0 elsewhere.
%

switch caseLoad
    case 1
        f = 1;
    case 2
        f = 0;
    case 3
        if norm(x-[.5 .5])<0.1
            f = 1;
        else
            f = 0;
        end
    case 4
        f = -4;
end
