function DibujaVec(X,Vectores, str, fact0 )
% DibujaVec(X,Vectores,str) dibuja un campo vectorial
% X:        coordenadas de los puntos base
% Vectores: vector en cada punto base
% str     : color para el dibujo (opcional)

if nargin > 2
    simb = str;
else
    simb = 'k-';
end
if nargin > 3
    fact = fact0;
else
    fact = 0.1;
end

Vmax = max (sqrt(Vectores(:,1).^2+Vectores(:,2).^2));
Dmax=min(max(X(:,1))-min(X(:,1)),max(X(:,2))-min(X(:,2)));
scalefac = Dmax/Vmax*fact;

npunts = size(X,1);

seno = sin(pi/8);
cose = cos(pi/8);

for i=1:npunts
   x = X(i,1);
   y = X(i,2);
   u = Vectores(i,1)*scalefac;
   v = Vectores(i,2)*scalefac;
   hu1 = (x+u) - (cose*u+seno*v)/6; 
   hv1 = (y+v) - (-seno*u+cose*v)/6; 
   hu2 = (x+u) - (cose*u-seno*v)/6; 
   hv2 = (y+v) - (seno*u+cose*v)/6; 
   hold on
   plot([x,x+u],[y,y+v],simb,[hu1,x+u,hu2],[hv1,y+v,hv2],simb,'LineWidth',2)
   hold off
end 
axis('equal');
%axis('off');

